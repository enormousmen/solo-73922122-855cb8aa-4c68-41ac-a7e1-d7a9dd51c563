"""SQLite metadata layer for the versioned VCF store.

All multi-step state changes happen inside BEGIN IMMEDIATE transactions so
that two app instances sharing the same database file are fully serialized.
This is what makes "only one concurrent replacement succeeds" work.
"""
from __future__ import annotations

import os
import sqlite3
import time
from contextlib import contextmanager

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (
    k TEXT PRIMARY KEY,
    v TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS versions (
    version   INTEGER PRIMARY KEY,
    chrom     TEXT NOT NULL,
    path      TEXT NOT NULL,      -- directory with variants.vcf.gz + .tbi, relative to data dir
    n_records INTEGER NOT NULL,
    created_at REAL NOT NULL
);
CREATE TABLE IF NOT EXISTS refs (
    id        INTEGER PRIMARY KEY AUTOINCREMENT,
    version   INTEGER NOT NULL,
    instance  TEXT NOT NULL,
    created_at REAL NOT NULL
);
CREATE TABLE IF NOT EXISTS instances (
    id        TEXT PRIMARY KEY,
    last_seen REAL NOT NULL
);
"""


class Conflict(Exception):
    """Raised when base_version does not match the current version."""

    def __init__(self, current: int):
        super().__init__(f"current version is {current}")
        self.current = current


class Store:
    def __init__(self, db_path: str):
        self.db_path = db_path
        parent = os.path.dirname(db_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        con = sqlite3.connect(db_path, timeout=30.0)
        try:
            # WAL: readers never block the single writer and vice versa.
            con.execute("PRAGMA journal_mode=WAL")
            con.execute("PRAGMA synchronous=NORMAL")
            con.executescript(SCHEMA)
            con.execute("INSERT OR IGNORE INTO meta(k, v) VALUES ('current_version', '0')")
            con.commit()
        finally:
            con.close()

    def _connect(self) -> sqlite3.Connection:
        con = sqlite3.connect(self.db_path, timeout=30.0)
        con.execute("PRAGMA busy_timeout=30000")
        con.row_factory = sqlite3.Row
        return con

    @contextmanager
    def _txn(self):
        """BEGIN IMMEDIATE transaction: serialized against all other writers."""
        con = self._connect()
        try:
            con.execute("BEGIN IMMEDIATE")
            yield con
            con.commit()
        except BaseException:
            con.rollback()
            raise
        finally:
            con.close()

    # ------------------------------------------------------------------ reads

    def current_version(self) -> int:
        con = self._connect()
        try:
            row = con.execute("SELECT v FROM meta WHERE k = 'current_version'").fetchone()
            return int(row["v"])
        finally:
            con.close()

    def version_info(self):
        """Row of the currently published version, or None if nothing published."""
        con = self._connect()
        try:
            return con.execute(
                "SELECT v.version, v.chrom, v.path, v.n_records, v.created_at "
                "FROM meta m JOIN versions v ON v.version = CAST(m.v AS INTEGER) "
                "WHERE m.k = 'current_version'"
            ).fetchone()
        finally:
            con.close()

    def all_version_paths(self) -> set:
        con = self._connect()
        try:
            return {r["path"] for r in con.execute("SELECT path FROM versions")}
        finally:
            con.close()

    def live_instances(self, now: float, dead_after: float) -> set:
        con = self._connect()
        try:
            cutoff = now - dead_after
            return {
                r["id"]
                for r in con.execute("SELECT id FROM instances WHERE last_seen >= ?", (cutoff,))
            }
        finally:
            con.close()

    # ------------------------------------------------------------------ publish

    def publish(self, base_version: int, chrom: str, path: str, n_records: int,
                now: float | None = None) -> int:
        """Atomically make `path` the current version iff current == base_version."""
        now = time.time() if now is None else now
        with self._txn() as con:
            cur = int(con.execute("SELECT v FROM meta WHERE k = 'current_version'").fetchone()["v"])
            if cur != base_version:
                raise Conflict(cur)
            new_version = cur + 1
            con.execute(
                "INSERT INTO versions(version, chrom, path, n_records, created_at) "
                "VALUES (?, ?, ?, ?, ?)",
                (new_version, chrom, path, n_records, now),
            )
            con.execute("UPDATE meta SET v = ? WHERE k = 'current_version'", (str(new_version),))
            return new_version

    # ------------------------------------------------------------------ refs

    def pin(self, instance_id: str, now: float | None = None):
        """Pin the current version for an in-flight query.

        Returns (version, chrom, path, ref_id) or None when nothing is published.
        The read of the current pointer and the ref insert are one transaction,
        so the reaper can never delete a version between the two steps.
        """
        now = time.time() if now is None else now
        with self._txn() as con:
            cur = int(con.execute("SELECT v FROM meta WHERE k = 'current_version'").fetchone()["v"])
            if cur == 0:
                return None
            row = con.execute(
                "SELECT chrom, path FROM versions WHERE version = ?", (cur,)
            ).fetchone()
            if row is None:  # pragma: no cover - defensive, should not happen
                return None
            c = con.execute(
                "INSERT INTO refs(version, instance, created_at) VALUES (?, ?, ?)",
                (cur, instance_id, now),
            )
            return cur, row["chrom"], row["path"], c.lastrowid

    def unpin(self, ref_id: int) -> None:
        with self._txn() as con:
            con.execute("DELETE FROM refs WHERE id = ?", (ref_id,))

    # ------------------------------------------------------------------ heartbeat

    def heartbeat(self, instance_id: str, now: float | None = None) -> None:
        now = time.time() if now is None else now
        with self._txn() as con:
            con.execute(
                "INSERT INTO instances(id, last_seen) VALUES (?, ?) "
                "ON CONFLICT(id) DO UPDATE SET last_seen = excluded.last_seen",
                (instance_id, now),
            )

    def unregister(self, instance_id: str) -> None:
        with self._txn() as con:
            con.execute("DELETE FROM refs WHERE instance = ?", (instance_id,))
            con.execute("DELETE FROM instances WHERE id = ?", (instance_id,))

    # ------------------------------------------------------------------ reaper

    def reap_metadata(self, now: float, dead_after: float):
        """Drop refs of dead instances and delete rows of reclaimable versions.

        A version is reclaimable when it is not the current one and no live
        ref points at it. Returns [(version, path)] whose files the caller
        should delete (after commit, so a crash leaves only sweepable orphans).
        """
        with self._txn() as con:
            cutoff = now - dead_after
            con.execute(
                "DELETE FROM refs WHERE instance NOT IN "
                "(SELECT id FROM instances WHERE last_seen >= ?)",
                (cutoff,),
            )
            cur = int(con.execute("SELECT v FROM meta WHERE k = 'current_version'").fetchone()["v"])
            rows = con.execute(
                "SELECT version, path FROM versions WHERE version != ? AND NOT EXISTS "
                "(SELECT 1 FROM refs WHERE refs.version = versions.version)",
                (cur,),
            ).fetchall()
            con.executemany(
                "DELETE FROM versions WHERE version = ?", [(r["version"],) for r in rows]
            )
            return [(r["version"], r["path"]) for r in rows]

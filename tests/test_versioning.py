"""Versioning, optimistic concurrency, old-query protection, reclamation."""
import os
import sqlite3
import threading
import time

import pysam
from conftest import SAMPLE, make_config, upload, vcf

from app.db import Conflict, Store
from app.main import create_app
from app.vcfbuild import GZ_NAME


def _version_dir(app, version_path):
    return os.path.join(app.state.cfg.data_dir, version_path)


def test_concurrent_publish_only_one_wins(tmp_path):
    """Two writers racing with the same base_version: exactly one succeeds."""
    store = Store(os.path.join(str(tmp_path), "meta.sqlite3"))
    store.publish(0, "chr1", "versions/v1", 10)

    results = {"ok": 0, "conflict": 0}
    lock = threading.Lock()
    barrier = threading.Barrier(8)

    def race(i):
        barrier.wait()
        try:
            store.publish(1, "chr1", f"versions/racer{i}", 1)
            with lock:
                results["ok"] += 1
        except Conflict:
            with lock:
                results["conflict"] += 1

    threads = [threading.Thread(target=race, args=(i,)) for i in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert results == {"ok": 1, "conflict": 7}
    assert store.current_version() == 2


def test_old_query_reads_old_version_then_files_reclaimed(client, app):
    upload(client, SAMPLE, base_version=0)
    store = app.state.store
    v1_path = store.version_info()["path"]

    # simulate an in-flight query on version 1
    pinned = store.pin("test-inst")
    assert pinned[0] == 1
    ref_id = pinned[3]

    # replacement gets published while the old query is still open
    r = upload(client, vcf([(500, "A", "G")]), base_version=1)
    assert r.status_code == 201 and r.json()["version"] == 2

    # reaper must not touch the pinned old version
    app.state.reap_once()
    assert os.path.isdir(_version_dir(app, v1_path))

    # the in-flight query still reads version 1 data through its index
    tb = pysam.TabixFile(os.path.join(_version_dir(app, v1_path), GZ_NAME))
    try:
        poss = [ln.split("\t")[1] for ln in tb.fetch("chr1", 99, 106)]
    finally:
        tb.close()
    assert poss == ["100", "105"]

    # new queries see version 2
    q = client.get("/variants", params={"chrom": "chr1", "start": 1, "end": 1000})
    assert q.json()["version"] == 2
    assert [rec.split("\t")[1] for rec in q.json()["records"]] == ["500"]

    # old query finishes -> reaper removes version 1 files and row
    store.unpin(ref_id)
    app.state.reap_once()
    assert not os.path.exists(_version_dir(app, v1_path))
    assert store.current_version() == 2
    q = client.get("/variants", params={"chrom": "chr1", "start": 1, "end": 1000})
    assert q.json()["version"] == 2


def test_dead_instances_pins_do_not_block_reclaim(client, app):
    upload(client, SAMPLE, base_version=0)
    store = app.state.store
    v1_path = store.version_info()["path"]

    # a ref owned by an instance that died without unpinning
    now = time.time()
    store.heartbeat("ghost", now - 1000)  # long past dead_after
    with store._txn() as con:
        con.execute(
            "INSERT INTO refs(version, instance, created_at) VALUES (1, 'ghost', ?)",
            (now - 1000,),
        )
    upload(client, vcf([(500, "A", "G")]), base_version=1)
    app.state.reap_once()
    assert not os.path.exists(_version_dir(app, v1_path)), "dead ref must be ignored"


def test_crash_leftovers_are_swept_and_old_version_survives(client, app):
    upload(client, SAMPLE, base_version=0)
    cfg = app.state.cfg

    # half-published junk: a versions/ dir with data but no DB row
    junk = os.path.join(cfg.data_dir, "versions", "junk-half-set")
    os.makedirs(junk)
    open(os.path.join(junk, GZ_NAME), "wb").write(b"not really bgzf")
    # staging leftover of a dead instance
    dead_staging = os.path.join(cfg.data_dir, "staging", "dead-instance")
    os.makedirs(dead_staging)
    open(os.path.join(dead_staging, "upload.vcf"), "w").write("partial")
    # make them old enough to be sweepable
    old = time.time() - 10 * cfg.orphan_grace
    for p in (junk, dead_staging):
        os.utime(p, (old, old))

    app.state.reap_once()
    assert not os.path.exists(junk)
    assert not os.path.exists(dead_staging)

    # the published version is untouched and still answers queries
    q = client.get("/variants", params={"chrom": "chr1", "start": 100, "end": 103})
    assert q.status_code == 200
    assert [rec.split("\t")[1] for rec in q.json()["records"]] == ["100"]


def test_restart_with_same_data_dir_keeps_version_queryable(tmp_path):
    """Simulated process restart: new app object over the same data dir."""
    cfg1 = make_config(tmp_path, instance_id="inst-a")
    app1 = create_app(cfg1)
    from fastapi.testclient import TestClient

    with TestClient(app1) as c1:
        assert upload(c1, SAMPLE, base_version=0).status_code == 201

    # "restart": fresh app, same data dir, same instance id
    cfg2 = make_config(tmp_path, instance_id="inst-a")
    app2 = create_app(cfg2)
    with TestClient(app2) as c2:
        q = c2.get("/variants", params={"chrom": "chr1", "start": 102, "end": 103})
        assert q.status_code == 200
        assert q.json()["version"] == 1
        assert [rec.split("\t")[1] for rec in q.json()["records"]] == ["100"]


def test_two_instances_share_state(tmp_path):
    """Two app instances (as in docker compose) over one shared data dir."""
    from fastapi.testclient import TestClient

    app1 = create_app(make_config(tmp_path, instance_id="inst-1"))
    app2 = create_app(make_config(tmp_path, instance_id="inst-2"))
    with TestClient(app1) as c1, TestClient(app2) as c2:
        assert upload(c1, SAMPLE, base_version=0).status_code == 201
        # instance 2 sees instance 1's publish
        q = c2.get("/variants", params={"chrom": "chr1", "start": 105, "end": 105})
        assert q.json()["version"] == 1
        assert q.json()["count"] == 1
        # instance 2 replaces; instance 1 sees it
        assert upload(c2, vcf([(500, "A", "G")]), base_version=1).status_code == 201
        q = c1.get("/variants", params={"chrom": "chr1", "start": 1, "end": 1000})
        assert q.json()["version"] == 2
        assert [rec.split("\t")[1] for rec in q.json()["records"]] == ["500"]
        # a stale replacement attempt from instance 1 fails
        assert upload(c1, SAMPLE, base_version=1).status_code == 409

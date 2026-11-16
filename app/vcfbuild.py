"""VCF validation and indexed build (sort -> BGZF -> Tabix)."""
from __future__ import annotations

import os
import re

import pysam

# Strict reading of the spec: REF/ALT may contain A/C/G/T only (uppercase).
ALLELE_RE = re.compile(r"[ACGT]+")
CHROM_RE = re.compile(r"\S+")
# TBI binning supports positions up to 2^29-1 (covers all human chromosomes).
MAX_POS = 2**29 - 1
EXPECTED_HEADER = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

GZ_NAME = "variants.vcf.gz"
TBI_NAME = GZ_NAME + ".tbi"


class VCFError(ValueError):
    """The uploaded VCF violates the acceptance rules."""


def _fail(lineno: int, msg: str):
    raise VCFError(f"line {lineno}: {msg}")


def _parse(upload_path: str):
    """Parse and validate the upload. Returns (meta_lines, records, chrom)."""
    meta_lines = []
    records = []  # (pos, raw line)
    chrom = None
    saw_header = False

    try:
        fh = open(upload_path, "r", encoding="utf-8", errors="strict", newline=None)
    except OSError as exc:  # pragma: no cover - defensive
        raise VCFError(f"cannot read upload: {exc}")
    with fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("##"):
                if saw_header:
                    _fail(lineno, "meta line found after the #CHROM header")
                meta_lines.append(line)
                continue
            if line.startswith("#"):
                if saw_header:
                    _fail(lineno, "duplicate #CHROM header line")
                cols = line.split("\t")
                if cols != EXPECTED_HEADER:
                    _fail(
                        lineno,
                        "column header must be exactly '#CHROM POS ID REF ALT QUAL "
                        "FILTER INFO' (8 columns, no FORMAT/sample columns)",
                    )
                saw_header = True
                continue
            if not saw_header:
                _fail(lineno, "data line before the #CHROM header")
            fields = line.split("\t")
            if len(fields) != 8:
                _fail(lineno, f"expected 8 tab-separated fields, got {len(fields)}")
            c, pos_s, _id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            if not CHROM_RE.fullmatch(c):
                _fail(lineno, "invalid CHROM value")
            if chrom is None:
                chrom = c
            elif c != chrom:
                _fail(lineno, f"multiple chromosomes: expected '{chrom}', found '{c}'")
            if not pos_s.isdigit():
                _fail(lineno, "POS must be a positive integer")
            pos = int(pos_s)
            if pos < 1 or pos > MAX_POS:
                _fail(lineno, f"POS out of supported range 1..{MAX_POS}")
            if not ALLELE_RE.fullmatch(ref):
                _fail(lineno, "REF must contain only A/C/G/T")
            if "," in alt:
                _fail(lineno, "multi-allelic ALT is not allowed (biallelic only)")
            if not ALLELE_RE.fullmatch(alt):
                _fail(lineno, "ALT must contain only A/C/G/T")
            records.append((pos, line))

    if not saw_header:
        raise VCFError("missing #CHROM header line")
    if not records:
        raise VCFError("VCF contains no variant records")
    if chrom is None:  # pragma: no cover - records imply a chromosome
        raise VCFError("could not determine chromosome")
    return meta_lines, records, chrom


def validate_and_build(upload_path: str, staging_dir: str):
    """Validate the uploaded VCF and build sorted BGZF + Tabix index.

    Writes variants.vcf.gz and variants.vcf.gz.tbi into staging_dir.
    Returns (chrom, n_records). Raises VCFError on invalid input.
    """
    try:
        meta_lines, records, chrom = _parse(upload_path)
    except UnicodeDecodeError as exc:
        raise VCFError(f"upload is not valid UTF-8 text: {exc}")

    # Stable sort by POS (single chromosome, ties keep upload order).
    records.sort(key=lambda r: r[0])

    sorted_path = os.path.join(staging_dir, "sorted.vcf")
    with open(sorted_path, "w", encoding="utf-8", newline="\n") as out:
        for h in meta_lines:
            out.write(h + "\n")
        out.write("\t".join(EXPECTED_HEADER) + "\n")
        for _, line in records:
            out.write(line + "\n")

    gz_path = os.path.join(staging_dir, GZ_NAME)
    pysam.tabix_compress(sorted_path, gz_path, force=True)
    os.remove(sorted_path)
    pysam.tabix_index(gz_path, preset="vcf", force=True)

    # Best-effort durability of the finished file set.
    for p in (gz_path, gz_path + ".tbi"):
        try:
            fd = os.open(p, os.O_RDONLY)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)
        except OSError:
            pass
    return chrom, len(records)

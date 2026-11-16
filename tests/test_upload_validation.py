"""Upload validation rules."""
import pytest

from conftest import HEADER, SAMPLE, upload, vcf


def test_first_upload_requires_base_version_zero(client):
    r = upload(client, SAMPLE, base_version=3)
    assert r.status_code == 409
    assert "current version 0" in r.json()["detail"]


def test_replace_requires_fresh_base_version(client):
    assert upload(client, SAMPLE, base_version=0).status_code == 201
    # stale base version -> conflict
    r = upload(client, SAMPLE, base_version=0)
    assert r.status_code == 409
    assert "current version 1" in r.json()["detail"]
    # correct one succeeds
    r = upload(client, SAMPLE, base_version=1)
    assert r.status_code == 201 and r.json()["version"] == 2


BAD_VCFS = {
    "two_chromosomes": "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tG\t.\t.\t.\nchr2\t200\t.\tA\tG\t.\t.\t.\n",
    "sample_columns": "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\n",
    "format_only_nine_columns": "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
        "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\n",
    "ref_with_N": vcf([(100, "AN", "A")]),
    "ref_lowercase": vcf([(100, "acgt", "A")]),
    "alt_multiallelic": vcf([(100, "A", "G,T")]),
    "alt_symbolic": vcf([(100, "A", "<DEL>")]),
    "alt_star": vcf([(100, "A", "*")]),
    "alt_missing": vcf([(100, "A", ".")]),
    "pos_zero": vcf([(0, "A", "G")]),
    "pos_not_a_number": HEADER + "chr1\tabc\t.\tA\tG\t.\t.\t.\n",
    "no_chrom_header": "##fileformat=VCFv4.2\nchr1\t100\t.\tA\tG\t.\t.\t.\n",
    "no_records": HEADER,
    "empty_file": "",
    "seven_fields": HEADER + "chr1\t100\t.\tA\tG\t.\t.\n",
    "nine_fields_data": HEADER + "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\n",
}


@pytest.mark.parametrize("name", sorted(BAD_VCFS))
def test_invalid_vcfs_are_rejected(client, app, name):
    r = upload(client, BAD_VCFS[name], base_version=0)
    assert r.status_code == 400, f"{name}: {r.status_code} {r.text}"
    assert "invalid VCF" in r.json()["detail"]
    # nothing published, nothing left behind
    assert app.state.store.current_version() == 0
    import os
    assert os.listdir(os.path.join(app.state.cfg.data_dir, "versions")) == []


def test_failed_upload_keeps_old_version_queryable(client):
    upload(client, SAMPLE, base_version=0)
    r = upload(client, BAD_VCFS["ref_with_N"], base_version=1)
    assert r.status_code == 400
    q = client.get("/variants", params={"chrom": "chr1", "start": 100, "end": 100})
    assert q.status_code == 200 and q.json()["version"] == 1


def test_oversize_upload_rejected(client, app):
    limit = app.state.cfg.max_upload_bytes
    big = HEADER.encode() + b"chr1\t100\t.\tA\tG\t.\t.\t.\n" + b"#" * limit
    assert len(big) > limit
    r = upload(client, big, base_version=0)
    assert r.status_code == 413
    assert app.state.store.current_version() == 0


def test_binary_upload_rejected_as_invalid_vcf(client):
    import gzip

    r = upload(client, gzip.compress(SAMPLE.encode()), base_version=0)
    assert r.status_code == 400
    assert "invalid VCF" in r.json()["detail"]


def test_exactly_at_limit_accepted(client, app):
    from conftest import CHROM_LINE, META

    limit = app.state.cfg.max_upload_bytes
    record = "chr1\t100\t.\tA\tG\t.\t.\t.\n"
    # pad with a ##meta line (before the #CHROM line) to just below the limit
    pad = limit - len(META.encode()) - len(CHROM_LINE.encode()) - len(record.encode()) - 3
    body = META + "##" + "x" * pad + "\n" + CHROM_LINE + record
    assert len(body.encode()) <= limit
    r = upload(client, body, base_version=0)
    assert r.status_code == 201, r.text

"""API-level tests: upload, validation, interval queries via the index."""
import os

from conftest import SAMPLE, upload, vcf

from app.vcfbuild import GZ_NAME, TBI_NAME


def test_upload_then_query_flow(client):
    r = upload(client, SAMPLE, base_version=0)
    assert r.status_code == 201, r.text
    assert r.json() == {"version": 1, "chrom": "chr1", "n_records": 4}

    q = client.get("/variants", params={"chrom": "chr1", "start": 105, "end": 105})
    assert q.status_code == 200
    body = q.json()
    assert body["version"] == 1
    assert body["count"] == 1
    assert body["records"][0].split("\t")[1] == "105"


def test_published_files_and_index_exist(client, app):
    upload(client, SAMPLE, base_version=0)
    info = app.state.store.version_info()
    d = os.path.join(app.state.cfg.data_dir, info["path"])
    assert os.path.isfile(os.path.join(d, GZ_NAME))
    assert os.path.isfile(os.path.join(d, TBI_NAME))


def test_deletion_starting_before_interval_is_returned(client):
    """The core span requirement: del at POS=100 (REF len 4) covers 100-103."""
    upload(client, SAMPLE, base_version=0)
    for start, end in [(102, 103), (103, 103), (100, 100), (99, 102)]:
        q = client.get("/variants", params={"chrom": "chr1", "start": start, "end": end})
        poss = [int(r.split("\t")[1]) for r in q.json()["records"]]
        assert 100 in poss, f"interval [{start},{end}] must include the deletion at 100"

    q = client.get("/variants", params={"chrom": "chr1", "start": 104, "end": 104})
    assert q.json()["records"] == [], "span ends at 103; 104 must be empty"


def test_long_deletion_span(client):
    upload(client, SAMPLE, base_version=0)
    # deletion at 200 with REF=TAAAA spans 200-204
    q = client.get("/variants", params={"chrom": "chr1", "start": 204, "end": 250})
    poss = [int(r.split("\t")[1]) for r in q.json()["records"]]
    assert poss == [200]
    q = client.get("/variants", params={"chrom": "chr1", "start": 205, "end": 250})
    assert q.json()["records"] == []


def test_insertion_and_boundary_conditions(client):
    upload(client, SAMPLE, base_version=0)
    # insertion at 110 (REF len 1) only intersects [110, 110]
    q = client.get("/variants", params={"chrom": "chr1", "start": 110, "end": 110})
    assert [int(r.split("\t")[1]) for r in q.json()["records"]] == [110]
    q = client.get("/variants", params={"chrom": "chr1", "start": 111, "end": 111})
    assert q.json()["records"] == []
    # whole-range query returns everything, sorted
    q = client.get("/variants", params={"chrom": "chr1", "start": 1, "end": 1000000})
    poss = [int(r.split("\t")[1]) for r in q.json()["records"]]
    assert poss == [100, 105, 110, 200]


def test_unsorted_upload_is_sorted_and_queryable(client):
    body = vcf([(300, "A", "G"), (100, "C", "T"), (200, "G", "A")])
    r = upload(client, body, base_version=0)
    assert r.status_code == 201
    q = client.get("/variants", params={"chrom": "chr1", "start": 1, "end": 1000})
    poss = [int(rec.split("\t")[1]) for rec in q.json()["records"]]
    assert poss == [100, 200, 300]


def test_interval_and_chrom_validation(client):
    upload(client, SAMPLE, base_version=0)
    assert client.get("/variants", params={"chrom": "chr1", "start": 0, "end": 5}).status_code == 422
    assert client.get("/variants", params={"chrom": "chr1", "start": 10, "end": 5}).status_code == 400
    r = client.get("/variants", params={"chrom": "chr2", "start": 1, "end": 5})
    assert r.status_code == 400 and "chr1" in r.json()["detail"]


def test_query_without_data_is_404(client):
    assert client.get("/variants", params={"chrom": "chr1", "start": 1, "end": 5}).status_code == 404


def test_version_endpoint(client):
    assert client.get("/version").json() == {"published": False, "current_version": 0}
    upload(client, SAMPLE, base_version=0)
    body = client.get("/version").json()
    assert body["published"] is True
    assert body["current_version"] == 1
    assert body["chrom"] == "chr1"
    assert body["n_records"] == 4

import pytest
import requests

import sra_tools

DOCSET = b"""<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_PACKAGE_SET>
  <EXPERIMENT_PACKAGE>
    <EXPERIMENT accession="ERX1">
      <STUDY_REF accession="ERP1"/>
      <DESIGN><LIBRARY_DESCRIPTOR><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN>
      <PLATFORM><ILLUMINA><INSTRUMENT_MODEL>NextSeq 500</INSTRUMENT_MODEL></ILLUMINA></PLATFORM>
    </EXPERIMENT>
    <SAMPLE accession="ERS1"><TITLE>s1</TITLE></SAMPLE>
    <RUN_SET><RUN accession="ERR1" total_bases="600" total_spots="2" size="100"/></RUN_SET>
  </EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>
"""

EMPTY_SET = b'<?xml version="1.0" encoding="UTF-8"?>\n<EXPERIMENT_PACKAGE_SET>\n</EXPERIMENT_PACKAGE_SET>\n'

ID_LIST_EMPTY = (b'<?xml version="1.0" encoding="UTF-8" ?>\n<eFetchResult>\n'
                 b'\t<ERROR>ID list is empty! Possibly it has no correct IDs.</ERROR>\n</eFetchResult>\n')


class FakeResponse:
    def __init__(self, status, content, headers=None):
        self.status_code = status
        self.content = content
        self.headers = headers or {}


@pytest.fixture
def efetch(monkeypatch):
    """Serve a scripted sequence of responses (or exceptions) to requests.get."""
    script = []
    calls = []
    sleeps = []

    def fake_get(url, params=None, timeout=None):
        calls.append(params)
        item = script.pop(0)
        if isinstance(item, Exception):
            raise item
        return item

    monkeypatch.setattr(sra_tools.requests, "get", fake_get)
    monkeypatch.setattr(sra_tools.time, "sleep", sleeps.append)
    return script, calls, sleeps


def test_found_on_first_try(efetch):
    script, calls, _ = efetch
    script.append(FakeResponse(200, DOCSET))
    md = sra_tools.get_accession_metadata("ERR1", None)
    assert [r["run_id"] for r in md] == ["ERR1"]
    assert len(calls) == 1


def test_empty_answer_is_reasked(efetch, capsys):
    # The job 23963680 shape: a well-formed answer naming no runs, then the run.
    script, calls, _ = efetch
    script += [FakeResponse(200, EMPTY_SET), FakeResponse(200, DOCSET)]
    md = sra_tools.get_accession_metadata("ERR1", None)
    assert [r["run_id"] for r in md] == ["ERR1"]
    assert len(calls) == 2
    assert "names no runs" in capsys.readouterr().err


def test_unknown_id_returns_empty_after_reasking(efetch):
    # --metaonly callers still get [] for an id NCBI does not know.
    script, calls, _ = efetch
    script += [FakeResponse(400, ID_LIST_EMPTY)] * sra_tools.EMPTY_ANSWER_ATTEMPTS
    assert sra_tools.get_accession_metadata("ERR404", None) == []
    assert len(calls) == sra_tools.EMPTY_ANSWER_ATTEMPTS


def test_rate_limit_retried_honoring_retry_after(efetch):
    script, calls, sleeps = efetch
    script += [FakeResponse(429, b"slow down", {"Retry-After": "7"}), FakeResponse(200, DOCSET)]
    assert sra_tools.get_accession_metadata("ERR1", None)
    assert sleeps == [7.0]


def test_server_error_retried(efetch):
    script, calls, _ = efetch
    script += [FakeResponse(503, b"<html>busy</html>"), FakeResponse(200, DOCSET)]
    assert sra_tools.get_accession_metadata("ERR1", None)
    assert len(calls) == 2


def test_network_error_retried(efetch):
    script, calls, _ = efetch
    script += [requests.exceptions.ConnectionError("reset"), FakeResponse(200, DOCSET)]
    assert sra_tools.get_accession_metadata("ERR1", None)


def test_persistent_failure_raises_with_status_and_body(efetch):
    script, _, _ = efetch
    script += [FakeResponse(503, b"backend down")] * sra_tools.EFETCH_ATTEMPTS
    with pytest.raises(RuntimeError, match="HTTP 503: backend down"):
        sra_tools.get_accession_metadata("ERR1", None)


def test_non_retryable_status_raises_immediately(efetch):
    script, calls, _ = efetch
    script.append(FakeResponse(403, b"error code: 1010"))
    with pytest.raises(RuntimeError, match="HTTP 403"):
        sra_tools.get_accession_metadata("ERR1", None)
    assert len(calls) == 1


def test_unparseable_200_raises_after_reasking(efetch):
    script, calls, _ = efetch
    script += [FakeResponse(200, b"<html><body>oops")] * sra_tools.EMPTY_ANSWER_ATTEMPTS
    with pytest.raises(RuntimeError, match="unparseable"):
        sra_tools.get_accession_metadata("ERR1", None)
    assert len(calls) == sra_tools.EMPTY_ANSWER_ATTEMPTS


def test_download_refuses_empty_metadata(monkeypatch, tmp_path):
    monkeypatch.setattr(sra_tools, "get_accession_metadata", lambda acc, f: [])
    ran = []
    monkeypatch.setattr(sra_tools, "retry_subprocess_check_output", lambda *a: ran.append(a))
    with pytest.raises(RuntimeError, match="No SRA run metadata found for ERR1"):
        sra_tools.download_sra_data("fasterq-dump", str(tmp_path), "ERR1", False, False, None, None)
    assert ran == []   # never reaches a bare `prefetch`


def test_metaonly_still_writes_empty_list(monkeypatch, tmp_path):
    monkeypatch.setattr(sra_tools, "get_accession_metadata", lambda acc, f: [])
    out = tmp_path / "md.json"
    sra_tools.download_sra_data("fasterq-dump", None, "ERR404", True, False, str(out), None)
    assert out.read_text().strip() == "[]"


def test_subprocess_failure_raises_runtime_error(monkeypatch):
    monkeypatch.setattr(sra_tools.subprocess, "call", lambda cmd, stderr=None: 3)
    monkeypatch.setattr(sra_tools.time, "sleep", lambda s: None)
    with pytest.raises(RuntimeError, match="last exit code 3"):
        sra_tools.retry_subprocess_check_output(["prefetch", "ERR1"], 2, 0)


def test_gzip_step_runs(monkeypatch, tmp_path):
    # `fq_files.count > 0` used to raise TypeError on every --gzip download.
    monkeypatch.setattr(sra_tools, "retry_subprocess_check_output", lambda *a: None)
    (tmp_path / "SRR1.fastq").write_text("@r\nA\n+\nI\n")
    gz = []
    monkeypatch.setattr(sra_tools.subprocess, "call", lambda cmd: gz.append(cmd))
    md = [{"run_id": "SRR1", "platform_name": "ILLUMINA"}]
    sra_tools.download_fastq_files("fasterq-dump", str(tmp_path), md, gzip_output=True)
    assert gz and gz[0][0] == "gzip"

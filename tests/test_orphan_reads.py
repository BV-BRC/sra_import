import os

import pytest

import sra_tools


def write_fastq(path, n):
    with open(path, "w") as fh:
        for i in range(n):
            fh.write(f"@r{i}\nACGT\n+\nIIII\n")


def make_run(d, run, pairs, orphans):
    write_fastq(os.path.join(d, f"{run}_1.fastq"), pairs)
    write_fastq(os.path.join(d, f"{run}_2.fastq"), pairs)
    if orphans is not None:
        write_fastq(os.path.join(d, f"{run}.fastq"), orphans)


def test_orphans_below_limit_are_discarded(tmp_path):
    make_run(tmp_path, "ERR1", 1000, 3)
    assert sra_tools.discard_orphan_reads(str(tmp_path), "ERR1", total_spots=1003, max_fraction=0.05) == 3
    assert sorted(os.listdir(tmp_path)) == ["ERR1_1.fastq", "ERR1_2.fastq"]


def test_orphans_above_limit_refuse_and_keep_files(tmp_path):
    make_run(tmp_path, "ERR1", 90, 10)
    with pytest.raises(RuntimeError, match="10 of 100 spots"):
        sra_tools.discard_orphan_reads(str(tmp_path), "ERR1", total_spots=100, max_fraction=0.05)
    assert os.path.exists(tmp_path / "ERR1.fastq")


def test_total_spots_derived_from_files_when_metadata_lacks_it(tmp_path):
    make_run(tmp_path, "ERR1", 90, 10)
    with pytest.raises(RuntimeError, match="10 of 100 spots"):
        sra_tools.discard_orphan_reads(str(tmp_path), "ERR1", total_spots=None, max_fraction=0.05)


def test_clean_paired_run_untouched(tmp_path):
    make_run(tmp_path, "SRR1", 10, None)
    assert sra_tools.discard_orphan_reads(str(tmp_path), "SRR1", total_spots=10) == 0
    assert sorted(os.listdir(tmp_path)) == ["SRR1_1.fastq", "SRR1_2.fastq"]


def test_single_end_run_untouched(tmp_path):
    # --split-3 and --split-files both write a single-end run to <run>.fastq;
    # with no pair beside it, that file is the data and must survive.
    write_fastq(tmp_path / "SRR2.fastq", 10)
    assert sra_tools.discard_orphan_reads(str(tmp_path), "SRR2", total_spots=10) == 0
    assert os.listdir(tmp_path) == ["SRR2.fastq"]


def test_other_runs_in_same_directory_untouched(tmp_path):
    # An SRX/SRP accession downloads several runs into one directory.
    make_run(tmp_path, "SRR1", 100, 1)
    write_fastq(tmp_path / "SRR2.fastq", 10)
    sra_tools.discard_orphan_reads(str(tmp_path), "SRR1", total_spots=101)
    assert sorted(os.listdir(tmp_path)) == ["SRR1_1.fastq", "SRR1_2.fastq", "SRR2.fastq"]


def test_limit_from_environment(tmp_path, monkeypatch):
    make_run(tmp_path, "ERR1", 90, 10)
    monkeypatch.setenv("P3_SRA_MAX_ORPHAN_FRACTION", "1")
    assert sra_tools.discard_orphan_reads(str(tmp_path), "ERR1", total_spots=100) == 10


def test_default_limit(monkeypatch):
    monkeypatch.delenv("P3_SRA_MAX_ORPHAN_FRACTION", raising=False)
    assert sra_tools.max_orphan_fraction() == sra_tools.DEFAULT_MAX_ORPHAN_FRACTION


def test_fasterq_dump_uses_split_3(tmp_path, monkeypatch):
    calls = []

    def fake_retry(cmd, n, sleep):
        calls.append(cmd)
        if cmd[0] == "fasterq-dump":
            make_run(tmp_path, "ERR1", 1000, 2)

    monkeypatch.setattr(sra_tools, "retry_subprocess_check_output", fake_retry)
    md = [{"run_id": "ERR1", "platform_name": "ILLUMINA", "total_spots": 1002}]
    sra_tools.download_fastq_files("fasterq-dump", str(tmp_path), md, gzip_output=False)

    dump = [c for c in calls if c[0] == "fasterq-dump"][0]
    assert "--split-3" in dump and "--split-files" not in dump
    assert sorted(os.listdir(tmp_path)) == ["ERR1_1.fastq", "ERR1_2.fastq"]

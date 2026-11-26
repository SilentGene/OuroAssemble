"""workflow/scripts/resume_guard.sh: resume only when reads and options are unchanged."""

import os
import subprocess
from pathlib import Path

GUARD = Path(__file__).parents[1] / "ouroassemble" / "workflow" / "scripts" / "resume_guard.sh"


def guard(workdir: Path, reads: Path, *params: str) -> str:
    res = subprocess.run(
        ["bash", str(GUARD), str(workdir), str(reads), *params],
        capture_output=True, text=True, check=True,
    )  # fmt: skip
    return res.stdout.strip()


def test_fresh_then_resume(tmp_path):
    reads = tmp_path / "reads.fq.gz"
    reads.write_text("x")
    work = tmp_path / "asm"
    assert guard(work, reads, "-k 21") == "fresh"
    (work / "stage1.done").write_text("")
    assert guard(work, reads, "-k 21") == "resume"
    assert (work / "stage1.done").exists()


def test_changed_options_or_reads_clear_workdir(tmp_path):
    reads = tmp_path / "reads.fq.gz"
    reads.write_text("x")
    work = tmp_path / "asm"
    guard(work, reads, "-k 21")
    (work / "sub").mkdir()
    (work / "sub" / "stale").write_text("")
    assert guard(work, reads, "-k 25") == "fresh"
    assert not (work / "sub").exists()

    (work / "stale").write_text("")
    reads.write_text("new reads")
    os.utime(reads, (1, 1))
    assert guard(work, reads, "-k 25") == "fresh"
    assert not (work / "stale").exists()


def test_legacy_workdir_is_adopted(tmp_path):
    reads = tmp_path / "reads.fq.gz"
    reads.write_text("x")
    work = tmp_path / "asm"
    work.mkdir()
    (work / "overlap.000001.ovb").write_text("")
    assert guard(work, reads, "5m") == "resume"
    assert (work / "overlap.000001.ovb").exists()
    assert guard(work, reads, "5m") == "resume"


def test_refuses_dangerous_workdir(tmp_path):
    reads = tmp_path / "reads.fq.gz"
    reads.write_text("x")
    res = subprocess.run(["bash", str(GUARD), "/", str(reads)], capture_output=True)
    assert res.returncode == 2

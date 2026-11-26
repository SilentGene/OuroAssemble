import yaml
from typer.testing import CliRunner

from ouroassemble import DEFAULT_CONFIG
from ouroassemble.cli import app

runner = CliRunner()


def test_init_config(tmp_path):
    out = tmp_path / "c.yaml"
    res = runner.invoke(app, ["init-config", "-o", str(out)])
    assert res.exit_code == 0, res.output
    assert yaml.safe_load(out.read_text()) == yaml.safe_load(DEFAULT_CONFIG.read_text())


def test_run_requires_reads(tmp_path):
    res = runner.invoke(app, ["run", "-o", str(tmp_path / "o")])
    assert res.exit_code == 1


def test_run_rejects_unknown_assembler(tmp_path):
    reads = tmp_path / "r.fq.gz"
    reads.write_bytes(b"")
    res = runner.invoke(app, ["run", "-r", str(reads), "-a", "spades"])
    assert res.exit_code == 1


def test_dry_run_builds_dag(tmp_path):
    reads = tmp_path / "reads.fq.gz"
    reads.write_bytes(b"")
    out = tmp_path / "out"
    res = runner.invoke(
        app,
        ["run", "-r", str(reads), "-o", str(out), "-a", "metaflye,hicanu",
         "--skip-checkm2", "--skip-gtdbtk", "--skip-genomad", "--dry-run"],
    )  # fmt: skip
    assert res.exit_code == 0, res.output
    cfg = yaml.safe_load((out / "ouroassemble.config.yaml").read_text())
    assert cfg["assemblers"] == ["metaflye", "hicanu"]
    assert cfg["classify"]["skip_gtdbtk"] is True


def test_memory_cap_applies_to_every_job(tmp_path):
    import re
    import subprocess
    import sys

    from ouroassemble import SNAKEFILE

    reads = tmp_path / "reads.fq.gz"
    reads.write_bytes(b"")
    out = tmp_path / "out"
    res = runner.invoke(
        app,
        ["run", "-r", str(reads), "-o", str(out), "-m", "20",
         "--skip-checkm2", "--skip-gtdbtk", "--skip-genomad", "--dry-run"],
    )  # fmt: skip
    assert res.exit_code == 0, res.output
    cfg = yaml.safe_load((out / "ouroassemble.config.yaml").read_text())
    assert cfg["max_memory_gb"] == 20
    dry = subprocess.run(
        [sys.executable, "-m", "snakemake", "-s", str(SNAKEFILE), "-d", str(out),
         "--configfile", str(out / "ouroassemble.config.yaml"), "-n", "--cores", "1"],
        capture_output=True, text=True,
    )  # fmt: skip
    mems = [int(m) for m in re.findall(r"mem_mb=(\d+)", dry.stdout + dry.stderr)]
    assert mems and max(mems) == 20480
    assert 4000 in mems  # small jobs keep their own (lower) value

"""Command-line interface for OuroAssemble."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
from enum import StrEnum
from pathlib import Path
from typing import Annotated

import typer
import yaml
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from . import DEFAULT_CONFIG, PACKAGE_DIR, SNAKEFILE, __version__

ASSEMBLERS = ("metamdbg", "hifiasm_meta", "metaflye", "hicanu")
DOWNLOAD_SNAKEFILE = PACKAGE_DIR / "workflow" / "download_db.smk"

app = typer.Typer(
    name="ouroassemble",
    help="Recover and validate circular genomes from PacBio HiFi metagenomes.",
    no_args_is_help=True,
    rich_markup_mode="rich",
    add_completion=False,
    context_settings={"help_option_names": ["-h", "--help"]},
)
console = Console()
err = Console(stderr=True)


class Database(StrEnum):
    checkm2 = "checkm2"
    gtdbtk = "gtdbtk"
    genomad = "genomad"
    all = "all"


def fail(message: str, title: str = "Error") -> None:
    err.print(Panel(message, title=f"[bold red]{title}", border_style="red"))
    raise typer.Exit(code=1)


def deep_update(base: dict, override: dict) -> dict:
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = value
    return base


def snakemake_cmd(
    snakefile: Path,
    workdir: Path,
    configfile: Path | None,
    cores: int,
    profile: str | None,
    dry_run: bool,
    conda_prefix: str | None,
    extra: list[str],
) -> list[str]:
    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--directory",
        str(workdir),
        "--software-deployment-method",
        "conda",
        "--rerun-incomplete",
        "--printshellcmds",
        "--show-failed-logs",
    ]
    if configfile:
        cmd += ["--configfile", str(configfile)]
    if profile:
        cmd += ["--profile", profile]
    else:
        cmd += ["--cores", str(cores)]
    conda_prefix = conda_prefix or os.environ.get("OUROASSEMBLE_CONDA_PREFIX")
    if conda_prefix:
        cmd += ["--conda-prefix", str(Path(conda_prefix).expanduser().resolve())]
    if dry_run:
        cmd += ["--dry-run"]
    return cmd + list(extra)


JOBSCRIPT_TEMPLATE = """#!/bin/bash
# properties = {{properties}}
# Written by OuroAssemble: cluster jobs often start with a minimal PATH, so put the
# driver environment (python, snakemake, conda) first; Snakemake needs `conda` on the
# compute node to activate the per-rule environments.
export PATH="__BINDIR__:$PATH"
{{exec_job}}
"""


def write_jobscript(outdir: Path) -> Path:
    """Cluster job script that exposes the driver environment's bin/ on compute nodes."""
    path = outdir / ".ouroassemble_jobscript.sh"
    bindir = str(Path(sys.executable).parent)
    path.write_text(
        JOBSCRIPT_TEMPLATE.replace("{{", "{").replace("}}", "}").replace("__BINDIR__", bindir)
    )
    return path


def run_snakemake(cmd: list[str]) -> None:
    console.print(f"[dim]$ {shlex.join(cmd)}[/dim]")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        fail(f"Snakemake exited with code {result.returncode}. See the log above.", "Run failed")


@app.command(
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
    epilog="Arguments after [bold]--[/bold] are passed to Snakemake unchanged, "
    "e.g. [cyan]ouroassemble run -r reads.fq.gz -o out -- --keep-going[/cyan]",
)
def run(
    ctx: typer.Context,
    reads: Annotated[
        Path | None, typer.Option("--reads", "-r", help="PacBio HiFi reads (FASTQ/FASTA[.gz]).")
    ] = None,
    outdir: Annotated[Path, typer.Option("--outdir", "-o", help="Output directory.")] = Path(
        "ouroassemble_out"
    ),
    configfile: Annotated[
        Path | None, typer.Option("--config", "-c", help="YAML config (see init-config).")
    ] = None,
    sample: Annotated[str | None, typer.Option("--sample", "-s", help="Sample name.")] = None,
    threads: Annotated[int, typer.Option("--threads", "-t", help="Total cores (local).")] = 8,
    memory: Annotated[
        float | None,
        typer.Option(
            "--memory",
            "-m",
            min=1,
            help="Maximum RAM in GB. Caps every job's memory; when running locally the "
            "jobs running at the same time also stay below it in total.",
        ),
    ] = None,
    assemblers: Annotated[
        str | None,
        typer.Option(
            "--assemblers", "-a", help=f"Comma-separated subset of {','.join(ASSEMBLERS)}."
        ),
    ] = None,
    checkm2_db: Annotated[Path | None, typer.Option(help="CheckM2 .dmnd database.")] = None,
    gtdbtk_db: Annotated[Path | None, typer.Option(help="GTDB-Tk reference data dir.")] = None,
    genomad_db: Annotated[Path | None, typer.Option(help="geNomad database dir.")] = None,
    skip_checkm2: Annotated[bool, typer.Option(help="Skip CheckM2.")] = False,
    skip_gtdbtk: Annotated[bool, typer.Option(help="Skip GTDB-Tk.")] = False,
    skip_genomad: Annotated[bool, typer.Option(help="Skip geNomad.")] = False,
    profile: Annotated[
        str | None, typer.Option("--profile", "-p", help="Snakemake profile (e.g. cluster).")
    ] = None,
    conda_prefix: Annotated[
        str | None,
        typer.Option(help="Shared dir for rule conda envs [env: OUROASSEMBLE_CONDA_PREFIX]."),
    ] = None,
    dry_run: Annotated[bool, typer.Option("--dry-run", "-n", help="Only show jobs.")] = False,
    resume: Annotated[
        bool,
        typer.Option(
            "--resume",
            help="Continue an interrupted run (e.g. killed at walltime) in the same output "
            "directory: removes Snakemake's stale lock, keeps every finished step even if "
            "the config or OuroAssemble changed since (only missing files are made), and "
            "lets HiCanu / metaFlye / hifiasm-meta pick up where they stopped. Threads and "
            "memory may be changed. Make sure no other run uses this directory.",
        ),
    ] = False,
) -> None:
    """Run the full pipeline on one HiFi metagenome."""
    cfg: dict = {}
    if configfile:
        if not configfile.exists():
            fail(f"Config file not found: {configfile}")
        cfg = yaml.safe_load(configfile.read_text()) or {}

    overrides: dict = {}
    if reads:
        overrides["reads"] = str(reads)
    if sample:
        overrides["sample"] = sample
    if memory:
        overrides["max_memory_gb"] = memory
    # `threads: all` jobs (HiCanu by default) take every core of a local run
    overrides["local_cores"] = None if profile else threads
    if assemblers:
        chosen = [a.strip() for a in assemblers.split(",") if a.strip()]
        bad = [a for a in chosen if a not in ASSEMBLERS]
        if bad:
            fail(f"Unknown assembler(s): {', '.join(bad)}. Choose from {', '.join(ASSEMBLERS)}.")
        overrides["assemblers"] = chosen
    db_opts = {"checkm2": checkm2_db, "gtdbtk": gtdbtk_db, "genomad": genomad_db}
    dbs = {k: str(v.resolve()) for k, v in db_opts.items() if v}
    if dbs:
        overrides["databases"] = dbs
    skip_opts = {"checkm2": skip_checkm2, "gtdbtk": skip_gtdbtk, "genomad": skip_genomad}
    skips = {f"skip_{k}": True for k, v in skip_opts.items() if v}
    if skips:
        overrides["classify"] = skips
    cfg = deep_update(cfg, overrides)

    if not cfg.get("reads"):
        fail("No input reads. Use [bold]--reads[/bold] or set [bold]reads:[/bold] in the config.")
    reads_path = Path(cfg["reads"]).expanduser()
    if not reads_path.is_absolute() and configfile and not reads:
        reads_path = configfile.parent / reads_path
    reads_path = reads_path.resolve()
    if not reads_path.exists():
        fail(f"Reads not found: {reads_path}")
    cfg["reads"] = str(reads_path)

    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    run_config = outdir / "ouroassemble.config.yaml"
    run_config.write_text(yaml.safe_dump(cfg, sort_keys=False))

    summary = Table.grid(padding=(0, 2))
    summary.add_row("[bold]reads", str(reads_path))
    summary.add_row("[bold]output", str(outdir))
    summary.add_row("[bold]assemblers", ", ".join(cfg.get("assemblers", ASSEMBLERS)))
    summary.add_row(
        "[bold]execution", f"profile {profile}" if profile else f"local, {threads} cores"
    )
    if resume:
        summary.add_row("[bold]mode", "resume")
    max_mem = cfg.get("max_memory_gb")
    if max_mem:
        scope = "per job" if profile else "per job and in total"
        summary.add_row("[bold]memory", f"≤ {max_mem:g} GB {scope}")
    console.print(Panel(summary, title=f"[bold cyan]OuroAssemble {__version__}", expand=False))

    extra = list(ctx.args)
    if max_mem and not profile and "--resources" not in extra:
        # local: keep the summed memory of concurrently running jobs below the limit
        extra += ["--resources", f"mem_mb={int(max_mem * 1024)}"]
    if profile and not any(a in ("--jobscript", "--js") for a in extra):
        extra += ["--jobscript", str(write_jobscript(outdir))]
    if resume:
        if "--rerun-triggers" not in extra:
            # only missing/outdated files count, not code, params or software-env changes
            extra += ["--rerun-triggers", "mtime"]
        locks = outdir / ".snakemake" / "locks"
        if not dry_run and locks.is_dir() and any(locks.iterdir()):
            console.print("[yellow]Removing the lock left by the interrupted run")
            unlock = [SNAKEFILE, outdir, run_config, threads, None, False, conda_prefix]
            run_snakemake(snakemake_cmd(*unlock, ["--unlock"]))
    cmd = snakemake_cmd(
        SNAKEFILE, outdir, run_config, threads, profile, dry_run, conda_prefix, extra
    )
    run_snakemake(cmd)
    if not dry_run:
        console.print(
            Panel(
                f"Report: [bold]{outdir / 'final' / 'report.html'}[/bold]\n"
                f"Table:  {outdir / 'final' / 'circular_genomes.tsv'}",
                title="[bold green]Done",
                expand=False,
            )
        )


@app.command("init-config")
def init_config(
    output: Annotated[Path, typer.Option("--output", "-o", help="Where to write.")] = Path(
        "ouroassemble.config.yaml"
    ),
    force: Annotated[bool, typer.Option("--force", "-f", help="Overwrite.")] = False,
) -> None:
    """Write an annotated copy of the default configuration."""
    if output.exists() and not force:
        fail(f"{output} exists (use --force to overwrite).")
    shutil.copyfile(DEFAULT_CONFIG, output)
    console.print(f"[green]Wrote[/green] {output}")


@app.command("download-db")
def download_db(
    database: Annotated[Database, typer.Argument(help="Which database to fetch.")] = Database.all,
    directory: Annotated[Path, typer.Option("--dir", "-d", help="Target directory.")] = Path(
        "ouroassemble_db"
    ),
    threads: Annotated[int, typer.Option("--threads", "-t")] = 4,
    conda_prefix: Annotated[
        str | None, typer.Option(help="Shared dir for rule conda envs.")
    ] = None,
) -> None:
    """Download CheckM2, GTDB-Tk and/or geNomad databases (large: GTDB-Tk > 100 GB)."""
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    names = ["checkm2", "gtdbtk", "genomad"] if database == Database.all else [database.value]
    targets = [f"{n}.done" for n in names]
    cmd = snakemake_cmd(DOWNLOAD_SNAKEFILE, directory, None, threads, None, False, conda_prefix, [])
    run_snakemake(cmd + targets)

    paths = {}
    for n in names:
        marker = directory / f"{n}.done"
        paths[n] = marker.read_text().strip() if marker.exists() else "?"
    console.print(
        Panel(
            "Add to your config:\n\n[bold]databases:[/bold]\n"
            + "\n".join(f"  {k}: {v}" for k, v in paths.items()),
            title="[bold green]Databases ready",
            expand=False,
        )
    )


@app.command()
def version() -> None:
    """Print the version."""
    console.print(f"OuroAssemble {__version__}")


if __name__ == "__main__":
    app()

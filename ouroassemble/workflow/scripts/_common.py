"""Small dependency-free helpers shared by the OuroAssemble workflow scripts."""

from __future__ import annotations

import csv
import gzip
import logging
import sys
from collections.abc import Iterable, Iterator
from pathlib import Path

LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"


def setup_logging() -> logging.Logger:
    logging.basicConfig(stream=sys.stdout, level=logging.INFO, format=LOG_FORMAT)
    return logging.getLogger("ouroassemble")


def open_text(path: str | Path, mode: str = "rt"):
    """Open plain or gzip-compressed text files transparently."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode.replace("t", "") or "r")


def read_fasta(path: str | Path) -> Iterator[tuple[str, str, str]]:
    """Yield (id, description, sequence); description is the full header line."""
    name = desc = None
    chunks: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, "".join(chunks)
                desc = line[1:].strip()
                name = desc.split()[0] if desc else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, desc, "".join(chunks)


def write_fasta(records: Iterable[tuple[str, str]], path: str | Path, width: int = 80) -> int:
    """Write (header, sequence) records. Returns the number of records written."""
    n = 0
    with open(path, "w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i : i + width] + "\n")
            n += 1
    return n


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open_text(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(rows: Iterable[dict], path: str | Path, columns: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _fmt(row.get(k, "")) for k in columns})


def _fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.4g}" if abs(value) < 1e6 else f"{value:.0f}"
    return str(value)


def to_float(value, default: float | None = None) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def n50(lengths: list[int]) -> int:
    total = sum(lengths)
    acc = 0
    for length in sorted(lengths, reverse=True):
        acc += length
        if acc * 2 >= total:
            return length
    return 0

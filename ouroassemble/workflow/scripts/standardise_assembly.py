#!/usr/bin/env python3
"""Standardise one assembly and extract its circular contig candidates.

Assembler-specific circularity markers
  metamdbg      FASTA header contains ``circular=yes``
  hifiasm_meta  GFA ``S`` line whose contig ID ends with ``c`` (e.g. s1.ctg000001c);
                optional ``*.rescue.fa`` records with IDs ending in ``c`` (source=rescue)
  metaflye      ``assembly_info.txt`` column 4 (circ.) == ``Y``; column 5 (repeat) kept
  hicanu        FASTA header contains ``suggestCircular=yes``

Outputs
  contigs.fa               all contigs, IDs renamed to ``{assembler}__{original_id}``
  contigs.tsv              per-contig metadata (length, depth, circular flag, ...)
  circular_candidates.fa   circular contigs >= --min-circular-len
"""

from __future__ import annotations

import argparse
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

from _common import read_fasta, setup_logging, write_fasta, write_tsv

COLUMNS = [
    "contig_id",
    "assembler",
    "orig_id",
    "length",
    "depth",
    "circular",
    "candidate",
    "source",
    "repeat",
    "notes",
]
ASSEMBLERS = ("metamdbg", "hifiasm_meta", "metaflye", "hicanu")


@dataclass
class Contig:
    orig_id: str
    seq: str
    circular: bool
    depth: float | None = None
    source: str = "main"
    repeat: str = ""
    notes: dict[str, str] = field(default_factory=dict)


def _tags(description: str) -> dict[str, str]:
    """key=value tokens of a FASTA description."""
    return dict(tok.split("=", 1) for tok in description.split()[1:] if "=" in tok)


def parse_metamdbg(fasta: str) -> list[Contig]:
    contigs = []
    for name, desc, seq in read_fasta(fasta):
        tags = _tags(desc)
        depth = tags.get("coverage")
        contigs.append(
            Contig(
                orig_id=name,
                seq=seq,
                circular=tags.get("circular", "").lower() == "yes",
                depth=float(depth) if depth else None,
            )
        )
    return contigs


def parse_gfa_segments(gfa: str) -> list[tuple[str, str, dict[str, str]]]:
    segments = []
    with open(gfa) as fh:
        for line in fh:
            if not line.startswith("S\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            tags = {}
            for tag in parts[3:]:
                bits = tag.split(":", 2)
                if len(bits) == 3:
                    tags[bits[0]] = bits[2]
            segments.append((parts[1], parts[2], tags))
    return segments


def parse_hifiasm_meta(gfa: str, rescue: str | None = None) -> list[Contig]:
    contigs = []
    for name, seq, tags in parse_gfa_segments(gfa):
        if seq == "*":
            raise ValueError(f"GFA segment {name} has no sequence stored")
        depth = tags.get("rd") or tags.get("dp")
        contigs.append(
            Contig(
                orig_id=name,
                seq=seq,
                circular=name.endswith("c"),
                depth=float(depth) if depth else None,
            )
        )
    if rescue and Path(rescue).exists():
        seen = {c.orig_id for c in contigs}
        for name, _desc, seq in read_fasta(rescue):
            if name in seen:
                continue
            contigs.append(
                Contig(
                    orig_id=name,
                    seq=seq,
                    circular=name.endswith("c"),
                    source="rescue",
                    notes={"needs_extra_validation": "rescue"},
                )
            )
    return contigs


def parse_flye_info(info: str) -> dict[str, dict[str, str]]:
    rows = {}
    with open(info) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            rows[parts[0]] = {
                "length": parts[1],
                "cov": parts[2],
                "circ": parts[3],
                "repeat": parts[4],
                "mult": parts[5] if len(parts) > 5 else "",
            }
    return rows


def parse_metaflye(fasta: str, info: str) -> list[Contig]:
    meta = parse_flye_info(info)
    contigs = []
    for name, _desc, seq in read_fasta(fasta):
        row = meta.get(name)
        if row is None:
            raise ValueError(f"{name} missing from {info}")
        contigs.append(
            Contig(
                orig_id=name,
                seq=seq,
                circular=row["circ"].upper() == "Y",
                depth=float(row["cov"]) if row["cov"] else None,
                repeat=row["repeat"],
                notes={"mult": row["mult"]} if row["mult"] else {},
            )
        )
    return contigs


def parse_hicanu(fasta: str) -> list[Contig]:
    contigs = []
    for name, desc, seq in read_fasta(fasta):
        tags = _tags(desc)
        notes = {k: tags[k] for k in ("reads", "class", "suggestBubble", "trim") if k in tags}
        contigs.append(
            Contig(
                orig_id=name,
                seq=seq,
                circular=tags.get("suggestCircular", "").lower() == "yes",
                repeat=tags.get("suggestRepeat", ""),
                notes=notes,
            )
        )
    return contigs


def safe_id(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", text)


def standardise(assembler: str, contigs: list[Contig], min_circular_len: int) -> list[dict]:
    rows = []
    for c in contigs:
        rows.append(
            {
                **asdict(c),
                "contig_id": f"{assembler}__{safe_id(c.orig_id)}",
                "assembler": assembler,
                "length": len(c.seq),
                "candidate": c.circular and len(c.seq) >= min_circular_len,
                "notes": ";".join(f"{k}={v}" for k, v in c.notes.items()),
            }
        )
    return rows


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--assembler", required=True, choices=ASSEMBLERS)
    ap.add_argument("--fasta")
    ap.add_argument("--gfa")
    ap.add_argument("--info", help="metaFlye assembly_info.txt")
    ap.add_argument("--rescue", help="hifiasm-meta *.rescue.fa (optional)")
    ap.add_argument("--min-circular-len", type=int, default=1000)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-candidates", required=True)
    args = ap.parse_args(argv)
    log = setup_logging()

    if args.assembler == "metamdbg":
        contigs = parse_metamdbg(args.fasta)
    elif args.assembler == "hifiasm_meta":
        contigs = parse_hifiasm_meta(args.gfa, args.rescue)
    elif args.assembler == "metaflye":
        contigs = parse_metaflye(args.fasta, args.info)
    else:
        contigs = parse_hicanu(args.fasta)

    rows = standardise(args.assembler, contigs, args.min_circular_len)
    write_fasta(((r["contig_id"], r["seq"]) for r in rows), args.out_fasta)
    write_tsv(rows, args.out_tsv, COLUMNS)
    n_cand = write_fasta(
        ((r["contig_id"], r["seq"]) for r in rows if r["candidate"]), args.out_candidates
    )
    n_circ = sum(r["circular"] for r in rows)
    log.info(
        "%s: %d contigs, %d flagged circular, %d candidates >= %d bp",
        args.assembler,
        len(rows),
        n_circ,
        n_cand,
        args.min_circular_len,
    )


if __name__ == "__main__":
    main()

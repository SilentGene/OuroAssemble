#!/usr/bin/env python3
"""Merge the assemblies of all assemblers into one non-redundant assembly.

Contigs are only *selected*, never cut or joined, so no new chimeras can arise.

``pool``    collect contigs from every assembler and rank them
            tier 0  validated circular representative (junction PASS)
            tier 1  circular representative with status PENDING
            tier 2  linear contig (incl. circular candidates whose junction FAILed,
                    and circles not selected for classification)
            within a tier: fewer suspicious sites (clip pile-ups / coverage gaps from
            the HiFi re-mapping) -> longer -> assembler priority.
            Non-representative members of a dereplication cluster are left out.

``select``  greedy containment filter over an all-vs-all minimap2 PAF of the pool:
            walking down the ranking, a contig is dropped as redundant when
            >= --max-covered-frac of its length is covered, at >= --min-identity, by
            contigs already kept. Strain variants below --min-identity are kept.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from _common import n50, read_fasta, read_tsv, setup_logging, to_float, write_fasta, write_tsv

POOL_COLUMNS = [
    "rank",
    "contig_id",
    "assembler",
    "length",
    "topology",
    "tier",
    "tier_label",
    "status",
    "category",
    "cluster_id",
    "median_depth",
    "n_suspicious_sites",
]
SELECT_COLUMNS = POOL_COLUMNS + ["kept", "covered_frac", "redundant_with"]
TIER_LABEL = {0: "circular_validated", 1: "circular_pending", 2: "linear"}


# -----------------------------------------------------------------------------
# pool
# -----------------------------------------------------------------------------
def build_pool(
    contig_tables: dict[str, list[dict]],
    scans: dict[str, dict],
    classification: dict[str, dict],
    clusters: dict[str, dict],
    min_len: int,
    include_failed_circular: bool,
    priority: list[str],
) -> list[dict]:
    rows = []
    for asm, contigs in contig_tables.items():
        for c in contigs:
            cid = c["contig_id"]
            cls = classification.get(cid)
            scan = scans.get(cid, {})
            row = {
                "contig_id": cid,
                "assembler": asm,
                "length": int(c["length"]),
                "status": "",
                "category": "",
                "cluster_id": "",
                "median_depth": scan.get("median_depth", ""),
                "n_suspicious_sites": int(to_float(scan.get("n_suspicious_sites"), 0)),
            }
            if cls is not None:  # circular candidate
                row["length"] = int(cls["length"])  # closing overlap already trimmed
                row["status"] = cls["status"]
                row["category"] = cls.get("category", "")
                if cls.get("selected") == "yes":
                    member = clusters.get(cid, {})
                    if member.get("representative") != "yes":
                        continue  # redundant copy of another assembler's circle
                    row["cluster_id"] = member.get("cluster_id", "")
                    row["tier"] = 0 if cls["status"] == "PASS" else 1
                    row["topology"] = "circular"
                    # the junction-aware count is the relevant one for a validated circle
                    row["n_suspicious_sites"] = int(to_float(cls.get("n_suspicious_sites"), 0))
                elif include_failed_circular:
                    row["tier"] = 2
                    row["topology"] = "linear"
                    row["category"] = "circular_unverified"
                else:
                    continue
            else:
                row["tier"] = 2
                row["topology"] = "linear"
            if row["length"] < min_len and row["tier"] == 2:
                continue
            rows.append(row)

    def key(r):
        asm = r["assembler"]
        return (
            r["tier"],
            r["n_suspicious_sites"],
            -r["length"],
            priority.index(asm) if asm in priority else len(priority),
            r["contig_id"],
        )

    rows.sort(key=key)
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    return rows


# -----------------------------------------------------------------------------
# select
# -----------------------------------------------------------------------------
def parse_paf(path: str, min_identity: float, min_block: int):
    """{contig: [(other, start, end)]} for hits with identity >= min_identity."""
    hits: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 12:
                continue
            q, qs, qe, t, ts, te = f[0], int(f[2]), int(f[3]), f[5], int(f[7]), int(f[8])
            matches, block = int(f[9]), int(f[10])
            if q == t or block < min_block or matches / block < min_identity:
                continue
            hits[q].append((t, qs, qe))
            hits[t].append((q, ts, te))
    return hits


def union_length(intervals: list[tuple[int, int]]) -> int:
    total, last = 0, -1
    for s, e in sorted(intervals):
        s = max(s, last)
        if e > s:
            total += e - s
            last = e
    return total


def greedy_select(pool: list[dict], hits, max_covered_frac: float) -> list[dict]:
    kept: set[str] = set()
    for row in sorted(pool, key=lambda r: int(r["rank"])):
        cid, length = row["contig_id"], int(row["length"])
        by_other: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for other, s, e in hits.get(cid, []):
            if other in kept:
                by_other[other].append((s, e))
        all_iv = [iv for ivs in by_other.values() for iv in ivs]
        frac = min(1.0, union_length(all_iv) / length) if length else 0.0
        row["covered_frac"] = round(frac, 4)
        if by_other and frac >= max_covered_frac:
            row["kept"] = False
            row["redundant_with"] = max(by_other, key=lambda o: union_length(by_other[o]))
        else:
            row["kept"] = True
            row["redundant_with"] = ""
            kept.add(cid)
    return pool


def summarise(rows: list[dict]) -> dict:
    kept = [r for r in rows if r["kept"]]
    lengths = [int(r["length"]) for r in kept]
    return {
        "contigs": len(kept),
        "total_bp": sum(lengths),
        "n50": n50(lengths),
        "circular": sum(r["topology"] == "circular" for r in kept),
        "dropped_redundant": sum(not r["kept"] for r in rows),
    }


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def cmd_pool(args, log) -> None:
    contig_tables = {}
    for path in args.contigs:
        rows = read_tsv(path)
        asm = rows[0]["assembler"] if rows else Path(path).parent.name
        contig_tables[asm] = rows
    scans = {r["contig_id"]: r for f in args.scan for r in read_tsv(f)}
    classification = {r["contig_id"]: r for r in read_tsv(args.classification)}
    clusters = {r["contig_id"]: r for r in read_tsv(args.clusters)}
    pool = build_pool(
        contig_tables,
        scans,
        classification,
        clusters,
        args.min_contig_len,
        args.include_failed_circular,
        args.priority,
    )
    for r in pool:
        r["tier_label"] = TIER_LABEL[r["tier"]]
    write_tsv(pool, args.out_tsv, POOL_COLUMNS)

    needed = {r["contig_id"] for r in pool}
    seqs: dict[str, str] = {}
    for f in args.trimmed + args.fasta:  # trimmed (overlap-free) circles take precedence
        for name, _d, seq in read_fasta(f):
            if name in needed and name not in seqs:
                seqs[name] = seq
    write_fasta(((r["contig_id"], seqs[r["contig_id"]]) for r in pool), args.out_fasta)
    log.info(
        "pool: %d contigs (%s)",
        len(pool),
        {t: sum(r["tier"] == t for r in pool) for t in TIER_LABEL},
    )


def cmd_select(args, log) -> None:
    pool = read_tsv(args.pool)
    hits = parse_paf(args.paf, args.min_identity, args.min_block)
    rows = greedy_select(pool, hits, args.max_covered_frac)
    write_tsv(rows, args.out_tsv, SELECT_COLUMNS)

    seqs = {name: seq for name, _d, seq in read_fasta(args.fasta)}

    def records():
        for r in rows:
            if not r["kept"]:
                continue
            header = (
                f"{r['contig_id']} assembler={r['assembler']} length={r['length']} "
                f"topology={r['topology']}"
            )
            if r["topology"] == "circular":
                header += f" status={r['status']} category={r['category']}"
            yield header, seqs[r["contig_id"]]

    write_fasta(records(), args.out_fasta)
    log.info("merged assembly: %s", summarise(rows))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("pool", help="collect and rank contigs from all assemblers")
    p.add_argument("--contigs", nargs="+", required=True, help="03_contigs/*/contigs.tsv")
    p.add_argument("--fasta", nargs="+", required=True, help="03_contigs/*/contigs.fa")
    p.add_argument("--trimmed", nargs="+", required=True, help="04_candidates/*/trimmed.fa")
    p.add_argument("--scan", nargs="+", required=True, help="05_validation/*/contig_scan.tsv")
    p.add_argument("--classification", required=True)
    p.add_argument("--clusters", required=True)
    p.add_argument("--min-contig-len", type=int, default=1000)
    p.add_argument("--include-failed-circular", action="store_true")
    p.add_argument("--priority", nargs="*", default=[])
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-fasta", required=True)

    s = sub.add_parser("select", help="greedy non-redundant selection")
    s.add_argument("--pool", required=True)
    s.add_argument("--fasta", required=True)
    s.add_argument("--paf", required=True)
    s.add_argument("--min-identity", type=float, default=0.99)
    s.add_argument("--max-covered-frac", type=float, default=0.9)
    s.add_argument("--min-block", type=int, default=500)
    s.add_argument("--out-tsv", required=True)
    s.add_argument("--out-fasta", required=True)

    args = ap.parse_args(argv)
    log = setup_logging()
    (cmd_pool if args.cmd == "pool" else cmd_select)(args, log)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Merge candidate metadata from all assemblers and stage inputs for classification.

* candidates.tsv  – one row per circular candidate (all assemblers) with assembly,
                    trimming and junction-validation metrics, plus ``selected`` and
                    ``size_class``
* genomes/        – one FASTA per selected candidate >= --min-chrom-len (CheckM2/GTDB-Tk)
* genomad input   – selected candidates < --min-chrom-len (or all selected with
                    --genomad-on-all)
* selected.fa     – all selected (trimmed, original orientation) sequences
* derep_large.fa  – selected sequences >= --small-len (skani dereplication)
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from _common import read_fasta, read_tsv, setup_logging, write_fasta, write_tsv

CONTIG_COLS = ["assembler", "orig_id", "depth", "source", "repeat", "notes"]
TRIM_COLS = ["orig_length", "overlap_len", "overlap_method", "overlap_identity"]
JUNCTION_COLS = [
    "length",
    "junction_pos",
    "median_depth",
    "depth_at_junction",
    "unique_depth_at_junction",
    "spanning_reads",
    "span_ratio",
    "clipped_at_junction",
    "low_mapq_frac",
    "n_clip_sites",
    "n_low_cov_regions",
    "n_suspicious_sites",
    "status",
    "reason",
]
COLUMNS = (
    ["contig_id", "genome_id"]
    + CONTIG_COLS
    + TRIM_COLS
    + JUNCTION_COLS
    + ["selected", "size_class", "to_checkm2", "to_genomad"]
)


def genome_id(contig_id: str) -> str:
    """File-name safe identifier used for CheckM2 / GTDB-Tk."""
    return re.sub(r"[^A-Za-z0-9_-]", "_", contig_id)


def build_table(contigs, trims, junctions, min_chrom_len, classify_pending, genomad_on_all):
    trim_by = {r["contig_id"]: r for r in trims}
    junc_by = {r["contig_id"]: r for r in junctions}
    rows = []
    for c in contigs:
        if c.get("candidate") != "yes":
            continue
        cid = c["contig_id"]
        t, j = trim_by.get(cid, {}), junc_by.get(cid, {})
        row = {"contig_id": cid, "genome_id": genome_id(cid)}
        row.update({k: c.get(k, "") for k in CONTIG_COLS})
        row.update({k: t.get(k, "") for k in TRIM_COLS})
        row.update({k: j.get(k, "") for k in JUNCTION_COLS})
        if not row["length"]:
            row["length"] = t.get("length") or c.get("length")
        status = row["status"] or "NA"
        row["status"] = status
        selected = status == "PASS" or (status == "PENDING" and classify_pending)
        large = int(row["length"]) >= min_chrom_len
        row["selected"] = selected
        row["size_class"] = "large" if large else "small"
        row["to_checkm2"] = selected and large
        row["to_genomad"] = selected and (not large or genomad_on_all)
        rows.append(row)
    return rows


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--contigs", nargs="+", required=True)
    ap.add_argument("--trim", nargs="+", required=True)
    ap.add_argument("--junction", nargs="+", required=True)
    ap.add_argument("--trimmed", nargs="+", required=True)
    ap.add_argument("--min-chrom-len", type=int, default=500000)
    ap.add_argument("--small-len", type=int, default=20000)
    ap.add_argument("--classify-pending", action="store_true")
    ap.add_argument("--genomad-on-all", action="store_true")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-genomes", required=True)
    ap.add_argument("--out-genomad", required=True)
    ap.add_argument("--out-selected", required=True)
    ap.add_argument("--out-derep-large", required=True)
    args = ap.parse_args(argv)
    log = setup_logging()

    contigs = [r for f in args.contigs for r in read_tsv(f)]
    trims = [r for f in args.trim for r in read_tsv(f)]
    junctions = [r for f in args.junction for r in read_tsv(f)]
    rows = build_table(
        contigs, trims, junctions, args.min_chrom_len, args.classify_pending, args.genomad_on_all
    )
    write_tsv(rows, args.out_tsv, COLUMNS)

    seqs = {name: seq for f in args.trimmed for name, _d, seq in read_fasta(f)}
    by_id = {r["contig_id"]: r for r in rows}
    selected = [cid for cid, r in by_id.items() if r["selected"]]

    outdir = Path(args.out_genomes)
    outdir.mkdir(parents=True, exist_ok=True)
    for cid in selected:
        if by_id[cid]["to_checkm2"]:
            write_fasta([(cid, seqs[cid])], outdir / f"{by_id[cid]['genome_id']}.fa")
    write_fasta(((c, seqs[c]) for c in selected if by_id[c]["to_genomad"]), args.out_genomad)
    write_fasta(((c, seqs[c]) for c in selected), args.out_selected)
    write_fasta(
        ((c, seqs[c]) for c in selected if len(seqs[c]) >= args.small_len), args.out_derep_large
    )
    status_counts = {}
    for r in rows:
        status_counts[r["status"]] = status_counts.get(r["status"], 0) + 1
    log.info(
        "%d candidates %s; %d selected (%d large, %d to geNomad)",
        len(rows),
        status_counts,
        len(selected),
        sum(by_id[c]["to_checkm2"] for c in selected),
        sum(by_id[c]["to_genomad"] for c in selected),
    )


if __name__ == "__main__":
    main()

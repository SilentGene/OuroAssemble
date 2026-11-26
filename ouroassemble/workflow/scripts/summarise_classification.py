#!/usr/bin/env python3
"""Combine CheckM2, GTDB-Tk and geNomad results with the candidate table.

category (selected candidates only)
  chromosome      >= min_chrom_len, completeness >= X and contamination <= Y
  plasmid/virus   geNomad call (small candidates, or large with genomad_on_all)
  large_circular  >= min_chrom_len but not passing the CheckM2 thresholds
  other_circular  small, no geNomad call
unselected candidates get ``unvalidated``.
"""

from __future__ import annotations

import argparse
import gzip

from _common import read_tsv, setup_logging, to_float, write_tsv
from prepare_classification import COLUMNS as CANDIDATE_COLUMNS

EXTRA_COLUMNS = [
    "checkm2_completeness",
    "checkm2_contamination",
    "gtdb_classification",
    "gtdb_closest_reference",
    "gtdb_closest_ani",
    "ref_genome_length",
    "length_ratio",
    "genomad_class",
    "genomad_score",
    "genomad_taxonomy",
    "category",
]
COLUMNS = CANDIDATE_COLUMNS + EXTRA_COLUMNS


def strip_gtdb_prefix(accession: str) -> str:
    return accession[3:] if accession[:3] in ("RS_", "GB_") else accession


def genome_sizes(metadata_files: list[str], wanted: set[str]) -> dict[str, int]:
    """Stream GTDB metadata TSV(.gz) files and return genome_size for wanted accessions."""
    sizes: dict[str, int] = {}
    if not wanted:
        return sizes
    for path in metadata_files:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            i_acc, i_size = header.index("accession"), header.index("genome_size")
            for line in fh:
                parts = line.split("\t", i_size + 1)
                acc = strip_gtdb_prefix(parts[i_acc])
                if acc in wanted:
                    sizes[acc] = int(float(parts[i_size]))
                    if len(sizes) == len(wanted):
                        return sizes
    return sizes


def genomad_calls(plasmid_tsv: str | None, virus_tsv: str | None) -> dict[str, dict]:
    calls: dict[str, dict] = {}
    for kind, path, score_col in (
        ("plasmid", plasmid_tsv, "plasmid_score"),
        ("virus", virus_tsv, "virus_score"),
    ):
        if not path:
            continue
        for row in read_tsv(path):
            name = row["seq_name"].split("|")[0]
            score = to_float(row.get(score_col), 0.0)
            if name not in calls or score > calls[name]["genomad_score"]:
                calls[name] = {
                    "genomad_class": kind,
                    "genomad_score": score,
                    "genomad_taxonomy": row.get("taxonomy", ""),
                }
    return calls


def categorise(row: dict, min_chrom_len: int, min_comp: float, max_cont: float) -> str:
    if row.get("selected") != "yes":
        return "unvalidated"
    large = int(row["length"]) >= min_chrom_len
    comp = to_float(row.get("checkm2_completeness"))
    cont = to_float(row.get("checkm2_contamination"))
    if large and comp is not None and cont is not None and comp >= min_comp and cont <= max_cont:
        return "chromosome"
    if row.get("genomad_class") in ("plasmid", "virus"):
        return row["genomad_class"]
    return "large_circular" if large else "other_circular"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--checkm2")
    ap.add_argument("--gtdbtk")
    ap.add_argument("--plasmid")
    ap.add_argument("--virus")
    ap.add_argument("--gtdb-metadata", nargs="*", default=[])
    ap.add_argument("--min-chrom-len", type=int, default=500000)
    ap.add_argument("--min-completeness", type=float, default=90)
    ap.add_argument("--max-contamination", type=float, default=5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    log = setup_logging()

    rows = read_tsv(args.candidates)
    checkm2 = {r["Name"]: r for r in read_tsv(args.checkm2)} if args.checkm2 else {}
    gtdb = {r["user_genome"]: r for r in read_tsv(args.gtdbtk)} if args.gtdbtk else {}
    genomad = genomad_calls(args.plasmid, args.virus)

    refs = {
        strip_gtdb_prefix(r.get("closest_genome_reference", ""))
        for r in gtdb.values()
        if r.get("closest_genome_reference") not in (None, "", "N/A")
    }
    sizes = genome_sizes(args.gtdb_metadata, refs)
    log.info("reference genome sizes found for %d/%d references", len(sizes), len(refs))

    for row in rows:
        gid = row["genome_id"]
        if gid in checkm2:
            row["checkm2_completeness"] = checkm2[gid].get("Completeness", "")
            row["checkm2_contamination"] = checkm2[gid].get("Contamination", "")
        if gid in gtdb:
            g = gtdb[gid]
            ref = strip_gtdb_prefix(g.get("closest_genome_reference", "") or "")
            row["gtdb_classification"] = g.get("classification", "")
            row["gtdb_closest_reference"] = ref
            row["gtdb_closest_ani"] = g.get("closest_genome_ani", "")
            if ref in sizes:
                row["ref_genome_length"] = sizes[ref]
                row["length_ratio"] = round(int(row["length"]) / sizes[ref], 4)
        row.update(genomad.get(row["contig_id"], {}))
        row["category"] = categorise(
            row, args.min_chrom_len, args.min_completeness, args.max_contamination
        )

    write_tsv(rows, args.out, COLUMNS)
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["category"]] = counts.get(r["category"], 0) + 1
    log.info("categories: %s", counts)


if __name__ == "__main__":
    main()

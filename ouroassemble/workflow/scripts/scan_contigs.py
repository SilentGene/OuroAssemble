#!/usr/bin/env python3
"""Read-support scan of EVERY contig of one assembly (clip pile-ups, coverage gaps).

Uses the per-assembly validation alignments (all contigs of the assembly are in the
reference), so the result ranks linear contigs for the merged assembly the same way
circular candidates are checked.
"""

from __future__ import annotations

import argparse

from _common import setup_logging, write_tsv
from eval_junctions import Params, bam_alignments, scan

COLUMNS = [
    "contig_id",
    "length",
    "median_depth",
    "n_clip_sites",
    "n_low_cov_regions",
    "n_suspicious_sites",
]


def main(argv=None) -> None:
    import pysam

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bam", required=True)
    for name in ("min_clip", "clip_window", "clip_frac", "min_clip_reads", "cov_window",
                 "low_cov_frac", "min_depth"):  # fmt: skip
        default = getattr(Params(), name)
        ap.add_argument("--" + name.replace("_", "-"), type=type(default), default=default)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    log = setup_logging()
    p = Params(**{k: getattr(args, k) for k in vars(Params()) if hasattr(args, k)})

    rows = []
    with pysam.AlignmentFile(args.bam) as bam:
        for contig in bam.references:
            length = bam.get_reference_length(contig)
            rows.append({"contig_id": contig, **scan(bam_alignments(bam, contig), length, p)})
    write_tsv(rows, args.out, COLUMNS)
    log.info(
        "%d contigs scanned, %d with suspicious sites",
        len(rows),
        sum(r["n_suspicious_sites"] > 0 for r in rows),
    )


if __name__ == "__main__":
    main()

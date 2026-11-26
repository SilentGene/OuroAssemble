#!/usr/bin/env python3
"""Evaluate the circular junction of every rotated candidate with HiFi read alignments.

For each candidate the old end->start seam sits at ``junction_pos`` (J). We count

* spanning reads   – alignments with MAPQ >= min_mapq covering [J-flank, J+flank]
                     in one block, without an indel >= min_clip near J
* depth at J       – all alignments (and uniquely placed ones, MAPQ >= min_mapq)
* clipped at J     – alignments with a clip >= min_clip ending within +-clip_window of J
* suspicious sites – elsewhere in the contig: clip pile-ups (>= min_clip_reads and
                     >= clip_frac x local depth) and low-coverage windows
                     (< low_cov_frac x median depth), similar in spirit to
                     anvi-script-find-misassemblies

and assign a status:

PASS     spanning >= min_spanning and span_ratio >= min_span_ratio and no suspicious sites
PENDING  low coverage at J, repeat at J (many low-MAPQ reads), or good junction but
         suspicious sites elsewhere -> keep, needs manual review
FAIL     enough coverage but the junction is not supported by spanning reads
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass, field

import numpy as np
from _common import read_tsv, setup_logging, write_tsv

COLUMNS = [
    "contig_id",
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

EDGE_TOL = 50  # clips this close to a contig end are reads overhanging the end
CLIP_OPS = (4, 5)  # S, H
INDEL_OPS = (1, 2)  # I, D
REF_CONSUMING = (0, 2, 3, 7, 8)  # M D N = X


@dataclass
class Params:
    min_mapq: int = 20
    flank: int = 1000
    min_clip: int = 200
    clip_window: int = 100
    min_spanning: int = 3
    min_span_ratio: float = 0.5
    min_depth: float = 5.0
    repeat_low_mapq_frac: float = 0.5
    clip_frac: float = 0.3
    min_clip_reads: int = 3
    cov_window: int = 1000
    low_cov_frac: float = 0.1


@dataclass
class Aln:
    """Minimal alignment record (decoupled from pysam for testing)."""

    name: str
    ref_start: int
    ref_end: int
    mapq: int
    cigar: list[tuple[int, int]]
    supplementary: bool = False


@dataclass
class Result:
    metrics: dict
    sites: list[tuple[int, int, str]] = field(default_factory=list)


def large_indels(aln: Aln, min_len: int) -> list[int]:
    """Reference positions of insertions/deletions >= min_len."""
    pos, out = aln.ref_start, []
    for op, length in aln.cigar:
        if op in INDEL_OPS and length >= min_len:
            out.append(pos)
        if op in REF_CONSUMING:
            pos += length
    return out


def clip_events(aln: Aln, length: int, min_clip: int) -> list[int]:
    """Reference positions where a read is clipped (or broken by a large indel)."""
    events = []
    if aln.cigar:
        op, n = aln.cigar[0]
        if op in CLIP_OPS and n >= min_clip and aln.ref_start > EDGE_TOL:
            events.append(aln.ref_start)
        op, n = aln.cigar[-1]
        if op in CLIP_OPS and n >= min_clip and aln.ref_end < length - EDGE_TOL:
            events.append(aln.ref_end)
    events.extend(large_indels(aln, min_clip))
    return events


def cluster_positions(events: list[tuple[int, str]], window: int) -> list[tuple[int, int, set]]:
    """Greedy 1-D clustering of (pos, read) events; returns (start, end, reads)."""
    clusters: list[tuple[int, int, set]] = []
    for pos, read in sorted(events):
        if clusters and pos - clusters[-1][1] <= window:
            start, _end, reads = clusters[-1]
            reads.add(read)
            clusters[-1] = (start, pos, reads)
        else:
            clusters.append((pos, pos, {read}))
    return clusters


def clip_pileups(events, cov_all, length: int, p: Params, exclude: int | None = None):
    """Clip pile-up sites; clips within 2 x clip_window of ``exclude`` are returned separately."""
    near: set[str] = set()
    sites = []
    for start, end, reads in cluster_positions(events, p.clip_window):
        centre = (start + end) // 2
        if exclude is not None and abs(centre - exclude) <= 2 * p.clip_window:
            near |= reads
            continue
        local = max(1, int(cov_all[min(centre, length - 1)]))
        if len(reads) >= p.min_clip_reads and len(reads) >= p.clip_frac * local:
            sites.append((start, end + 1, f"clip_pileup:{len(reads)}"))
    return sites, near


def low_coverage_regions(cov_all, median: float, p: Params) -> list[tuple[int, int, str]]:
    """Merged windows whose mean depth is < low_cov_frac x median (if median >= min_depth)."""
    regions: list[tuple[int, int]] = []
    length = len(cov_all)
    if median >= p.min_depth:
        win = max(10, p.cov_window)
        threshold = p.low_cov_frac * median
        region = None
        for s in range(0, length, win):
            e = min(length, s + win)
            if cov_all[s:e].mean() < threshold:
                region = (region[0], e) if region else (s, e)
            elif region:
                regions.append(region)
                region = None
        if region:
            regions.append(region)
    return [(s, e, "low_coverage") for s, e in regions]


def scan(alns: list[Aln], length: int, p: Params) -> dict:
    """Whole-contig read-support scan without a junction (used for every contig)."""
    cov = np.zeros(length + 1, dtype=np.int32)
    events: list[tuple[int, str]] = []
    for a in alns:
        cov[a.ref_start] += 1
        cov[a.ref_end] -= 1
        events.extend((pos, a.name) for pos in clip_events(a, length, p.min_clip))
    cov = np.cumsum(cov)[:length]
    median = float(np.median(cov)) if length else 0.0
    clip_sites, _ = clip_pileups(events, cov, length, p)
    low_sites = low_coverage_regions(cov, median, p)
    return {
        "length": length,
        "median_depth": round(median, 2),
        "n_clip_sites": len(clip_sites),
        "n_low_cov_regions": len(low_sites),
        "n_suspicious_sites": len(clip_sites) + len(low_sites),
    }


def evaluate(alns: list[Aln], length: int, junction: int, p: Params) -> Result:
    cov_all = np.zeros(length + 1, dtype=np.int32)
    cov_hq = np.zeros(length + 1, dtype=np.int32)
    events: list[tuple[int, str]] = []
    spanning: set[str] = set()
    flank = max(0, min(p.flank, junction - 1, length - junction - 1))

    for a in alns:
        cov_all[a.ref_start] += 1
        cov_all[a.ref_end] -= 1
        if a.mapq >= p.min_mapq:
            cov_hq[a.ref_start] += 1
            cov_hq[a.ref_end] -= 1
        events.extend((pos, a.name) for pos in clip_events(a, length, p.min_clip))
        if (
            a.mapq >= p.min_mapq
            and a.ref_start <= junction - flank
            and a.ref_end >= junction + flank
            and not any(abs(x - junction) <= p.clip_window for x in large_indels(a, p.min_clip))
        ):
            spanning.add(a.name)

    cov_all = np.cumsum(cov_all)[:length]
    cov_hq = np.cumsum(cov_hq)[:length]
    j = min(max(junction, 0), length - 1)
    depth_j = int(cov_all[j])
    depth_hq = int(cov_hq[j])
    median = float(np.median(cov_all)) if length else 0.0

    clip_sites, clipped_j = clip_pileups(events, cov_all, length, p, exclude=junction)
    low_sites = low_coverage_regions(cov_all, median, p)

    n_span = len(spanning)
    span_ratio = n_span / depth_hq if depth_hq else 0.0
    low_mapq_frac = 1 - depth_hq / depth_j if depth_j else 0.0
    n_susp = len(clip_sites) + len(low_sites)

    if depth_j < p.min_depth:
        status, reason = "PENDING", "low_coverage_at_junction"
    elif n_span >= p.min_spanning and span_ratio >= p.min_span_ratio:
        if n_susp:
            status, reason = "PENDING", "suspicious_sites_elsewhere"
        else:
            status, reason = "PASS", "junction_spanned"
    elif low_mapq_frac >= p.repeat_low_mapq_frac:
        status, reason = "PENDING", "repeat_at_junction"
    else:
        status, reason = "FAIL", "junction_not_spanned"

    metrics = {
        "length": length,
        "junction_pos": junction,
        "median_depth": round(median, 2),
        "depth_at_junction": depth_j,
        "unique_depth_at_junction": depth_hq,
        "spanning_reads": n_span,
        "span_ratio": round(span_ratio, 3),
        "clipped_at_junction": len(clipped_j),
        "low_mapq_frac": round(low_mapq_frac, 3),
        "n_clip_sites": len(clip_sites),
        "n_low_cov_regions": len(low_sites),
        "n_suspicious_sites": n_susp,
        "status": status,
        "reason": reason,
    }
    return Result(metrics=metrics, sites=clip_sites + low_sites)


def bam_alignments(bam, contig: str) -> list[Aln]:
    out = []
    for r in bam.fetch(contig):
        if r.is_unmapped or r.is_secondary:
            continue
        out.append(
            Aln(
                name=r.query_name,
                ref_start=r.reference_start,
                ref_end=r.reference_end,
                mapq=r.mapping_quality,
                cigar=r.cigartuples or [],
                supplementary=r.is_supplementary,
            )
        )
    return out


def write_igv_batch(path: str, ref: str, bam: str, sites_bed: str, rows, sites) -> None:
    snapdir = os.path.join(os.path.dirname(os.path.abspath(path)), "igv_snapshots")
    with open(path, "w") as fh:
        fh.write("new\n")
        fh.write(f"genome {os.path.abspath(ref)}\n")
        fh.write(f"load {os.path.abspath(bam)}\n")
        fh.write(f"load {os.path.abspath(sites_bed)}\n")
        fh.write(f"snapshotDirectory {snapdir}\nmaxPanelHeight 600\n")
        for row in rows:
            j = int(row["junction_pos"])
            s, e = max(1, j - 3000), j + 3000
            fh.write(f"goto {row['contig_id']}:{s}-{e}\n")
            fh.write(f"snapshot {row['contig_id']}.junction.png\n")
        for contig, start, end, name in sites:
            fh.write(f"goto {contig}:{max(1, start - 3000)}-{end + 3000}\n")
            fh.write(f"snapshot {contig}.{start}.{name.split(':')[0]}.png\n")
        fh.write("exit\n")


def main(argv=None) -> None:
    import pysam

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--candidates", required=True, help="trim_rotate.tsv")
    for name, default in vars(Params()).items():
        ap.add_argument("--" + name.replace("_", "-"), type=type(default), default=default)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-sites", required=True)
    ap.add_argument("--out-igv", required=True)
    ap.add_argument("--out-igv-bed", required=True)
    args = ap.parse_args(argv)
    log = setup_logging()
    p = Params(**{k: getattr(args, k) for k in vars(Params())})

    rows, all_sites = [], []
    with pysam.AlignmentFile(args.bam) as bam:
        for cand in read_tsv(args.candidates):
            contig = cand["contig_id"]
            length = bam.get_reference_length(contig)
            res = evaluate(bam_alignments(bam, contig), length, int(cand["junction_pos"]), p)
            rows.append({"contig_id": contig, **res.metrics})
            all_sites.extend((contig, s, e, n) for s, e, n in res.sites)
            log.info(
                "%s: %s (%s) spanning=%d depth=%d",
                contig,
                res.metrics["status"],
                res.metrics["reason"],
                res.metrics["spanning_reads"],
                res.metrics["depth_at_junction"],
            )

    write_tsv(rows, args.out_tsv, COLUMNS)
    with open(args.out_sites, "w") as fh:
        for contig, s, e, name in all_sites:
            fh.write(f"{contig}\t{s}\t{e}\t{name}\n")
    with open(args.out_igv_bed, "w") as fh:
        for row in rows:
            j = int(row["junction_pos"])
            fh.write(f"{row['contig_id']}\t{max(0, j - 1)}\t{j + 1}\tjunction:{row['status']}\n")
        for contig, s, e, name in all_sites:
            fh.write(f"{contig}\t{s}\t{e}\t{name}\n")
    write_igv_batch(args.out_igv, args.ref, args.bam, args.out_igv_bed, rows, all_sites)
    log.info("evaluated %d candidates, %d suspicious sites", len(rows), len(all_sites))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Trim duplicated end overlaps of circular candidates and rotate them.

1. Look for a copy of the contig start at the contig end (a "closing overlap"
   kept by some assemblers, e.g. HiCanu). An exact suffix/prefix match is tried
   first, then an approximate one with minimap2 (mappy, asm5 preset). If found,
   one copy is removed from the end.
2. Rotate at the midpoint: A + B -> B + A. The original end->start seam is now
   inside the sequence at ``junction_pos`` (0-based: first base of the old start).
3. Build the per-assembly validation reference: every non-candidate contig plus
   the rotated candidates, so reads can still map to background contigs.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

from _common import read_fasta, setup_logging, write_fasta, write_tsv

COLUMNS = [
    "contig_id",
    "orig_length",
    "overlap_len",
    "overlap_method",
    "overlap_identity",
    "length",
    "junction_pos",
]


@dataclass
class Overlap:
    length: int = 0
    method: str = "none"
    identity: float | None = None


def exact_overlap(seq: str, max_len: int) -> int:
    """Length of the longest proper prefix of ``seq`` that is also its suffix (<= max_len)."""
    max_len = min(max_len, len(seq) // 2)
    if max_len <= 0:
        return 0
    s = seq[:max_len] + "\x00" + seq[-max_len:]
    pi = [0] * len(s)
    for i in range(1, len(s)):
        k = pi[i - 1]
        while k and s[i] != s[k]:
            k = pi[k - 1]
        if s[i] == s[k]:
            k += 1
        pi[i] = k
    return pi[-1]


def approx_overlap(seq: str, max_len: int, min_overlap: int, min_identity: float) -> Overlap:
    """Approximate closing overlap: the start of ``seq`` aligned to the end of ``seq``."""
    import mappy  # imported lazily so the exact path works without mappy

    span = min(max_len, len(seq) // 2)
    if span < min_overlap:
        return Overlap()
    prefix, tail = seq[:span], seq[-span:]
    aligner = mappy.Aligner(seq=tail, preset="asm5")
    if not aligner:
        return Overlap()
    tol = 50
    best = Overlap()
    for hit in aligner.map(prefix):
        if hit.strand != 1 or hit.q_st > tol or hit.r_en < len(tail) - tol:
            continue
        identity = hit.mlen / hit.blen if hit.blen else 0.0
        length = len(tail) - hit.r_st + hit.q_st
        if length >= min_overlap and identity >= min_identity and length > best.length:
            best = Overlap(length=length, method="minimap2", identity=round(identity, 5))
    return best


def find_overlap(
    seq: str, max_len: int, min_overlap: int, min_identity: float, max_frac: float = 0.5
) -> Overlap:
    """Closing overlap to trim (at most ``max_frac`` x length).

    Assembler artefacts are near-identical copies of the same consensus (HiCanu even writes
    small circles as ~2 full copies), whereas genuine tandem repeats (e.g. the two trpEG
    copies of Buchnera pTrp, 99.7 %) diverge: ``min_identity`` (default 0.999) separates them.
    """
    limit = int(len(seq) * max_frac)
    k = exact_overlap(seq, min(max_len, limit))
    if k >= min_overlap:
        return Overlap(length=k, method="exact", identity=1.0)
    ov = approx_overlap(seq, min(max_len, limit), min_overlap, min_identity)
    return ov if ov.length <= limit else Overlap()


def rotate(seq: str) -> tuple[str, int]:
    """A + B -> B + A, cut at the midpoint. Returns (rotated, junction_pos)."""
    mid = len(seq) // 2
    return seq[mid:] + seq[:mid], len(seq) - mid


def process(seq: str, max_len: int, min_overlap: int, min_identity: float, max_frac: float = 0.5):
    ov = find_overlap(seq, max_len, min_overlap, min_identity, max_frac)
    trimmed = seq[: len(seq) - ov.length] if ov.length else seq
    rotated, junction = rotate(trimmed)
    return ov, trimmed, rotated, junction


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--overlap-search", type=int, default=50000)
    ap.add_argument("--min-overlap", type=int, default=100)
    ap.add_argument("--min-identity", type=float, default=0.999)
    ap.add_argument("--max-overlap-frac", type=float, default=0.5)
    ap.add_argument("--out-trimmed", required=True)
    ap.add_argument("--out-rotated", required=True)
    ap.add_argument("--out-ref", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-bed", required=True)
    args = ap.parse_args(argv)
    log = setup_logging()

    rows, trimmed_recs, rotated_recs = [], [], []
    rotated_by_id = {}
    for name, _desc, seq in read_fasta(args.candidates):
        ov, trimmed, rotated, junction = process(
            seq, args.overlap_search, args.min_overlap, args.min_identity, args.max_overlap_frac
        )
        if ov.length:
            log.info("%s: trimmed %d bp closing overlap (%s)", name, ov.length, ov.method)
        rows.append(
            {
                "contig_id": name,
                "orig_length": len(seq),
                "overlap_len": ov.length,
                "overlap_method": ov.method,
                "overlap_identity": ov.identity,
                "length": len(trimmed),
                "junction_pos": junction,
            }
        )
        trimmed_recs.append((name, trimmed))
        rotated_recs.append((name, rotated))
        rotated_by_id[name] = rotated

    write_fasta(trimmed_recs, args.out_trimmed)
    write_fasta(rotated_recs, args.out_rotated)
    write_tsv(rows, args.out_tsv, COLUMNS)
    with open(args.out_bed, "w") as bed:
        for r in rows:
            j = r["junction_pos"]
            bed.write(f"{r['contig_id']}\t{max(0, j - 1)}\t{j + 1}\tjunction\n")

    def reference():
        for name, _desc, seq in read_fasta(args.contigs):
            yield name, rotated_by_id.get(name, seq)

    n = write_fasta(reference(), args.out_ref)
    log.info("%d candidates processed; validation reference has %d contigs", len(rows), n)


if __name__ == "__main__":
    main()

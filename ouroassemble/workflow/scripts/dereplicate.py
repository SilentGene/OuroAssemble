#!/usr/bin/env python3
"""Dereplicate circular genomes recovered by several assemblers.

Pairs with ANI >= --min-ani and aligned fraction >= --min-af in BOTH directions are
linked (single linkage). Sequences >= --small-len are compared with skani (edges
pre-computed by ``skani triangle -E``); pairs involving a shorter sequence are
compared here with minimap2 (mappy, asm5) against the doubled target sequence so
that differing rotations do not matter.

Representative ranking (best first)
  1. junction status            PASS > PENDING
  2. CheckM2 quality            completeness - 5 x contamination
  3. spanning reads at junction
  4. fewer suspicious sites
  5. length closest to the cluster median
  6. assembler priority (--priority)
"""

from __future__ import annotations

import argparse
import statistics

from _common import read_fasta, read_tsv, setup_logging, to_float, write_tsv

COLUMNS = [
    "contig_id",
    "cluster_id",
    "cluster_size",
    "representative",
    "representative_id",
    "members",
]
STATUS_RANK = {"PASS": 0, "PENDING": 1}


class UnionFind:
    def __init__(self, items):
        self.parent = {i: i for i in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra

    def groups(self):
        out: dict[str, list[str]] = {}
        for item in self.parent:
            out.setdefault(self.find(item), []).append(item)
        return list(out.values())


def skani_edges(path: str, min_ani: float, min_af: float) -> list[tuple[str, str]]:
    edges = []
    for row in read_tsv(path):
        ani = to_float(row.get("ANI"), 0)
        af_r = to_float(row.get("Align_fraction_ref"), 0)
        af_q = to_float(row.get("Align_fraction_query"), 0)
        if ani >= min_ani and af_r >= min_af and af_q >= min_af:
            edges.append((row["Ref_name"].split()[0], row["Query_name"].split()[0]))
    return edges


class CircularIndex:
    """Lazily built minimap2 indexes of doubled (circular) target sequences."""

    def __init__(self, seqs: dict[str, str]):
        self.seqs = seqs
        self._cache: dict[str, object] = {}

    def aligner(self, name: str):
        if name not in self._cache:
            import mappy

            seq = self.seqs[name]
            self._cache[name] = mappy.Aligner(seq=seq + seq, preset="asm5")
        return self._cache[name]

    def ani_af(self, query: str, target: str) -> tuple[float, float]:
        """(ANI %, aligned fraction of query %) of sequence ``query`` vs circular ``target``."""
        qseq = self.seqs[query]
        aligner = self.aligner(target)
        if not aligner:
            return 0.0, 0.0
        intervals, mlen, blen = [], 0, 0
        for hit in aligner.map(qseq):
            intervals.append((hit.q_st, hit.q_en))
            mlen += hit.mlen
            blen += hit.blen
        if not intervals or not blen:
            return 0.0, 0.0
        covered, last_end = 0, 0
        for st, en in sorted(intervals):
            st = max(st, last_end)
            if en > st:
                covered += en - st
                last_end = en
        return 100 * mlen / blen, 100 * min(covered, len(qseq)) / len(qseq)


def small_edges(seqs: dict[str, str], small_len: int, min_ani: float, min_af: float):
    """Links between pairs where at least one sequence is < small_len."""
    index = CircularIndex(seqs)
    ids = sorted(seqs, key=lambda x: (len(seqs[x]), x))
    edges = []
    for i, a in enumerate(ids):
        la = len(seqs[a])
        for b in ids[i + 1 :]:
            lb = len(seqs[b])
            if la / lb < min_af / 100:
                break  # sorted by length: no longer b can reach the AF threshold
            if la >= small_len and lb >= small_len:
                continue  # handled by skani
            ani1, af1 = index.ani_af(a, b)
            if ani1 < min_ani or af1 < min_af:
                continue
            ani2, af2 = index.ani_af(b, a)
            if ani2 >= min_ani and af2 >= min_af:
                edges.append((a, b))
    return edges


def rank_key(row: dict, median_len: float, priority: list[str]):
    comp = to_float(row.get("checkm2_completeness"), 0.0)
    cont = to_float(row.get("checkm2_contamination"), 0.0)
    asm = row.get("assembler", "")
    return (
        STATUS_RANK.get(row.get("status"), 9),
        -(comp - 5 * cont),
        -int(to_float(row.get("spanning_reads"), 0)),
        int(to_float(row.get("n_suspicious_sites"), 0)),
        abs(int(row["length"]) - median_len),
        priority.index(asm) if asm in priority else len(priority),
        row["contig_id"],
    )


def cluster(rows: list[dict], edges, priority: list[str]) -> list[dict]:
    by_id = {r["contig_id"]: r for r in rows}
    uf = UnionFind(by_id)
    for a, b in edges:
        if a in by_id and b in by_id:
            uf.union(a, b)
    clusters = []
    for members in uf.groups():
        median_len = statistics.median(int(by_id[m]["length"]) for m in members)
        members.sort(key=lambda m: rank_key(by_id[m], median_len, priority))
        clusters.append(members)
    clusters.sort(key=lambda ms: (-int(by_id[ms[0]]["length"]), ms[0]))
    out = []
    for n, members in enumerate(clusters, 1):
        for m in members:
            out.append(
                {
                    "contig_id": m,
                    "cluster_id": f"C{n:04d}",
                    "cluster_size": len(members),
                    "representative": m == members[0],
                    "representative_id": members[0],
                    "members": ",".join(members),
                }
            )
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--classification", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--skani", required=True)
    ap.add_argument("--min-ani", type=float, default=99)
    ap.add_argument("--min-af", type=float, default=90)
    ap.add_argument("--small-len", type=int, default=20000)
    ap.add_argument("--priority", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    log = setup_logging()

    rows = [r for r in read_tsv(args.classification) if r.get("selected") == "yes"]
    seqs = {name: seq for name, _d, seq in read_fasta(args.fasta)}
    edges = skani_edges(args.skani, args.min_ani, args.min_af)
    edges += small_edges(seqs, args.small_len, args.min_ani, args.min_af)
    out = cluster(rows, edges, args.priority)
    write_tsv(out, args.out, COLUMNS)
    n_clusters = len({r["cluster_id"] for r in out})
    log.info("%d genomes -> %d clusters (%d links)", len(rows), n_clusters, len(edges))


if __name__ == "__main__":
    main()

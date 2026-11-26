#!/usr/bin/env python3
"""Simulate PacBio HiFi-like reads from circular reference sequences.

Reads are sampled uniformly around each circle (wrapping over the origin), on a
random strand, with log-normal lengths (capped at the molecule length for small
circles such as plasmids) and ~0.1 % substitution/indel errors (Q30 qualities).

usage: simulate_hifi.py OUT.fq.gz GENOME.fa:DEPTH [GENOME.fa:DEPTH ...]
"""

import gzip
import math
import random
import sys

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name:
        yield name, "".join(chunks).upper()


def mutate(seq, rng, rate=0.001):
    out = []
    for base in seq:
        r = rng.random()
        if r < rate / 3:
            out.append(rng.choice("ACGT".replace(base, "") or "A"))
        elif r < 2 * rate / 3:
            continue  # deletion
        elif r < rate:
            out.append(base + rng.choice("ACGT"))  # insertion
        else:
            out.append(base)
    return "".join(out)


def main():
    out, specs = sys.argv[1], sys.argv[2:]
    rng = random.Random(7)
    n = 0
    with gzip.open(out, "wt") as fq:
        for spec in specs:
            path, depth = spec.rsplit(":", 1)
            for name, seq in read_fasta(path):
                glen = len(seq)
                circle = seq + seq
                mean_len = 12000
                target = float(depth) * glen
                produced = 0
                while produced < target:
                    length = int(rng.lognormvariate(math.log(mean_len), 0.35))
                    length = max(600, min(length, glen))
                    start = rng.randrange(glen)
                    frag = circle[start : start + length]
                    if rng.random() < 0.5:
                        frag = frag.translate(COMP)[::-1]
                    frag = mutate(frag, rng)
                    n += 1
                    fq.write(f"@read{n} {name}\n{frag}\n+\n{'?' * len(frag)}\n")
                    produced += length
                print(f"{name}\t{glen}\t{depth}x", file=sys.stderr)
    print(f"{n} reads written to {out}", file=sys.stderr)


if __name__ == "__main__":
    main()

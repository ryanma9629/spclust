#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
"""
The **fixed generator** for benchmark inputs.

Why it must be fixed: runtime depends strongly on the similarity distribution. With the
same object and edge counts, a range of 0..100 (most edges below the threshold, few
merges) and 1..65535 (almost everything merges) measured over 3x apart. Change the seed
or the range and the figures from before and after are not comparable. So the parameters
are pinned below, and the numbers in tests/bench_baseline.txt are only comparable while
they stay the same; changing them requires re-recording the baseline
(python3 tests/bench.py --update).

Usage:
    python3 tests/gen_bench.py <case> > /tmp/case.csv     # generate one case
    python3 tests/gen_bench.py --list                      # list all cases
    python3 tests/gen_bench.py --all <dir>                 # generate all into a directory
"""

import os
import random
import sys

# (name, objects n, mean degree deg, sim min, sim max, seed)
# edges ~ n * deg / 2. The upper bound is the reader's limit: sim in [0, 65535].
CASES = [
    ("bench_small",  2000, 20,     1, 65535, 1),   # fast, used for many repetitions
    ("bench_mid",   20000, 20,     1, 65535, 2),   # most edges merge (heavy)
    ("bench_large", 50000, 20,     1, 65535, 3),   # larger, to see the scaling trend
    ("bench_lowval", 20000, 20,    0,   100, 4),   # most edges below the threshold (light)
    ("bench_single", 50000, 20,    1, 65535, 5),   # for single linkage (the Kruskal path)
]


def case_params(name):
    for c in CASES:
        if c[0] == name:
            return c
    raise SystemExit(f"unknown case {name!r}; use --list to see the available ones")


def generate(name):
    """Generate edge-list text (objects named o<i>, pairs sorted, so it is reproducible)."""
    _, n, deg, lo, hi, seed = case_params(name)
    rng = random.Random(seed)
    target = n * deg // 2
    pairs = set()
    # Draw distinct pairs through a set: with far fewer edges than n^2 this is fastest
    # and safest.
    while len(pairs) < target:
        a = rng.randrange(n)
        b = rng.randrange(n)
        if a == b:
            continue
        if a > b:
            a, b = b, a
        pairs.add((a, b))
    # Values are drawn after the pairs, so which edge gets which value is deterministic too.
    lines = []
    for i, j in sorted(pairs):
        lines.append(f"o{i},o{j},{rng.randint(lo, hi)}")
    return "\n".join(lines) + "\n"


def main(argv):
    if len(argv) < 2 or argv[1] in ("-h", "--help"):
        print(__doc__)
        return 0
    if argv[1] == "--list":
        for name, n, deg, lo, hi, seed in CASES:
            print(f"{name:<14} n={n:<7} deg={deg:<3} sim in [{lo},{hi}] seed={seed:<3} "
                  f"edges ~{n * deg // 2}")
        return 0
    if argv[1] == "--all":
        outdir = argv[2]
        os.makedirs(outdir, exist_ok=True)
        for name, *_ in CASES:
            path = os.path.join(outdir, name + ".csv")
            with open(path, "w") as f:
                f.write(generate(name))
            print(f"wrote {path} ({os.path.getsize(path)} bytes)")
        return 0
    sys.stdout.write(generate(argv[1]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

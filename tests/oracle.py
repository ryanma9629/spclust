#!/usr/bin/env python3
"""
An independent reference implementation (oracle) for the differential test.

Deliberately written from the *definition* rather than copying the Lance-Williams
recurrence and cross-linked list of the core, so the two sides are computed
independently and must agree exactly.

Model: similarity matrix S, an absent pair counts as 0 (implicit zero). For clusters
A and B the similarity is

    average : sum_{a in A, b in B} S(a,b) / (|A| * |B|)   denominator is all |A|*|B| pairs
    max     : max_{a in A, b in B} S(a,b)                 missing counts as 0
    min     : min_{a in A, b in B} S(a,b)                 a candidate only when > 0

These correspond to spclust's linkage ids 0 / 1 / 2 and all three match the sparseAHC
documentation (average uses the full pair count, min/max ignore missing edges).

Every step merges the globally most similar pair of clusters; if that value is below
min_sim, it stops.

Three modes:
  brute - recompute from the raw dense matrix by definition at every step. Most
          trustworthy, but O(n^4), so only for small inputs.
  fast  - maintain exact aggregates incrementally, then still divide by |A|*|B| by
          definition each step. Suitable for large inputs.
  trunc - mirror spclust's **pre-fix** unsigned-short integer semantics. Not an
          independent oracle but a regression guard: if spclust ever went back to
          integer truncation, the differential test uses this to classify the
          deviation as TRUNC and fail.

The runner first runs brute on every small case and cross-checks that fast agrees.
Everything uses exact Fraction arithmetic (fractional input such as 87.5 is allowed),
so no floating-point error is introduced.

Command line (for debugging):

    python3 tests/oracle.py sim.csv --linkage average --min-sim 50
"""

import json
import re
import sys
from fractions import Fraction

LINKAGES = ("average", "max", "min")

# spclust linkage id -> the linkage name used here
LINKAGE_BY_MTD = {0: "average", 1: "max", 2: "min"}


def to_fraction(s):
    """Parse a string into an exact fraction. Returns None when invalid."""
    try:
        return Fraction(s.strip())  # "90" / "87.5" / "175/2"
    except (ValueError, ZeroDivisionError):
        pass
    try:
        # Tolerate exponent notation like "1e-05"; nan / inf are rejected here
        return Fraction(*float(s).as_integer_ratio())
    except (ValueError, OverflowError):
        return None


def load_edges(path):
    """Read an edge list; returns (order, edges, warnings).

    order holds object names in order of first appearance - the order the CLI assigns
    ids in, which the merge comparison relies on.

    Parsing follows readdb.cpp: '#' comments, comma or whitespace separation, an
    optional header, and self-loops dropped. Duplicate pairs are kept; in the matrix
    """
    with open(path) as f:
        lines = f.read().splitlines()

    order, seen, edges, warnings = [], set(), [], []
    for lineno, raw in enumerate(lines, 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "," in line:
            fields = [x.strip().strip('"') for x in line.split(",")]
        else:
            fields = line.split()
        if len(fields) < 3:
            warnings.append(f"{path}:{lineno}: fewer than 3 fields, skipped")
            continue

        a, b, s = fields[0], fields[1], fields[2]
        if not a or not b:
            warnings.append(f"{path}:{lineno}: empty object name, skipped")
            continue
        v = to_fraction(s)
        if v is None:
            # Only a header at the very beginning is tolerated
            if not edges and [x.lower() for x in fields[:3]] == ["obj1", "obj2", "sim"]:
                continue
            warnings.append(f"{path}:{lineno}: invalid sim value {s!r}, skipped")
            continue
        if not 0 <= v <= 65535:
            warnings.append(f"{path}:{lineno}: sim={v} out of range, skipped")
            continue
        if a == b:
            warnings.append(f"{path}:{lineno}: self-loop {a!r}, skipped")
            continue

        for x in (a, b):
            if x not in seen:
                seen.add(x)
                order.append(x)
        edges.append((a, b, v))

    return order, edges, warnings


def _build_dense(order, edges):
    """The raw similarity matrix (brute mode), missing pairs are 0. Returns (S, duplicate_count).

    When a pair appears several times with different values the **maximum** wins - the
    actual semantics of spclust (its value comes from get_max(), the row maximum).
    """
    n = len(order)
    idx = {o: i for i, o in enumerate(order)}
    S = [[0] * n for _ in range(n)]
    pairs, dup = set(), 0
    for a, b, v in edges:
        i, j = idx[a], idx[b]
        if i == j:
            continue
        if i > j:
            i, j = j, i
        if (i, j) in pairs:
            dup += 1
        pairs.add((i, j))
        if v > S[i][j]:                     # duplicates take the maximum (not last-wins)
            S[i][j] = S[j][i] = v
    return S, dup


def _canon(a, b):
    """Canonicalize a cluster pair into a fixed order for exact comparison."""
    return (a, b) if min(a) < min(b) else (b, a)


def hac(order, edges, linkage="average", min_sim=50, mode="brute"):
    """Run the reference HAC and return a result dict.

    mode='brute' - recompute from the raw matrix by definition each step (O(n^4))
    mode='fast'  - maintain exact aggregates incrementally (O(n^3))
    mode='trunc' - mirror spclust's integer semantics: truncate the link value each step.
                   This is not an independent oracle but a hypothesis test: if spclust's
                   deviation came *entirely* from unsigned-short truncation, it must match
    """
    if linkage not in LINKAGES:
        raise ValueError(f"unknown linkage {linkage!r}, expect one of {LINKAGES}")
    if mode not in ("brute", "fast", "trunc"):
        raise ValueError(f"unknown mode {mode!r}")

    n = len(order)
    S, duplicate_pairs = _build_dense(order, edges)

    members = {i: [i] for i in range(n)}
    live = list(range(n))

    # Accumulators for fast / trunc: (i,j) -> aggregate
    #   fast  : sum for average, max/min for the others
    #   trunc : the truncated integer link value, as spclust stored it
    acc = {}
    if mode in ("fast", "trunc"):
        for i in range(n):
            for j in range(i + 1, n):
                acc[(i, j)] = S[i][j]

    def value(i, j):
        """Exact link similarity between clusters i and j."""
        A, B = members[i], members[j]
        if mode == "brute":
            if linkage == "average":
                # denominator is the full |A|*|B| pair count (missing pairs contribute 0)
                return Fraction(sum(S[a][b] for a in A for b in B), len(A) * len(B))
            if linkage == "max":
                # missing pairs count as 0, so the lower bound is 0
                return Fraction(max(S[a][b] for a in A for b in B), 1)
            # min: missing pairs count as 0, so the result can be 0 (= not a candidate)
            return Fraction(min(S[a][b] for a in A for b in B), 1)

        key = (i, j) if i < j else (j, i)
        a = acc[key]
        if mode == "fast" and linkage == "average":
            # divide by |A|*|B| again each step, by definition
            return Fraction(a, len(members[i]) * len(members[j]))
        return Fraction(a, 1)

    merges = []
    tie_steps, trunc_steps = set(), set()
    step = 0
    next_id = n

    while len(live) > 1:
        # Taking the argmax of (Fraction, i, j) with max() cannot really go wrong, unlike
        # a hand-written cross-multiplication (the easiest place to swap ids and values).
        cand = []
        for x in range(len(live)):
            for y in range(x + 1, len(live)):
                i, j = live[x], live[y]
                v = value(i, j)
                if linkage == "min" and not (v > 0):
                    continue  # under min, only a result > 0 is a candidate
                cand.append((v, i, j))
        if not cand:
            break

        best_val, i, j = max(cand)
        #A similarity <= 0 stops merging - aligned with spclust's get_max() (0 is not a candidate).
        #So min_sim=0 gives the **full merge tree**; min_sim>0 merely stops early (equivalent
        if best_val < min_sim or best_val <= 0:
            break

        best_floor = best_val.numerator // best_val.denominator
        # Tie for the maximum? (the order is unspecified; both sides may legitimately diverge)
        if sum(1 for v, _, _ in cand if v == best_val) > 1:
            tie_steps.add(step)
        # Tie after truncation? (an unsigned short collapses 88.5 and 88.4 into one value)
        if any(v != best_val and v.numerator // v.denominator == best_floor
               for v, _, _ in cand):
            trunc_steps.add(step)

        A, B = frozenset(members[i]), frozenset(members[j])
        merges.append((_canon(A, B), best_val))

        k = next_id
        next_id += 1
        si, sj = len(members[i]), len(members[j])
        members[k] = members[i] + members[j]

        for m in live:
            if m == i or m == j:
                continue
            if mode in ("fast", "trunc"):
                ki = (i, m) if i < m else (m, i)
                kj = (j, m) if j < m else (m, j)
                kk = (k, m) if k < m else (m, k)
                if linkage == "average":
                    if mode == "trunc":
                        # Mirror spclust:
                        #   sim = (sim1*|G1| + sim2*|G2|) / (|G1|+|G2|)   <- integer division
                        # then stored in an unsigned short, so every step truncates
                        acc[kk] = (acc[ki] * si + acc[kj] * sj) // (si + sj)
                    else:
                        acc[kk] = acc[ki] + acc[kj]
                elif linkage == "max":
                    acc[kk] = acc[ki] if acc[ki] > acc[kj] else acc[kj]
                else:
                    acc[kk] = acc[ki] if acc[ki] < acc[kj] else acc[kj]

        live = [m for m in live if m != i and m != j] + [k]
        for m in (i, j):
            for key in [kk for kk in acc if m in kk]:
                del acc[key]
        step += 1

    order_of_names = [order[i] for i in range(n)]
    # The partition uses object indices so it compares directly with parse_spclust()
    partition = frozenset(frozenset(members[i]) for i in live)

    return {
        "objects": n,
        "order": order_of_names,
        "linkage": linkage,
        "mode": mode,
        "min_sim": min_sim,
        "merge_count": len(merges),
        "merges": merges,
        "tie_steps": sorted(tie_steps),
        "trunc_steps": sorted(trunc_steps),
        "duplicate_pairs": duplicate_pairs,
        "groups": sorted(
            (sorted(order_of_names[x] for x in members[i]) for i in live),
            key=lambda g: (-len(g), g),
        ),
        "partition": partition,
    }


def cut(result, min_sim):
    """Cut the full merge tree at min_sim, returning a flat partition (frozenset of

    frozensets of object indices). Similarity is non-increasing along the merge sequence, so
    cutting means taking the prefix with sim >= min_sim and joining the two sub-clusters of
    each merge. This matches spclust (a union-find there; joining member sets is blunter).
    The result passed in must be the **full tree** from min_sim=0.
    """
    live = set(frozenset([i]) for i in range(result["objects"]))
    for (A, B), v in result["merges"]:
        if v < min_sim:
            break                      # end of the prefix (the sequence is non-increasing)
        live.discard(A)
        live.discard(B)
        live.add(A | B)
    return frozenset(live)


# ---------------------------------------------------------------------------
# Parsing spclust output
# ---------------------------------------------------------------------------

RE_GROUP = re.compile(r"^Group (\d+) \(size=(\d+), sim=([0-9.eE+-]+)\): (.*)$")
RE_SINGLE = re.compile(r"^Unclustered \((\d+)\): (.*)$")
RE_OBJECTS = re.compile(r"^Objects: (\d+), groups: (\d+)$")


def parse_spclust(text, order):
    """Parse the CLI output into a structure isomorphic to hac().

    The merge sequence comes from **--matrix** (the public linkage matrix: left right
    height size) using scipy's node ids: 0..n-1 are objects, n+i is merge i, so each node's
    member set can be rebuilt row by row. The flat partition comes from the trailing lines.

    This used to parse the DEBUG=2 merge log, coupling the test to the format of one debug
    line (one changed line could make the runner silently under-count); it now uses the

    public interface. order is the first-appearance order of objects - the order the CLI
    assigns ids in, so leaf id i corresponds to order[i].
    """
    n = len(order)
    members = {i: [i] for i in range(n)}
    merges = []
    groups, singles = [], []
    objects = None

    for line in text.splitlines():
        s2 = line.strip()
        if not s2 or s2.startswith("#"):
            continue

        m = RE_OBJECTS.match(s2)
        if m:
            objects = int(m.group(1))
            continue

        m = RE_GROUP.match(s2)
        if m:
            groups.append(frozenset(order.index(x) for x in m.group(4).split()))
            continue

        m = RE_SINGLE.match(s2)
        if m:
            for x in m.group(2).split():
                singles.append(frozenset([order.index(x)]))
            continue

        f = s2.split()
        if len(f) == 4:
            try:
                left, right = int(f[0]), int(f[1])
            except ValueError:
                continue
            height = to_fraction(f[2])
            if height is None:
                continue
            node = n + len(merges)
            merges.append((_canon(frozenset(members[left]), frozenset(members[right])), height))
            members[node] = members[left] + members[right]

    if objects is None:
        raise ValueError("no 'Objects: N, groups: K' line; the CLI output format may have changed")

    partition = frozenset(groups + singles)
    return {
        "objects": objects,
        "merges": merges,
        "groups": groups,
        "partition": partition,
    }


def parse_spclust_partition(text, order):
    return parse_spclust(text, order)["partition"]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv):
    import argparse

    ap = argparse.ArgumentParser(description="independent reference implementation for spclust")
    ap.add_argument("csv")
    ap.add_argument("--linkage", default="average", choices=LINKAGES)
    ap.add_argument("--min-sim", type=int, default=50)
    ap.add_argument("--mode", default="brute", choices=("brute", "fast", "trunc"))
    ap.add_argument("--cut", type=float, default=None,
                    help="cut the result at this threshold (defaults to --min-sim)")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv[1:])

    order, edges, warnings = load_edges(args.csv)
    for w in warnings:
        print(f"warning: {w}", file=sys.stderr)

    r = hac(order, edges, args.linkage, args.min_sim, args.mode)

    if args.json:
        print(json.dumps(
            {
                "objects": r["objects"],
                "merge_count": r["merge_count"],
                "merges": [[sorted(a), sorted(b), str(v)] for (a, b), v in r["merges"]],
                "tie_steps": r["tie_steps"],
                "trunc_steps": r["trunc_steps"],
                "duplicate_pairs": r["duplicate_pairs"],
                "groups": r["groups"],
            },
            ensure_ascii=False, indent=2))
        return 0

    cut_at = args.min_sim if args.cut is None else args.cut
    part = cut(r, cut_at)
    groups = sorted((sorted(order[i] for i in m) for m in part), key=lambda g: (-len(g), g))
    print(f"oracle({args.linkage}/{args.mode}, full tree of {r['merge_count']} merges, "
          f"cut at {cut_at}): {r['objects']} objects, {len(groups)} groups")
    for g in groups:
        if len(g) >= 2:
            print(f"  size={len(g):<4} {g}")
    singles = sorted(x for g in groups if len(g) == 1 for x in g)
    if singles:
        print(f"  Unclustered ({len(singles)}): {singles}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

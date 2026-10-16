#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
#
# scipy interop test: hand the `spclust_cli --scipy` export to the **real scipy**.

# Why this layer exists: every previous test used an oracle we wrote ourselves. Here the ruler
# is scipy - a completely differently sourced, widely used implementation. It checks two things:
#   1. the format really is valid (scipy's is_valid_linkage / to_tree / fcluster accept it);
#   2. the semantics really match (on cases **without ties**, row for row identical to scipy).

# The key to that correspondence: this library's "an absent edge is similarity 0" is exactly
# "distance = C - similarity" (C = the largest similarity), so average linkage equals scipy's
# method='average' applied to the same dense distance matrix.

# About ties: equal similarities never guaranteed a merge order (the contract says so). Note that
# **distinct input values can still tie** - derived intermediate similarities collide (measured:
# 15_many_merges.csv does exactly that). So the row-by-row comparison runs only when the case and
# linkage are tie-free (asked of the oracle); otherwise only format validity is checked.

# Run: python3 tests/test_scipy.py --bin ./spclust
# Without scipy the whole test skips and returns 0, so CI can run without it.

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
import oracle  # noqa: E402

CASES = sorted(str(p) for p in (HERE / "cases").glob("*.csv")) + [str(REPO / "sim.csv")]

# This library's linkage -> scipy's method. Note the third pair is **deliberately different**:
#   Linkage::MinNonzero takes the minimum only over pairs where both sides have an edge, while
#   scipy's complete reaches the maximum distance (= similarity 0) on missing pairs. On dense
#   input they agree (measured); on sparse input they do not.
# (spclust linkage id, the name used by the oracle, scipy's method name)
MATCHING = [(0, "average", "average"), (1, "max", "single")]
KNOWN_DIFFERENT = (2, "complete")


def parse_matrix(text):
    """Parse the --scipy output (# starts a comment; inf is allowed)."""
    rows = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) != 4:
            continue
        try:
            rows.append([float(x) for x in f])
        except ValueError:
            continue
    return rows


def run_scipy(bin_path, case, linkage, extra=()):
    proc = subprocess.run(
        [str(bin_path), str(case), "--linkage", str(linkage), "--scipy", *extra],
        cwd=str(REPO), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"spclust exited with {proc.returncode}:\n{proc.stderr}")
    return parse_matrix(proc.stdout), proc.stderr


def dense_from_case(case):
    """(object name order, dense similarity matrix, largest similarity). Missing pairs are 0,
    matching this library's semantics."""
    import numpy as np
    order, edges, _ = oracle.load_edges(case)
    n = len(order)
    idx = {o: i for i, o in enumerate(order)}
    S = np.zeros((n, n))
    for (a, b, s) in edges:
        i, j = idx[a], idx[b]
        S[i, j] = S[j, i] = max(S[i, j], s)
    return order, S, (float(S.max()) if n else 0.0)


def n_forest_skip(case, linkage_name, n):
    """A tree that is not complete (fewer than n-1 merges) is not compared row by row -
    the order of bridging rows is arbitrary."""
    order, edges, _ = oracle.load_edges(case)
    res = oracle.hac(order, edges, linkage_name, 0,
                     "brute" if len(order) <= 24 else "fast")
    return len(res["merges"]) != n - 1


def has_ties(case, linkage_name):
    order, edges, _ = oracle.load_edges(case)
    res = oracle.hac(order, edges, linkage_name, 0,
                     "brute" if len(order) <= 24 else "fast")
    return len(res["tie_steps"]) > 0


def main():
    ap = argparse.ArgumentParser(description="scipy interop test")
    ap.add_argument("--bin", default="./spclust")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    bin_path = Path(args.bin)
    if not bin_path.is_absolute():
        bin_path = (Path.cwd() / bin_path).resolve()

    try:
        import numpy as np
        from scipy.cluster.hierarchy import (linkage as scipy_linkage, fcluster,
                                             is_valid_linkage, is_monotonic, to_tree)
    except ImportError as e:
        print(f"skipping the scipy interop test: {e}")
        print("(install numpy/scipy to run it, e.g. python3 -m venv /tmp/v && /tmp/v/bin/pip install numpy scipy)")
        return 0

    print("=" * 78)
    print(f"spclust x scipy interop    {bin_path}    scipy {__import__('scipy').__version__}"
          f"    {len(CASES)} cases")
    print("=" * 78)

    bad = []
    n_exact = n_tie = n_valid = 0

    for case in CASES:
        name = Path(case).name
        order, S, C = dense_from_case(case)
        n = len(order)
        if n < 2:
            continue
        D = C - S
        np.fill_diagonal(D, 0.0)
        cond = D[np.triu_indices(n, 1)]

        for mtd, oracle_name, method in MATCHING:
            ties = has_ties(case, oracle_name)
            rows, stderr = run_scipy(bin_path, case, mtd, ("--scipy-bridge",))
            Z = np.array(rows) if rows else np.empty((0, 4))

            # (1) Format validity: scipy itself must accept it
            if not is_valid_linkage(Z, throw=False, warning=False) or not is_monotonic(Z):
                bad.append(f"{name} linkage={mtd}: scipy considers the linkage matrix invalid")
                continue
            try:
                to_tree(Z)
                fcluster(Z, 2, criterion='maxclust')
            except Exception as e:                                    # pragma: no cover
                bad.append(f"{name} linkage={mtd}: a downstream scipy function failed: {e}")
                continue
            n_valid += 1

            # (2) Semantics: when the tree is **complete** (our merge count == n-1) and there are
            # no ties, the rows must match scipy exactly. A forest (components that never merge,
            # with bridging rows appended by us) is out of scope - the bridge distances are right,
            # but the order in which components merge is arbitrary.
            Z_sc = scipy_linkage(cond, method=method)
            if len(Z) != len(Z_sc):
                bad.append(f"{name} linkage={mtd}: rows {len(Z)} != scipy {len(Z_sc)}")
                continue
            if ties:
                n_tie += 1
                if args.verbose:
                    print(f"  (tie, skipping the row-by-row comparison) {name} linkage={mtd}")
                continue
            if n_forest_skip(case, oracle_name, n):
                if args.verbose:
                    print(f"  (forest, skipping the row-by-row comparison) {name} linkage={mtd}")
                continue
            same = True
            for i in range(len(Z)):
                if {int(Z[i][0]), int(Z[i][1])} != {int(Z_sc[i][0]), int(Z_sc[i][1])} \
                   or abs(Z[i][2] - Z_sc[i][2]) > 1e-9 or int(Z[i][3]) != int(Z_sc[i][3]):
                    bad.append(f"{name} linkage={mtd}: row {i} differs "
                               f"{Z[i].tolist()} vs {Z_sc[i].tolist()}")
                    same = False
                    break
            if same:
                n_exact += 1

    # (3) Forests: an incomplete matrix looks to scipy like "a smaller n" (it cannot detect the
    #     missing rows), so what is checked here is the real guarantee: it must warn, and after
    #     --scipy-bridge the row count is exactly n-1 and the last size is exactly n.
    forest_checked = 0
    for case in CASES:
        order, edges, _ = oracle.load_edges(case)
        n = len(order)
        if n < 2:
            continue
        res = oracle.hac(order, edges, "min", 0, "brute" if n <= 24 else "fast")
        if len(res["merges"]) == n - 1:
            continue                      # complete tree, this check does not apply
        forest_checked += 1
        plain, stderr_plain = run_scipy(bin_path, case, 2, ())
        if "forest" not in stderr_plain:
            bad.append(f"{Path(case).name}: the incomplete tree did not warn on stderr")
        bridged, _ = run_scipy(bin_path, case, 2, ("--scipy-bridge",))
        if len(bridged) != n - 1:
            bad.append(f"{Path(case).name}: after bridging there should be n-1={n-1} rows, got {len(bridged)}")
        elif int(bridged[-1][3]) != n:
            bad.append(f"{Path(case).name}: after bridging the last size should be {n}, "
                       f"got {int(bridged[-1][3])}")
        if len(plain) >= len(bridged):
            bad.append(f"{Path(case).name}: the row count did not change after bridging")
        if forest_checked >= 4:
            break
    if forest_checked == 0:
        print("  (no forest cases, skipping the bridging check)")
    else:
        print(f"  {forest_checked} forest cases: all warned, and bridging gives exactly n-1 rows")
    # (4) Known difference: MinNonzero vs scipy complete.
    #     On a dense matrix they must agree; on a sparse one they may differ (and should, or the
    #     check would be pointless).
    dense_same = dense_total = sparse_same = sparse_total = 0
    for case in CASES:
        order, S, C = dense_from_case(case)
        n = len(order)
        if n < 3:
            continue
        D = C - S
        np.fill_diagonal(D, 0.0)
        cond = D[np.triu_indices(n, 1)]
        rows, _ = run_scipy(bin_path, case, 2, ("--scipy-bridge",))
        Z = np.array(rows) if rows else np.empty((0, 4))
        Z_sc = scipy_linkage(cond, method="complete")
        if Z.shape != Z_sc.shape:
            continue
        same = all({int(Z[i][0]), int(Z[i][1])} == {int(Z_sc[i][0]), int(Z_sc[i][1])}
                   and abs(Z[i][2] - Z_sc[i][2]) <= 1e-9 for i in range(len(Z)))
        # dense = every pair has an edge (no zero in the upper triangle)
        dense = bool((np.triu(S) > 0).sum() == n * (n - 1) // 2)
        if dense:
            dense_total += 1
            dense_same += 1 if same else 0
        else:
            sparse_total += 1
            sparse_same += 1 if same else 0
    if dense_total and dense_same != dense_total:
        bad.append(f"MinNonzero and scipy complete must agree on dense input, "
                   f"but {dense_same}/{dense_total} agree")
    print(f"  dense input, MinNonzero == scipy complete: {dense_same}/{dense_total}")
    print(f"  sparse input, MinNonzero == scipy complete: {sparse_same}/{sparse_total}"
          f" (expected to differ: scipy reaches the maximum distance on missing pairs)")

    print("-" * 78)
    print(f"format valid in {n_valid} combinations; row-for-row identical to scipy in {n_exact};"
          f" skipped because of ties in {n_tie}")
    if bad:
        print("failure detail:")
        for b in bad[:12]:
            print("  ✗", b)
    print(f"RESULT: {'PASS' if not bad else 'FAIL'}")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
#
# Python binding tests. The point of the binding is that a Python user gets the C++
# implementation without hand-writing ctypes, so the tests check the two things such a user
# would notice: the results agree with the same library driven from the command line, and the
# things that are easy to get wrong in a binding (dtypes, strides, lifetimes, the GIL, error
# reporting) actually work.
#
# Skips and returns 0 when numpy or the shared library is missing, so `make test` works
# without either.
#
# Run: python3 tests/test_python.py [--lib /path/to/libspclust.so] [--verbose]

import argparse
import array
import ctypes
import random
import subprocess
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "python"))

try:
    import numpy as np
except ImportError as exc:                                   # pragma: no cover
    print(f"skipping the Python binding tests: {exc}")
    sys.exit(0)

import oracle                                                # noqa: E402

CASES = sorted((HERE / "cases").glob("*.csv")) + [REPO / "sim.csv"]

# The command-line binary to compare against. The Makefile build puts it next to the
# sources, a CMake build puts it in the build directory, so the caller passes --cli
# (the CMake test does) and this is only the default.
CLI = REPO / "spclust"


def load_module(lib):
    import spclust
    spclust._capi._lib = None
    try:
        spclust._capi._lib = spclust._capi.load_library(lib)
    except Exception as exc:                                 # pragma: no cover
        print(f"skipping the Python binding tests: {exc}")
        sys.exit(0)
    return spclust


def edges_of(case):
    """(rows, cols, sims, names) with ids in first-appearance order, exactly like the CLI."""
    order, edges, _ = oracle.load_edges(str(case))
    idx = {}
    rows, cols, sims = [], [], []
    for (a, b, s) in edges:
        for name in (a, b):
            if name not in idx:
                idx[name] = len(idx)
        rows.append(idx[a]); cols.append(idx[b]); sims.append(float(s))
    names = [None] * len(idx)
    for name, i in idx.items():
        names[i] = name
    return np.array(rows, dtype=np.uint32), np.array(cols, dtype=np.uint32), \
        np.array(sims, dtype=np.float64), names


def cli_matrix(case, linkage):
    out = subprocess.run([str(CLI), str(case), "--floor", "0", "--matrix",
                          "--linkage", str(linkage)], capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        f = line.split()
        if len(f) == 4:
            try:
                rows.append(tuple(float(x) for x in f))
            except ValueError:
                pass
    return rows


def same_tree(got, want):
    """Compare a binding tree against the CLI's --matrix rows.

    The CLI prints 12 significant digits, so heights are compared with a relative
    tolerance - the ids and sizes, which are what the structure is, exactly.
    """
    if len(got) != len(want):
        return False
    for g, w in zip(got, want):
        if int(g[0]) != int(w[0]) or int(g[1]) != int(w[1]) or int(g[3]) != int(w[3]):
            return False
        if abs(g[2] - w[2]) > 1e-9 * max(1.0, abs(w[2])):
            return False
    return True


def cli_groups(case, cut, linkage):
    out = subprocess.run([str(CLI), str(case), "--cut", str(cut),
                          "--linkage", str(linkage)], capture_output=True, text=True).stdout
    groups = []
    for line in out.splitlines():
        if line.startswith("Group ") and "):" in line:
            groups.append(frozenset(line.split("):", 1)[1].split()))
        elif line.startswith("Unclustered ("):
            groups += [frozenset([x]) for x in line.split("):", 1)[1].split()]
    return frozenset(groups)


def main():
    global CLI
    ap = argparse.ArgumentParser(description="Python binding tests")
    ap.add_argument("--lib", default=None, help="shared library (default: auto-locate)")
    ap.add_argument("--cli", default=str(CLI),
                    help="spclust command-line binary to cross-check against")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    # Path("./spclust") normalizes to "spclust" (dropping the slash) and subprocess would
    # then search PATH, so make it absolute before anything else - the same trap test_cli.py
    # documents.
    CLI = Path(args.cli)
    if not CLI.is_absolute():
        CLI = (Path.cwd() / CLI).resolve()
    if not CLI.exists():
        # Failing loudly matters here: it would otherwise look like a passing run that
        # checked nothing, which is exactly how the CMake build's different binary
        # location stayed hidden until CI.
        print(f"cannot find the spclust CLI at {CLI}; pass --cli <path>")
        return 1

    sp = load_module(args.lib)
    ok = 0
    bad = []

    def check(cond, what):
        nonlocal ok
        if cond:
            ok += 1
        else:
            bad.append(what)
            print(f"  X {what}")

    print("=" * 78)
    print(f"spclust Python binding    {sp.library_path()}    version {sp.version()}")
    print("=" * 78)

    # ---- 1. agreement with the CLI, tree and partition, all linkages and cases ----
    n_cmp = 0
    for linkage, name in ((0, "average"), (1, "max-similarity"), (2, "min-nonzero")):
        for case in CASES:
            rows, cols, sims, names = edges_of(case)
            h = sp.HAC(linkage=name)
            h.add_edges(rows, cols, sims)
            h.build()
            want = cli_matrix(case, linkage)
            got = h.tree()
            if not same_tree(got, want):
                check(False, f"{case.name} linkage={name}: tree differs from the CLI")
                continue
            n_cmp += 1
            # the flat partition, by name, must match too
            cut = 50.0
            c = h.cut(cut)
            got_groups = frozenset(frozenset(names[int(o)] for o in g) for g in c.groups)
            if got_groups != cli_groups(case, cut, linkage):
                check(False, f"{case.name} linkage={name}: partition differs from the CLI")
                continue
            n_cmp += 1
    check(n_cmp == 3 * len(CASES) * 2, f"agreement with the CLI on {n_cmp} comparisons")

    # ---- 2. batch feeding == per-edge feeding ----
    rows, cols, sims, _ = edges_of(CASES[0])
    a, b = sp.HAC(), sp.HAC()
    a.add_edges(rows, cols, sims)
    for i in range(len(rows)):
        b.add_edge(int(rows[i]), int(cols[i]), float(sims[i]))
    a.build(); b.build()
    check(np.array_equal(a.tree(), b.tree()), "add_edges matches add_edge in a loop")

    # ---- 3. reproducibilty: shuffled edges, same tree ----
    perm = np.arange(len(rows))
    random.Random(1).shuffle(perm)
    c = sp.HAC()
    c.add_edges(rows[perm], cols[perm], sims[perm])
    c.build()
    check(np.array_equal(a.tree(), c.tree()), "shuffled edges give the same tree")

    # ---- 4. input flexibility: lists, array.array, memoryview, dtypes, strides ----
    variants = {
        "python lists": (list(map(int, rows)), list(map(int, cols)), list(map(float, sims))),
        "array.array": (array.array("I", map(int, rows)), array.array("I", map(int, cols)),
                        array.array("d", map(float, sims))),
        "memoryview": (memoryview(array.array("I", map(int, rows))),
                       memoryview(array.array("I", map(int, cols))),
                       memoryview(array.array("d", map(float, sims)))),
        "int64 + float32 numpy": (rows.astype(np.int64), cols.astype(np.int64),
                                  sims.astype(np.float32)),
        "int32 numpy": (rows.astype(np.int32), cols.astype(np.int32), sims),
    }
    for what, (ra, rb, rs) in variants.items():
        h = sp.HAC()
        h.add_edges(ra, rb, rs)
        h.build()
        check(np.array_equal(h.tree(), a.tree()), f"accepts {what}")
    # non-contiguous (strided) numpy arrays
    big_r = np.zeros(len(rows) * 2, dtype=np.uint32); big_r[::2] = rows
    big_c = np.zeros(len(cols) * 2, dtype=np.uint32); big_c[::2] = cols
    big_s = np.zeros(len(sims) * 2, dtype=np.float64); big_s[::2] = sims
    h = sp.HAC(); h.add_edges(big_r[::2], big_c[::2], big_s[::2]); h.build()
    check(np.array_equal(h.tree(), a.tree()), "accepts non-contiguous arrays")

    # ---- 5. cuts: labels_ invariants, all three modes, agreement with --matrix semantics ----
    h = sp.HAC()
    h.add_edges(rows, cols, sims)
    h.build()
    n = h.num_objects
    for cut in (0.0, 50.0, 1e9):
        cu = h.cut(cut)
        check(len(cu.labels_) == n, f"labels_ has one entry per object (cut={cut})")
        check(sorted(int(x) for x in cu.labels_) == sorted(range(n)) or
              len(set(int(x) for x in cu.labels_)) == cu.n_groups,
              f"labels_ covers every object (cut={cut})")
        check(sum(cu.sizes()) == n, f"cluster sizes sum to the object count (cut={cut})")
        check(all(cu.cohesion[g] == cu.cohesion[g] for g in range(cu.n_groups)),
              f"cohesion is defined for every cluster (cut={cut})")
    one = h.cut_n_clusters(1)
    check(one.n_groups == 1 and one.sizes() == [n], "cut_n_clusters(1) is one cluster")
    big = sp.HAC()
    br, bc, bs, _ = edges_of(REPO / "sim.csv")          # 16 objects, a full tree
    big.add_edges(br, bc, bs)
    big.build()
    check(big.cut_n_clusters(3).n_groups == 3, "cut_n_clusters(3) gives 3 clusters")
    check(big.cut_n_clusters(big.num_objects).n_groups == big.num_objects,
          "cut_n_clusters(n) leaves everything separate")
    m2 = h.cut_min_size(2)
    check(all(s >= 2 for s in m2.sizes()), "cut_min_size(2) has no cluster below 2")

    # ---- 6. the export goes through the binding unchanged ----
    # Semantics against scipy are already covered by tests/test_scipy.py; what this layer
    # has to show is that the binding passes the matrix through unchanged, so it is compared
    # against the CLI's own --scipy output rather than against scipy.
    n_exp = 0
    for linkage, name in ((0, "average"), (1, "max-similarity")):
        for case in CASES:
            r2, c2, s2, _ = edges_of(case)
            h2 = sp.HAC(linkage=name)
            h2.add_edges(r2, c2, s2)
            h2.build()
            got = h2.to_scipy(bridge=True)
            cli = subprocess.run([str(CLI), str(case), "--scipy", "--scipy-bridge",
                                  "--linkage", str(linkage)],
                                 capture_output=True, text=True).stdout
            want = []
            for line in cli.splitlines():
                f = line.split()
                if len(f) == 4:
                    try:
                        want.append([float(x) for x in f])
                    except ValueError:
                        pass
            if len(want) != len(got) or any(
                    int(g[0]) != int(w[0]) or int(g[1]) != int(w[1]) or
                    int(g[3]) != int(w[3]) or
                    not (abs(g[2] - w[2]) <= 1e-9 * max(1.0, abs(w[2])) or
                         (g[2] != g[2] and w[2] != w[2]))          # both NaN/inf
                    for g, w in zip(got, want)):
                check(False, f"{case.name} linkage={name}: to_scipy differs from the CLI --scipy")
                continue
            n_exp += 1
    check(n_exp == 2 * len(CASES), f"to_scipy matches the CLI on {n_exp} matrices")

    # and scipy itself still accepts what the binding produces
    try:
        from scipy.cluster.hierarchy import fcluster, is_valid_linkage
        # Use a case whose tree is complete (every object ends in one cluster). A forest -
        # sim.csv is one - cannot be compared this way: the bridging rows the export appends
        # are not merges this library ever made, so scipy reaches one cluster while
        # cut_n_clusters(1) legitimately stops at the components.
        chosen = None
        for case in CASES:
            ra, rb, rs, _ = edges_of(case)
            probe = sp.HAC()
            probe.add_edges(ra, rb, rs)
            probe.build()
            if probe.cut_n_clusters(1).n_groups == 1:
                chosen = (case, ra, rb, rs)
                break
        check(chosen is not None, "found a case with a complete tree")
        case, r2, c2, s2 = chosen
        h2 = sp.HAC()
        h2.add_edges(r2, c2, s2)
        h2.build()
        Z = h2.to_scipy(bridge=True)
        check(is_valid_linkage(Z, throw=False, warning=False),
              "scipy accepts a matrix produced through the binding")
        # fcluster works in distances; cutting by cluster count sidesteps the mapping.
        for ncl in (1, 2, 3):
            labels = h2.cut_n_clusters(ncl).labels_
            sc = fcluster(Z, ncl, criterion="maxclust")
            pa = {frozenset(np.where(labels == g)[0]) for g in set(labels.tolist())}
            pb = {frozenset(np.where(sc == g)[0]) for g in set(sc.tolist())}
            if pa != pb:
                check(False, f"fcluster(maxclust={ncl}) disagrees with cut_n_clusters({ncl})")
                break
        else:
            check(True, "cut_n_clusters agrees with scipy fcluster(maxclust)")
    except ImportError:
        print("  (scipy not installed: check that scipy accepts the export)")

    # ---- 7. error reporting ----
    h = sp.HAC()
    def raises(fn, needle, kinds=(sp.SpclustError,)):
        """True when fn() raises one of `kinds` mentioning `needle`.

        The binding draws the line where Python does: arguments it validates itself raise
        ValueError, while anything the library rejects raises SpclustError.
        """
        try:
            fn()
        except kinds as exc:
            return needle in str(exc)
        except Exception:
            return False
        return False
    check(raises(lambda: h.add_edge(3, 3, 5.0), "self-loop"), "self-loop reports itself")
    check(raises(lambda: h.add_edge(0, 1, -1.0), "non-negative"),
          "a negative single edge reports the reason")
    check(raises(lambda: h.add_edges(np.array([0, 1]), np.array([1, 2]),
                                     np.array([1.0, -2.0])), "index 1"),
          "a bad edge inside a batch names its index")
    check(raises(lambda: h.add_edges(np.array([0, 1]), np.array([1, 2]),
                                     np.array([np.nan, 1.0])), "index 0"),
          "NaN inside a batch names its index")
    check(raises(lambda: h.add_edges(np.array([0, 1]), np.array([1]),
                                     np.array([1.0, 2.0])), "length mismatch", (ValueError,)),
          "length mismatch is reported")
    check(raises(lambda: sp.HAC(linkage="nope"), "unknown linkage", (ValueError,)),
          "bad linkage name")
    check(raises(lambda: a.to_scipy(distance="nope"), "unknown distance", (ValueError,)),
          "bad distance mapping")

    # ---- 8. lifetimes and independence ----
    x, y = sp.HAC(), sp.HAC()
    x.add_edges(rows, cols, sims); y.add_edges(rows, cols, sims)
    x.build(); y.build()
    del y
    check(np.array_equal(x.tree(), a.tree()), "one handle can be freed while another lives")
    x.close(); x.close()                                   # idempotent
    check(True, "close() is idempotent")
    with sp.HAC() as z:
        z.add_edges(rows, cols, sims)
        z.build()
        check(z.num_objects == len({*map(int, rows), *map(int, cols)}), "context manager works")
    r = sp.HAC()
    r.add_edges(rows, cols, sims); r.build()
    first = r.tree()
    r.reset()
    check(r.num_objects == 0, "reset clears the object count")
    r.add_edges(rows, cols, sims); r.build()
    check(np.array_equal(r.tree(), first), "reset then reuse gives the same tree")

    # ---- 9. the GIL is released while building ----
    import gen_bench
    text = gen_bench.generate("bench_mid")            # n=20000, 200k edges, ~0.1 s
    pairs = [line.split(",") for line in text.strip().splitlines()]
    name_id = {}
    g_a, g_b, g_s = [], [], []
    for x, y, v in pairs:
        for nm in (x, y):
            name_id.setdefault(nm, len(name_id))
        g_a.append(name_id[x]); g_b.append(name_id[y]); g_s.append(float(v))
    big_a = np.array(g_a, dtype=np.uint32)
    big_b = np.array(g_b, dtype=np.uint32)
    big_s = np.array(g_s, dtype=np.float64)
    ticks = [0]
    stop = threading.Event()

    def ticker():
        while not stop.is_set():
            ticks[0] += 1
            time.sleep(0.001)

    th = threading.Thread(target=ticker, daemon=True)
    th.start()
    t0 = time.perf_counter()
    hb = sp.HAC(linkage="average")
    hb.add_edges(big_a, big_b, big_s)
    hb.build()
    dt = time.perf_counter() - t0
    stop.set(); th.join()
    if dt > 0.05:
        check(ticks[0] > 5, f"another thread made progress during build ({ticks[0]} ticks in {dt*1000:.0f} ms)")
    else:
        print(f"  (build was too fast ({dt*1000:.0f} ms) to test GIL release)")

    print("-" * 78)
    print(f"PASS={ok} FAIL={len(bad)}")
    if bad:
        print("failures:")
        for b in bad[:10]:
            print("  -", b)
    print(f"RESULT: {'PASS' if not bad else 'FAIL'}")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())

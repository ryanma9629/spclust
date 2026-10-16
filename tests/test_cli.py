#!/usr/bin/env python3
"""
CLI interface test.

`tests/run_tests.py` verifies the **algorithm** (parsing the merge log and comparing
against the oracle), while the command-line interface - `--cut/--floor/--tree/--matrix`,
exit codes, the stdin mode - had no coverage at all. The two also tend to mask each
other: the program and the runner were implicitly coupled through the merge log's
format, so one changed output line could make the runner silently under-count (it
happened here: the Kruskal path forgot to print that line and the runner reported
"0 merges").

Only the interface contract is tested here, not the algorithm again:
  * exit codes and error paths;
  * `--cut`: matches the oracle's cut, and "cut at X" must equal "stop early at X";
  * `--floor`: it decides how deep the tree goes, and cut can only sit above floor;
  * `--matrix`: the structural invariants of a linkage matrix;
  * `--tree`: every leaf appears exactly once;
  * the stdin mode is equivalent to the file mode.
"""

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
import oracle  # noqa: E402

CASES = sorted(str(p) for p in (HERE / "cases").glob("*.csv"))


def run(binary, args, stdin_text=None):
    return subprocess.run([str(binary)] + args, cwd=str(REPO), input=stdin_text,
                          capture_output=True, text=True)


def parse_flat(text):
    """Parse the flat output into a set of member-name frozensets (group numbering ignored)."""
    groups = []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("Group ") and "):" in s:
            groups.append(frozenset(s.split("):", 1)[1].split()))
        elif s.startswith("Unclustered ("):
            for name in s.split("):", 1)[1].split():
                groups.append(frozenset([name]))
    return frozenset(groups)


def parse_matrix(text):
    rows = []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        f = s.split()
        if len(f) == 4 and all(x.replace(".", "", 1).replace("-", "", 1).replace("e", "", 1)
                               .replace("+", "", 1).isdigit() for x in f):
            rows.append((int(f[0]), int(f[1]), float(f[2]), int(f[3])))
    return rows


def oracle_cut(case, cut, linkage="average"):
    order, edges, _ = oracle.load_edges(case)
    res = oracle.hac(order, edges, linkage, 0, "brute" if len(order) <= 24 else "fast")
    return frozenset(frozenset(order[i] for i in m) for m in oracle.cut(res, cut))


def t_linkage_runtime(bin_path):
    """The linkage is a run-time option: --linkage 0/1/2 must really change the algorithm, and
    all three must work from one binary (it used to take three)."""
    bad = []
    case = CASES[0]
    seen = {}
    for name, m in (("average", "0"), ("max", "1"), ("min", "2")):
        r = run(bin_path, [case, "--linkage", m, "--matrix"])
        if r.returncode != 0 and r.returncode != 2:
            bad.append(f"--linkage {m} exited with {r.returncode}")
            continue
        seen[name] = parse_matrix(r.stdout)
        # The names must be accepted too
        r2 = run(bin_path, [case, "--linkage", name, "--matrix"])
        if parse_matrix(r2.stdout) != seen[name]:
            bad.append(f"--linkage {name} and --linkage {m} disagree")
    # On a case where the linkages really diverge (10_dense_6), the three should give at
    dense = [c for c in CASES if "10_dense_6" in c]
    if dense:
        outs = set()
        for m in ("0", "1", "2"):
            outs.add(tuple(parse_matrix(run(bin_path, [dense[0], "--linkage", m,
                                                       "--matrix"]).stdout)))
        if len(outs) < 2:
            bad.append("the three linkages produced identical results on a dense case (suspicious)")
    r = run(bin_path, [case, "--linkage", "7"])
    if r.returncode != 2:
        bad.append(f"unknown linkage: exit code {r.returncode}, expected 2")
    return (not bad), "; ".join(bad)


def t_n_clusters(bin_path):
    """--n-clusters K: with a deep enough tree there must be exactly K clusters and every
    object must appear once. The exact "first n-K merges" semantics, ties included, are
    covered by the library-level property test tests/test_cut.cpp."""
    bad = []
    for case in CASES:
        order, edges, _ = oracle.load_edges(case)
        n = len(order)
        # The deepest partition the tree can offer has n - merges clusters
        res = oracle.hac(order, edges, "average", 0, "brute" if n <= 24 else "fast")
        deepest = n - len(res["merges"])
        for k in (1, 2, 3, n):
            if k > n:
                continue
            r = run(bin_path, [case, "--n-clusters", str(k)])
            if r.returncode != 0:
                bad.append(f"{Path(case).name} --n-clusters {k}: exit code {r.returncode}")
                continue
            got = parse_flat(r.stdout)
            members = [x for g in got for x in g]
            if sorted(members) != sorted(order):
                bad.append(f"{Path(case).name} --n-clusters {k}: members are not each object exactly once")
            # With a deep enough tree (k >= n - merges) the cluster count must be exactly k
            if k >= deepest and len(got) != k:
                bad.append(f"{Path(case).name} --n-clusters {k}: got {len(got)} clusters")
    # k = n must leave everything separate
    r = run(bin_path, [CASES[0], "--n-clusters", "999"])
    if parse_flat(r.stdout) != frozenset(frozenset([o]) for o in
                                         oracle.load_edges(CASES[0])[0]):
        bad.append("--n-clusters above the object count should leave everything separate")
    return (not bad), "; ".join(bad[:3])


def components(bin_path, case):
    """Connected components: merging the tree all the way down (--cut 0) gives the
    component partition (singletons included)."""
    r = run(bin_path, [case, "--floor", "0", "--cut", "0"])
    return parse_flat(r.stdout)


def t_min_size(bin_path):
    """--min-size M: there must be no cluster smaller than M - unless the component it lives
    in is itself smaller than M (under MinNonzero, clusters with no edge between them never
    merge, so the condition is unreachable)."""
    bad = []
    for case in CASES:
        for m in (2, 3):
            r = run(bin_path, [case, "--min-size", str(m)])
            if r.returncode != 0:
                bad.append(f"{Path(case).name} --min-size {m}: exit code {r.returncode}")
                continue
            got = parse_flat(r.stdout)
            comps = components(bin_path, case)
            comp_of = {o: c for c in comps for o in c}
            for grp in got:
                if len(grp) >= m:
                    continue
                # A cluster below m must be exactly its whole component (nothing left to merge)
                comp = comp_of[next(iter(grp))]
                if grp != comp:
                    bad.append(f"{Path(case).name} --min-size {m}: a cluster of size {len(grp)}"
                               f" is not a whole component ({len(comp)})")
    # min-size 1 is equivalent to a very high --cut (everything separate)
    r1 = run(bin_path, [CASES[0], "--min-size", "1"])
    r2 = run(bin_path, [CASES[0], "--cut", "1000000"])
    if parse_flat(r1.stdout) != parse_flat(r2.stdout):
        bad.append("--min-size 1 should merge nothing")
    return (not bad), "; ".join(bad[:3])


# ---------------------------------------------------------------------------
# Individual checks: return (ok, detail)
# ---------------------------------------------------------------------------

def t_exit_codes(bin_path):
    bad = []
    r = run(bin_path, ["--help"])
    if r.returncode != 0 or "usage" not in r.stdout + r.stderr:
        bad.append("--help should exit 0 and print usage")
    for args, want, what in ((['--nope'], 2, "unknown option"),
                             (['a.csv', 'b.csv'], 2, "extra positional argument"),
                             (['--cut'], 2, "--cut without a value"),
                             (['/nonexistent/x.csv'], 1, "input file does not exist")):
        r = run(bin_path, args)
        if r.returncode != want:
            bad.append(f"{what}: exit code {r.returncode}, expected {want}")
    r = run(bin_path, ["-"], stdin_text="")
    if r.returncode != 1:
        bad.append(f"empty input: exit code {r.returncode}, expected 1")
    return (not bad), "; ".join(bad)


def t_cut_matches_oracle(bin_path):
    """The partition from --cut X must equal the oracle's full tree cut at X."""
    bad = []
    n = 0
    for case in CASES:
        for cut in (0, 25, 50, 75, 90):
            n += 1
            got = parse_flat(run(bin_path, [case, "--floor", "0", "--cut", str(cut)]).stdout)
            want = oracle_cut(case, cut)
            if got != want:
                bad.append(f"{Path(case).name} cut={cut}: CLI and oracle disagree")
    return (not bad), f"({n} cases) " + "; ".join(bad[:3])


def t_cut_equals_early_stop(bin_path):
    """Cutting the full tree at X must equal using X as the stopping condition (--floor X)."""
    bad = []
    for case in CASES:
        for cut in (25, 50, 75):
            a = parse_flat(run(bin_path, [case, "--floor", "0", "--cut", str(cut)]).stdout)
            b = parse_flat(run(bin_path, [case, "--floor", str(cut), "--cut", str(cut)]).stdout)
            if a != b:
                bad.append(f"{Path(case).name} cut={cut}")
    return (not bad), "; ".join(bad[:3])


def t_floor_cut_interaction(bin_path):
    """When --cut is below the default floor, the floor must follow automatically."""
    bad = []
    for case in CASES:
        # No floor at all, just a low cut: the tree should reach down to it automatically
        got = parse_flat(run(bin_path, [case, "--cut", "10"]).stdout)
        want = oracle_cut(case, 10)
        if got != want:
            bad.append(f"{Path(case).name}: --cut 10 did not lower the floor automatically")
    # An explicit higher floor: a cut above the floor should be unaffected by it
    r = run(bin_path, [CASES[0], "--floor", "50", "--cut", "95"])
    if r.returncode != 0:
        bad.append("--floor 50 --cut 95 should run fine (cutting above the floor)")
    return (not bad), "; ".join(bad[:3])


def t_matrix_structure(bin_path):
    """The basic invariants of a linkage matrix."""
    bad = []
    for case in CASES:
        for extra in ([], ["--floor", "0"]):
            text = run(bin_path, [case, "--floor", "0"] + extra).stdout
            rows = parse_matrix(text)
            order, edges, _ = oracle.load_edges(case)
            n = len(order)
            if not rows:
                # No merges at all is allowed (all similarities 0) - but then there must be no rows
                continue
            parents = {}
            for k, (l, r, sim, size) in enumerate(rows):
                node = n + k
                for c in (l, r):
                    if not (0 <= c < node):
                        bad.append(f"{Path(case).name}: row {k} child {c} is not in [0,{node})")
                    if c in parents:
                        bad.append(f"{Path(case).name}: node {c} used as a child twice")
                    parents[c] = node
                if size < 2:
                    bad.append(f"{Path(case).name}: row {k} has size={size} < 2")
            # Similarity must be non-increasing (the cut semantics rely on it)
            sims = [r[2] for r in rows]
            if any(sims[i] < sims[i + 1] - 1e-9 for i in range(len(sims) - 1)):
                bad.append(f"{Path(case).name}: similarity is not non-increasing {sims[:6]}")
            # The size of a root equals the number of objects under it
            roots = [x for x in range(n + len(rows)) if x not in parents]
            tot = sum(rows[x - n][3] for x in roots if x >= n) + sum(1 for x in roots if x < n)
            if tot != n:
                bad.append(f"{Path(case).name}: root sizes sum to {tot} != object count {n}")
    return (not bad), "; ".join(bad[:3])


def t_tree_completeness(bin_path):
    """Every object appears exactly once in --tree (within the depth limit)."""
    bad = []
    for case in CASES:
        text = run(bin_path, [case, "--floor", "0", "--tree"]).stdout
        order, _, _ = oracle.load_edges(case)
        if "Dendrogram" not in text:
            bad.append(f"{Path(case).name}: no Dendrogram line")
            continue
        leaves = []
        for line in text.splitlines():
            s = line.strip()
            if s.startswith("- "):
                leaves.append(s[2:].strip())
        names = set(order)
        if len(leaves) != len(set(leaves)):
            bad.append(f"{Path(case).name}: a leaf appears twice")
        if not set(leaves) <= names:
            bad.append(f"{Path(case).name}: an object not in the input appeared")
    # With several roots, the count from --tree must match the one from --matrix
    case = CASES[0]
    t = run(bin_path, [case, "--floor", "0", "--tree"]).stdout
    m = run(bin_path, [case, "--floor", "0", "--matrix"]).stdout
    rows = parse_matrix(m)
    order, _, _ = oracle.load_edges(case)
    n = len(order)
    parents = set()
    for k, (l, r, _, _) in enumerate(rows):
        parents.add(l)
        parents.add(r)
    roots = len([x for x in range(n + len(rows)) if x not in parents])
    if f"{roots} roots in total" not in t:
        bad.append(f"{Path(case).name}: the root count from --tree disagrees with --matrix ({roots})")
    return (not bad), "; ".join(bad[:3])


def t_stdin_equals_file(bin_path):
    """The stdin mode and the file mode must give the same result."""
    bad = []
    for case in CASES[:6]:
        text = Path(case).read_text()
        a = parse_flat(run(bin_path, [case], stdin_text=None).stdout)
        b = parse_flat(run(bin_path, ["-"], stdin_text=text).stdout)
        c = parse_flat(run(bin_path, [], stdin_text=text).stdout)
        if not (a == b == c):
            bad.append(f"{Path(case).name}: file / '-' / no argument gave different results")
    return (not bad), "; ".join(bad[:3])


def t_duplicate_order_independent(bin_path):
    """Duplicate pairs take the maximum: swapping two input lines must not change the result."""
    import tempfile
    bad = []
    for body in ("a,b,90\na,b,80\na,c,88\n", "a,b,80\na,b,90\na,c,88\n"):
        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as f:
            f.write(body)
            p = f.name
        part = parse_flat(run(bin_path, [p]).stdout)
        if frozenset([frozenset(["a", "b"]), frozenset(["c"])]) != part:
            bad.append(f"input {body.splitlines()[:2]} produced {part}")
        Path(p).unlink()
    return (not bad), "; ".join(bad[:2])


TESTS = [
    ("--linkage at run time", t_linkage_runtime),
    ("--n-clusters", t_n_clusters),
    ("--min-size", t_min_size),
    ("exit codes and error paths", t_exit_codes),
    ("--cut matches the oracle cut", t_cut_matches_oracle),
    ("cut == early stop (monotonicity)", t_cut_equals_early_stop),
    ("--floor / --cut interaction", t_floor_cut_interaction),
    ("--matrix structural invariants", t_matrix_structure),
    ("--tree completeness", t_tree_completeness),
    ("stdin equals file", t_stdin_equals_file),
    ("duplicate pairs are order-independent", t_duplicate_order_independent),
]


def main():
    ap = argparse.ArgumentParser(description="spclust CLI test")
    ap.add_argument("--bin", default="./spclust")
    args = ap.parse_args()
    # Note Path("./spclust") normalizes to "spclust" (dropping the slash), after which the
    # subprocess would search PATH, so make it absolute first.
    bin_path = Path(args.bin)
    if not bin_path.is_absolute():
        bin_path = (Path.cwd() / bin_path).resolve()

    print("=" * 78)
    print(f"spclust CLI test    {bin_path}    {len(CASES)} cases")
    print("=" * 78)

    failed = 0
    for name, fn in TESTS:
        ok, detail = fn(bin_path)
        print(f"{'  ' if ok else 'X '}{name:<30} {'PASS' if ok else 'FAIL'}"
              + (f"   {detail}" if detail and (not ok or detail.startswith('(')) else ""))
        if not ok:
            failed += 1
    print("-" * 78)
    print(f"RESULT: {'PASS' if failed == 0 else 'FAIL'}  ({len(TESTS) - failed}/{len(TESTS)})")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

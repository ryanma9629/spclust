#!/usr/bin/env python3
"""
spclust differential test.

The core is hand-written malloc/free with four-way pointers and a cross-linked list,
which is exactly the kind of code that breaks at the edges, and whose correctness

cannot be established by reading it. So this does a three-way cross-check:

  1. oracle(brute) - by definition, recomputed from the dense matrix at every step. Most trustworthy, O(n^4).
  2. oracle(fast)  - incremental aggregates via the identity. brute/fast must agree first,
                     and only then is it trusted on large inputs.
  3. spclust       - the subject under test.

All three are independent implementations. spclust and the oracle must agree merge by merge
(member sets plus similarity), except for two explicable, legitimate divergences:
  * a tie for the maximum (the merge order is unspecified)
  * a tie created by truncating similarities to integers

Known, accepted divergences live in tests/known_divergences.txt and must state a reason.

There is also a regression guard: if spclust ever went back to storing similarities as unsigned
short (truncating every step, with accumulating error) the deviation is classified TRUNC and

fails - that bug is fixed and must not come back.
    make test
    python3 tests/run_tests.py --bin-prefix ./spclust-linkage- --min-sim 50
    python3 tests/run_tests.py --only average --verbose
"""

import argparse
import random
import subprocess
import sys
import tempfile
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))

import oracle  # noqa: E402

LINKAGE_BY_MTD = {0: "average", 1: "max", 2: "min"}
MTD_BY_LINKAGE = {v: k for k, v in LINKAGE_BY_MTD.items()}

# Above this size the O(n^4) brute oracle is not run
BRUTE_MAX_N = 24

# Random cases: (objects, edge density, seed)
RANDOM_CONFIGS = [
    (10, 0.30, 1), (10, 0.30, 2), (10, 0.30, 3),
    (18, 0.45, 1), (18, 0.45, 2), (18, 0.45, 3),
    (24, 0.65, 1), (24, 0.65, 2),
    (12, 1.00, 1), (12, 1.00, 2),
    (44, 0.25, 1), (44, 0.80, 1),
]

OUTCOME_ORDER = ("PASS", "WARN", "KNOWN", "TRUNC", "DRIFT", "XPASS", "FAIL", "ORACLE")

# Relative/absolute tolerance: the oracle uses exact fractions, spclust uses double, and
# linkage values accumulate in floating point. 1e-9 absorbs rounding while being far
# smaller than any real deviation (1 unit at magnitude ~100 is 1e-2).
#Thresholds for the multi-threshold cut check (integers and fractions)
CUT_LEVELS = (1, 10, 25, 50, 66.6, 75, 85, 90, 95, 99)

REL_TOL = 1e-9
ABS_TOL = 1e-9


def obj_names(member_ids, order):
    return [order[i] for i in sorted(member_ids)]


def gen_random_case(rng, n, density, tag):
    """Build one random sparse case. Similarities are distinct, to avoid ties as far as possible."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if rng.random() < density]
    if not pairs:
        return None
    k = len(pairs)
    if k <= 299:
        # Starting at 0 lets edges with sim=0 (the "implicit zero" path) be covered by
        # random cases too. The range matters: the use-after-free in
        # tests/cases/22_zero_row_max.csv could never be generated from the old 1..299 range
        vals = rng.sample(range(0, 300), k)
    else:
        vals = rng.sample(range(1, 65535), k)
    lines = [
        f"# random case n={n} density={density} tag={tag} (all similarities distinct)",
    ]
    lines += [f"r{i:02d},r{j:02d},{v}" for (i, j), v in zip(pairs, vals)]
    return "\n".join(lines) + "\n"


def load_known_divergences():
    path = HERE / "known_divergences.txt"
    known = {}
    if not path.exists():
        return known
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 2)
        if len(parts) < 3:
            raise SystemExit(f"known_divergences.txt is malformed (expected 'case linkage reason'): {raw!r}")
        known[(parts[0], parts[1])] = parts[2]
    return known


def collect_cases(args):
    cases = []
    for p in sorted((HERE / "cases").glob("*.csv")):
        cases.append((p.name, p.read_text(), None))
    if args.no_random:
        return cases
    for n, dens, seed in RANDOM_CONFIGS:
        if args.quick and n > BRUTE_MAX_N:
            continue
        rng = random.Random(n * 100000 + int(dens * 100) * 1000 + seed)
        tag = f"n{n}d{int(dens * 100)}s{seed}"
        text = gen_random_case(rng, n, dens, tag)
        if text:
            cases.append((f"rand_{tag}.csv", text, tag))
    return cases


def run_spclust(bin_path, case_path, cut, linkage=None):
    proc = subprocess.run(
        # --floor 0: compute the full tree so it can be compared merge by merge with the
        # oracle's full tree (stronger than comparing only above the threshold).
        # --cut: the cut position being verified.
        # --matrix: take the merges from the public linkage-matrix interface (no DEBUG-log parsing).
        [str(bin_path), str(case_path), "--floor", "0", "--cut", str(cut), "--matrix"]
        + ([] if linkage is None else ["--linkage", str(linkage)]),
        cwd=str(REPO), capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"spclust exited with {proc.returncode}:\n{proc.stderr}")
    return proc.stdout


def values_equal(oracle_val, spclust_val):
    """The oracle uses exact fractions, spclust uses double - compare by relative error."""
    a, b = float(oracle_val), float(spclust_val)
    return abs(a - b) <= ABS_TOL + REL_TOL * max(abs(a), abs(b))


def merge_multiset(res):
    """View the merge sequence as a multiset of merged set pairs (order ignored)."""
    out = []
    for (A, B), _v in res["merges"]:
        a, b = tuple(sorted(A)), tuple(sorted(B))
        out.append((a, b) if a < b else (b, a))
    return sorted(out)


def trunc_explains(o, o_trunc, got, o_trunc_cut):
    """Whether spclust's deviation from the exact definition is fully explained by unsigned-short

    truncation. Two effects are involved:
      a) each value is truncated (losing up to 1 per step, accumulating)
      b) truncation can collapse two distinct values (262.5 and 262.0) into one, creating a tie
         whose resolution is an implementation detail

    So it requires:
      1) up to the first truncation-induced tie, the truncated mirror matches spclust exactly
      2) both sides performed the same merges as a multiset
      3) the final partitions agree
    """
    ts = o["trunc_steps"]
    om_list, gm_list = o_trunc["merges"], got["merges"]
    if ts:
        # Only exact agreement up to the first truncation-induced tie is required; after that,
        limit = min(min(ts), len(om_list), len(gm_list))
    else:
        # a different order from that tie is allowed. With no such tie, all steps must agree.
        if len(om_list) != len(gm_list):
            return False
        limit = len(om_list)
    for k in range(limit):
        if om_list[k][0] != gm_list[k][0] \
                or not values_equal(om_list[k][1], gm_list[k][1]):
            return False
    return (merge_multiset(o_trunc) == merge_multiset(got)
            and o_trunc_cut == got["partition"])


def compare(o, got, order, o_cut):
    """Return (outcome, detail).

    o is the oracle's **full merge tree** (min_sim=0).
    o_cut is o cut at the current threshold; spclust prints its partition at the same threshold,
    so the two must agree. The merge sequences themselves are compared in full - stronger
    """
    om, gm = o["merges"], got["merges"]

    # Find the first divergence
    k = None
    for i in range(min(len(om), len(gm))):
        if om[i][0] != gm[i][0]:
            k = i
            break
    if k is None and len(om) != len(gm):
        k = min(len(om), len(gm))

    if k is None:
        drift = [i for i in range(len(om)) if not values_equal(om[i][1], gm[i][1])]
        if o_cut != got["partition"]:
            return "FAIL", "all merges agree but the partitions at this threshold differ"
        if drift:
            return "DRIFT", f"{len(drift)} steps differ in similarity, first at step {drift[0]}: " \
                           f"oracle={float(om[drift[0]][1]):.6f} spclust={float(gm[drift[0]][1]):.6f}"
        return "PASS", ""

    ties = sorted(set(o["tie_steps"]) | set(o["trunc_steps"]))
    tie_explains = any(t <= k for t in ties)

    detail = []
    if k < len(om) and k < len(gm):
        A, B = om[k][0]
        C, D = gm[k][0]
        detail.append(f"step {k}: oracle merged {obj_names(A, order)}|{obj_names(B, order)}"
                      f"(sim={float(om[k][1]):.6f}) "
                      f"vs spclust {obj_names(C, order)}|{obj_names(D, order)}(sim={float(gm[k][1]):.6f})")
    else:
        detail.append(f"stop position differs from step {k}: oracle {len(om)} merges, spclust {len(gm)}")
    detail.append("oracle partition=" + str(sorted(sorted(obj_names(g, order)) for g in o_cut)))
    detail.append("spclust partition=" + str(sorted(sorted(obj_names(g, order)) for g in got["partition"])))
    if ties:
        detail.append(f"tie/truncation at step {ties}")

    return ("WARN" if tie_explains else "FAIL"), "; ".join(detail)


def main(argv=None):
    ap = argparse.ArgumentParser(description="spclust differential test")
    ap.add_argument("--bin", default="./spclust",
                    help="binary under test, default ./spclust (pick the linkage with --linkage)")
    ap.add_argument("--min-sim", type=int, default=50,
                    help="must match the MIN_SIM_CLUST the binary was compiled with")
    ap.add_argument("--only", choices=sorted(LINKAGE_BY_MTD.values()),
                    help="test only one linkage")
    ap.add_argument("--case", help="only run cases whose filename contains this substring")
    ap.add_argument("--quick", action="store_true", help="skip the large cases")
    ap.add_argument("--no-random", action="store_true", help="only run the hand-written cases in tests/cases")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)

    linkages = [args.only] if args.only else sorted(LINKAGE_BY_MTD.values())

    # Check the binary
    binp = REPO / args.bin
    if not binp.exists():
        raise SystemExit(f"cannot find {binp}; run `make` first (or just use `make test`)")

    known = load_known_divergences()
    cases = collect_cases(args)
    if args.case:
        cases = [c for c in cases if args.case in c[0]]
    if not cases:
        raise SystemExit("no matching cases")

    print("=" * 78)
    print(f"spclust differential test   {len(cases)} cases   linkage {linkages}   min_sim={args.min_sim}")
    print("=" * 78)

    counts = {lk: {} for lk in linkages}
    failures = []
    notes = []
    selfcheck = 0
    cut_checks = 0
    tmpdir = tempfile.TemporaryDirectory()
    tmp = Path(tmpdir.name)

    for name, text, _tag in cases:
        case_path = tmp / name
        case_path.write_text(text)

        order, edges, warnings = oracle.load_edges(str(case_path))
        for w in warnings:
            notes.append(f"{name}: {w}")
        if len(order) < 2:
            notes.append(f"{name}: fewer than 2 objects, skipped")
            continue

        for lk in linkages:
            # Oracle self-check 1: on small inputs run both brute (recompute by definition each
            # step) and fast (incremental). They must agree, otherwise the harness itself is at
            # fault rather than spclust. The **full tree** (min_sim=0) is compared.
            use_brute = len(order) <= BRUTE_MAX_N
            mode = "brute" if use_brute else "fast"
            o = oracle.hac(order, edges, lk, 0, mode)
            o_cut = oracle.cut(o, args.min_sim)
            outcome, detail = None, ""
            if use_brute:
                o_fast = oracle.hac(order, edges, lk, 0, "fast")
                if o["merges"] != o_fast["merges"] or o["partition"] != o_fast["partition"]:
                    outcome, detail = "ORACLE", "the brute and fast oracle implementations disagree"
                else:
                    selfcheck += 1

            # Oracle self-check 2: monotonicity of cutting.
            # "cut the full tree at X" must equal "run with X as the stopping condition" -
            # which relies on similarity being non-increasing along the merge sequence (all three
            # linkages are reducible). If it failed, the cut semantics would be broken - and
            o_stop = oracle.hac(order, edges, lk, args.min_sim, mode)
            if o_stop["partition"] != o_cut:
                outcome, detail = "ORACLE", "cutting and early stopping gave different partitions (monotonicity?)"
            else:
                selfcheck += 1

            # Truncated mirror: check whether spclust's deviation from the exact definition
            o_trunc = oracle.hac(order, edges, lk, 0, "trunc")
            o_trunc_cut = oracle.cut(o_trunc, args.min_sim)

            if outcome is None:
                got = oracle.parse_spclust(run_spclust(binp, case_path, args.min_sim,
                                                       MTD_BY_LINKAGE[lk]), order)
                outcome, detail = compare(o, got, order, o_cut)

                if outcome != "PASS" and o_trunc["merges"] != o["merges"] \
                        and trunc_explains(o, o_trunc, got, o_trunc_cut):
                    # The deviation is explained by integer truncation at each step, i.e. the old
                    outcome = "TRUNC"
                    detail = "[REGRESSION] the deviation is explained by unsigned-short truncation; " \
                             "similarities should be double: " + detail
                elif (name, lk) in known and outcome in ("FAIL", "DRIFT", "WARN"):
                    detail = f"[known] {known[(name, lk)]} -- {detail}"
                    outcome = "KNOWN"
                elif (name, lk) in known and outcome == "PASS":
                    outcome = "XPASS"
                    detail = f"[known divergence is gone] {known[(name, lk)]}"

                # Multi-threshold cut consistency: since both sides have the **full tree**,
                # cut at several thresholds and compare - a direct test of --cut.
                # Only when the merge sequences agree exactly (PASS): with a tied maximum the order
                # may legitimately differ, and so may the resulting partitions.
                if outcome == "PASS":
                    got_tree = {"objects": o["objects"], "merges": got["merges"]}
                    for t in CUT_LEVELS:
                        if oracle.cut(got_tree, t) != oracle.cut(o, t):
                            outcome, detail = "FAIL", f"partitions differ at cut={t}"
                            break
                    else:
                        cut_checks += len(CUT_LEVELS)

            counts[lk][outcome] = counts[lk].get(outcome, 0) + 1
            if outcome in ("FAIL", "ORACLE", "TRUNC"):
                failures.append((lk, name, detail))

            flag = {"PASS": "  ", "WARN": "~ ", "KNOWN": "k ", "TRUNC": "t ", "DRIFT": "! ",
                    "XPASS": "? ", "FAIL": "X ", "ORACLE": "X "}[outcome]
            print(f"{flag}[{lk:^7}] {name:<26} {outcome:<7} merges={o['merge_count']:<3} "
                  f"groups={len(o['groups'])}"
                  + (f"  {detail}" if args.verbose or outcome in ("FAIL", "ORACLE", "DRIFT", "TRUNC") else ""))
            if o["duplicate_pairs"] and outcome in ("PASS", "TRUNC"):
                notes.append(f"{name}/{lk}: input has {o['duplicate_pairs']} duplicate pairs")

    tmpdir.cleanup()

    print("-" * 78)
    total = {k: 0 for k in OUTCOME_ORDER}
    for lk in linkages:
        row = "  ".join(f"{k}={counts[lk].get(k, 0)}"
                        for k in OUTCOME_ORDER
                        if counts[lk].get(k, 0))
        print(f"{lk:^7}: {row}")
        for k, v in counts[lk].items():
            total[k] += v
    print("-" * 78)
    print(f"oracle self-checks (brute/fast agree + cut == early stop): {selfcheck}")
    print(f"multi-threshold cut consistency checks: {cut_checks}")
    print(f"total: PASS={total['PASS']} WARN={total['WARN']} KNOWN={total['KNOWN']} "
          f"TRUNC={total['TRUNC']} DRIFT={total['DRIFT']} XPASS={total['XPASS']} "
          f"FAIL={total['FAIL']} ORACLE={total['ORACLE']}")

    if notes:
        print()
        print("notes:")
        for n in notes:
            print(f"  - {n}")

    if failures:
        print()
        print("failure detail:")
        for lk, name, detail in failures:
            print(f"  [{lk}] {name}")
            for part in detail.split("; "):
                print(f"      {part}")

    # XPASS (a known divergence that disappeared) is not a failure: it depends on the
    # threshold - e.g. a case that merges nothing at some threshold makes duplicate-pair
    # divergences vanish. It is still printed so stale records get cleaned up. TRUNC does
    ok = not failures and total["DRIFT"] == 0
    print()
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

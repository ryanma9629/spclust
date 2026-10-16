#!/usr/bin/env python3
"""
Performance regression guard.

Nothing prevented this repository's performance gains (about 39x time, -55% memory in total)
from rotting away, and machine noise alone is enough to fool a single measurement - it fooled
this project twice (too few repetitions) and once (attribution from sampling). So the
**protocol** and the **baseline** are pinned down here.
Protocol (pinned here, do not bypass):
  * inputs come from the fixed generator in tests/gen_bench.py (its parameters are pinned);
  * each case is repeated N times and the median taken; **if min / lower quartile / median
    differ by more than 15% it is reported as noisy** and no conclusion is drawn;
  * comparison against tests/bench_baseline.txt only **warns** over tolerance; only a FATAL
    factor fails (numbers drift across machines and load; --strict fails on tolerance).

Usage:
    make bench                          # run once and compare with the baseline
    make bench BENCH_ARGS=--quick       # skip the largest cases
    python3 tests/bench.py --update     # re-record the baseline (new machine, or a real change)
"""

import argparse
import os
import platform
import re
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
import gen_bench  # noqa: E402

BASELINE = HERE / "bench_baseline.txt"

# Tolerances: time is loose (load-sensitive), memory is tight because it is stable
TOL_TIME = 1.25
TOL_RSS = 1.15
# Beyond this factor it is a real regression (even without --strict); 3x is far above noise
FATAL_TIME = 3.0
FATAL_RSS = 2.0
CANARY_CASE = "bench_small"
CANARY_TOL = 0.30

# Machine-state canary: bench_small is the fastest fixed-input case, so it probes "how fast is
# this machine right now". It caught a real trap during development: a baseline recorded on a
# busy/hot machine, rerun idle, showed everything "2x faster" - which looked like a huge
# optimisation but was the environment. The reverse false-alarms everything as slower. So a
# canary deviation beyond this ratio means no conclusion is drawn at all (a canary only
# catches gross changes; see docs/testing.md).

# Repetitions per case: fewer for the large ones, to bound the total time. They are slow, and
# with 7 repetitions the gap between min and median reached 23% (something else was running),
# so a few extra seconds are worth reaching the true minimum.
REPS = {"bench_small": 15, "bench_mid": 13, "bench_large": 11,
        "bench_lowval": 13, "bench_single": 11}

# Which binary to run, and with which arguments
PLAN = [
    ("bench_small",  "spclust",   ["--floor", "50"]),
    ("bench_mid",    "spclust",   ["--floor", "50"]),
    ("bench_large",  "spclust",   ["--floor", "50"]),
    ("bench_lowval", "spclust",   ["--floor", "50"]),
    ("bench_single", "linkage1",  ["--floor", "50", "--linkage", "1"]),  # single linkage uses Kruskal
]


def machine_fingerprint():
    try:
        chip = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                              capture_output=True, text=True).stdout.strip()
    except Exception:
        chip = ""
    if not chip:
        chip = platform.processor() or "unknown"
    return f"{chip} / {platform.system()} {platform.release()}"


def peak_rss_mb(cmd):
    """Peak RSS of a child process, in MB. Uses time -l on macOS and time -v on Linux."""
    out = subprocess.run(["/usr/bin/time", "-l"] + cmd, capture_output=True, text=True).stderr
    m = re.search(r"(\d+)\s+maximum resident set size", out)
    if m:
        return int(m.group(1)) / 1048576.0
    out = subprocess.run(["/usr/bin/time", "-v"] + cmd, capture_output=True, text=True).stderr
    m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", out)
    if m:
        return int(m.group(1)) / 1024.0
    return None


def measure(binary, case_path, extra, reps):
    """Run repeatedly; returns (min, lower quartile, median, peak RSS in MB)."""
    times = []
    for _ in range(reps):
        t0 = time.perf_counter()
        proc = subprocess.run([str(binary), str(case_path)] + extra, capture_output=True)
        times.append(time.perf_counter() - t0)
        if proc.returncode != 0:
            raise SystemExit(f"{binary} exited with {proc.returncode}: {proc.stderr[:200]}")
    times.sort()
    return times[0], times[len(times) // 4], statistics.median(times), peak_rss_mb(
        [str(binary), str(case_path)] + extra)


def load_baseline(path):
    if not path.exists():
        return None, {}
    lines = path.read_text().splitlines()
    fp = ""
    data = {}
    for line in lines:
        line = line.strip()
        if line.startswith("# machine:"):
            fp = line.split(":", 1)[1].strip()
        elif line and not line.startswith("#"):
            parts = line.split()
            if len(parts) == 4:
                data[(parts[0], parts[1])] = (float(parts[2]), float(parts[3]))
    return fp, data


def main():
    ap = argparse.ArgumentParser(description="spclust performance regression guard")
    ap.add_argument("--bin", default="./spclust", help="binary to measure")
    ap.add_argument("--baseline", default=str(BASELINE))
    ap.add_argument("--update", action="store_true", help="re-record the baseline")
    ap.add_argument("--strict", action="store_true", help="fail on tolerance (default only warns)")
    ap.add_argument("--quick", action="store_true", help="skip the largest cases")
    args = ap.parse_args()

    baseline_path = Path(args.baseline)

    # Path("./x") normalizes to "x" and the subprocess would then search PATH, so make it absolute
    def absbin(p):
        p = Path(p)
        return p if p.is_absolute() else (Path.cwd() / p).resolve()
    # The linkage is a run-time choice, so all three share one binary (column 2 of PLAN is just a
    bins = {"spclust": absbin(args.bin), "linkage1": absbin(args.bin)}
    for name, b in bins.items():
        if not b.exists():
            raise SystemExit(f"cannot find {b}; run `make` first (or just use `make bench`)")

    old_fp, old = load_baseline(baseline_path)
    fp = machine_fingerprint()
    if old_fp and old_fp != fp and not args.update:
        print(f"note: the baseline machine is [{old_fp}], this one is [{fp}] - not directly comparable.\n")

    print("=" * 78)
    print(f"spclust performance guard    machine: {fp}")
    print("=" * 78)

    rows = []
    noisy = []
    exits = []
    tmpdir = tempfile.TemporaryDirectory()
    tmp = Path(tmpdir.name)

    for case, which, extra in PLAN:
        if args.quick and case in ("bench_large", "bench_single"):
            continue
        path = tmp / (case + ".csv")
        path.write_text(gen_bench.generate(case))
        reps = REPS.get(case, 9)
        mn, q1, md, rss = measure(bins[which], path, extra, reps)

        spread = (md - mn) / md if md > 0 else 0.0
        flag = ""
        if spread > 0.15:
            flag = f"  noisy (min/median differ by {spread*100:.0f}%)"
            noisy.append(case)

        key = (case, which)
        if args.update:
            row = ("record", mn, q1, md, rss, None, None)
        elif key in old:
            o_t, o_r = old[key]
            rt = md / o_t
            rr = (rss / o_r) if (rss and o_r) else None
            verdict = "ok"
            if rt > FATAL_TIME or (rr and rr > FATAL_RSS):
                verdict = "REGRESSION"
                exits.append(case)
            elif rt > TOL_TIME or (rr and rr > TOL_RSS):
                verdict = "slower"
                if args.strict:
                    exits.append(case)
            row = (verdict, mn, q1, md, rss, rt, rr)
        else:
            row = ("new", mn, q1, md, rss, None, None)
        rows.append((case, which, reps, row, flag))

    print(f"{'case':<14}{'binary':<10}{'n':>7}{'median_s':>9}{'min_s':>9}{'q1_s':>10}{'RSS MB':>9}"
          f"{'vs base':>9}{'RSS x':>8}")
    for case, which, reps, row, flag in rows:
        n = gen_bench.case_params(case)[1]
        verdict, mn, q1, md, rss, rt, rr = row
        rts = f"{rt:.2f}x" if rt else "-"
        rrs = f"{rr:.2f}x" if rr else "-"
        print(f"{case:<14}{which:<10}{n:>7}{md:>9.4f}{mn:>9.4f}{q1:>10.4f}{rss:>9.1f}"
              f"{rts:>9}{rrs:>8}  {verdict}{flag}")

    if args.update:
        with open(baseline_path, "w") as f:
            f.write("# spclust performance baseline, produced by `python3 tests/bench.py --update`\n")
            f.write("# comparable only when the tests/gen_bench.py parameters and the machine are the same\n")
            f.write(f"# machine: {fp}\n")
            f.write("# protocol: N repetitions per case, median taken; min/q1/median differing >15% = noisy\n")
            f.write("# columns: case binary median_wall_s peak_RSS_MB\n")
            for case, which, reps, row, flag in rows:
                _, _, _, md, rss, _, _ = row
                f.write(f"{case:<14} {which:<10} {md:.4f} {rss:.1f}\n")
        print(f"\nbaseline written to: {baseline_path}")
        if noisy:
            print(f"note: {len(noisy)} cases were noisy this run ({', '.join(noisy)}) -"
                  f" those baseline numbers are imprecise; re-record them on an idle machine")
        return 0

    print("-" * 78)

    # Canary check: if the machine state differs too much from the baseline, no ratio below counts
    canary_ratio = None
    for case, which, reps, row, flag in rows:
        if case == CANARY_CASE and row[5]:
            canary_ratio = row[5]
    if canary_ratio is not None and abs(canary_ratio - 1.0) > CANARY_TOL:
        # rt = this median / baseline median; below 1 means less time, i.e. **faster**
        direction = "faster" if canary_ratio < 1 else "slower"
        print(f"no conclusion: the machine state differs from the baseline - canary {CANARY_CASE} is {direction} "
              f"{abs(canary_ratio - 1) * 100:.0f}% ({canary_ratio:.2f}x)."
              f"this time. Rerun on an idle machine, or re-record with `make bench-update`.")
        return 0

    if noisy:
        print(f"warning: {len(noisy)} noisy cases ({', '.join(noisy)}) - not treated as conclusions,"
              f" rerun them on an idle machine")
    if exits:
        print(f"FAILED: {len(exits)} cases out of tolerance ({', '.join(exits)})")
        return 1
    print("conclusion: consistent with the baseline (over-tolerance cases are marked slower; only "
          f"{FATAL_TIME:.0f}x time / {FATAL_RSS:.0f}x memory fails)")
    if canary_ratio:
        print(f"(canary {CANARY_CASE} {canary_ratio:.2f}x, tolerance plus/minus {CANARY_TOL * 100:.0f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

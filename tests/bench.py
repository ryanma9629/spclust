#!/usr/bin/env python3
"""
性能回归护栏。

这个仓库的性能收益（累计约 39x 时间 / -55% 内存）没有任何东西防止它退化，
而机器噪声又足以骗过单次测量 —— 本项目开发过程中就被骗过两次（重复次数不够）
和一次（采样归因）。所以把**协议**和**基线**固化下来。

协议（写死在这里，不要绕过）：
  * 输入来自 tests/gen_bench.py 的固定生成器（参数写死，见那个文件）
  * 每个用例重复 N 次，取中位数；**最小 / 下四分位 / 中位数三者差异超过 15% 就报"噪声大"**，
    这时不给出结论（说明机器上还有别的东西在跑）
  * 与 tests/bench_baseline.txt 比对，超出容差只**警告**；超过 FATAL 倍才失败
    （默认不硬失败，因为跨机器/负载的数字本来就会漂；`--strict` 可按容差硬失败）

用法:
    make bench                          # 跑一遍并与基线比较
    make bench BENCH_ARGS=--quick       # 跳过最大的用例
    python3 tests/bench.py --update     # 重新记录基线（改变机器或有意改性能后用）
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

# 容差：时间放宽一些（受负载影响大），内存很稳所以收紧
TOL_TIME = 1.25
TOL_RSS = 1.15
# 超过这个倍数视为真正的回归（哪怕没开 --strict）；3 倍远大于任何噪声
FATAL_TIME = 3.0
FATAL_RSS = 2.0

# 机器状态哨兵：bench_small 是固定输入的最快用例，把它当成"这台机器现在多快"的探针。
# 开发过程中踩过一次：基线是在机器忙/热的时候记的，之后空闲时跑，所有用例都显示"快 2 倍"，
# 看上去像大幅优化，其实是环境变了。反之（空闲时记、忙时跑）会全线误报"变慢"。
# 所以哨兵偏差超过这个比例时，本次与基线的比较一律不作结论。
CANARY_CASE = "bench_small"
CANARY_TOL = 0.30

# 每个用例跑几次（大的少跑几次，控制总时长）
# 大的用例本身慢，重复次数只能少一些；但实测 7 次时 min 与中位的差能到 23%（机器上
# 还有别的东西在跑），所以宁可多花几秒也要把最小值探到底。
REPS = {"bench_small": 15, "bench_mid": 13, "bench_large": 11,
        "bench_lowval": 13, "bench_single": 11}

# 跑哪个二进制、用什么参数
PLAN = [
    ("bench_small",  "spclust",   ["--floor", "50"]),
    ("bench_mid",    "spclust",   ["--floor", "50"]),
    ("bench_large",  "spclust",   ["--floor", "50"]),
    ("bench_lowval", "spclust",   ["--floor", "50"]),
    ("bench_single", "linkage1",  ["--floor", "50"]),   # 单链接走 Kruskal 路径
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
    """子进程峰值 RSS（MB）。macOS 用 time -l，Linux 用 time -v。"""
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
    """重复执行，返回 (最小, 下四分位, 中位, 峰值RSS_MB)。"""
    times = []
    for _ in range(reps):
        t0 = time.perf_counter()
        proc = subprocess.run([str(binary), str(case_path)] + extra, capture_output=True)
        times.append(time.perf_counter() - t0)
        if proc.returncode != 0:
            raise SystemExit(f"{binary} 退出码 {proc.returncode}: {proc.stderr[:200]}")
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
        if line.startswith("# 机器:"):
            fp = line.split(":", 1)[1].strip()
        elif line and not line.startswith("#"):
            parts = line.split()
            if len(parts) == 4:
                data[(parts[0], parts[1])] = (float(parts[2]), float(parts[3]))
    return fp, data


def main():
    ap = argparse.ArgumentParser(description="spclust 性能回归护栏")
    ap.add_argument("--bin", default="./spclust", help="默认二进制")
    ap.add_argument("--bin-linkage1", default="./spclust-linkage-50-1",
                    help="单链接（Kruskal）用的二进制")
    ap.add_argument("--baseline", default=str(BASELINE))
    ap.add_argument("--update", action="store_true", help="重新记录基线")
    ap.add_argument("--strict", action="store_true", help="超出容差即失败（默认只警告）")
    ap.add_argument("--quick", action="store_true", help="跳过最大的用例")
    args = ap.parse_args()

    baseline_path = Path(args.baseline)

    # Path("./x") 会规范化成 "x"，子进程就会去 PATH 找 —— 必须先变绝对路径
    def absbin(p):
        p = Path(p)
        return p if p.is_absolute() else (Path.cwd() / p).resolve()
    bins = {"spclust": absbin(args.bin), "linkage1": absbin(args.bin_linkage1)}
    for name, b in bins.items():
        if not b.exists():
            raise SystemExit(f"找不到 {b}；先 `make linkages`（或直接用 `make bench`）")

    old_fp, old = load_baseline(baseline_path)
    fp = machine_fingerprint()
    if old_fp and old_fp != fp and not args.update:
        print(f"注意: 基线的机器是 [{old_fp}]，当前是 [{fp}] —— 数字不可直接比。\n")

    print("=" * 78)
    print(f"spclust 性能护栏    机器: {fp}")
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
            flag = f"  噪声大(min/中位差 {spread*100:.0f}%)"
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
                verdict = "回归!"
                exits.append(case)
            elif rt > TOL_TIME or (rr and rr > TOL_RSS):
                verdict = "变慢"
                if args.strict:
                    exits.append(case)
            row = (verdict, mn, q1, md, rss, rt, rr)
        else:
            row = ("新增", mn, q1, md, rss, None, None)
        rows.append((case, which, reps, row, flag))

    print(f"{'用例':<14}{'二进制':<10}{'n':>7}{'中位s':>9}{'最小s':>9}{'下四分位':>10}{'RSS MB':>9}"
          f"{'vs基线':>9}{'RSS比':>8}")
    for case, which, reps, row, flag in rows:
        n = gen_bench.case_params(case)[1]
        verdict, mn, q1, md, rss, rt, rr = row
        rts = f"{rt:.2f}x" if rt else "-"
        rrs = f"{rr:.2f}x" if rr else "-"
        print(f"{case:<14}{which:<10}{n:>7}{md:>9.4f}{mn:>9.4f}{q1:>10.4f}{rss:>9.1f}"
              f"{rts:>9}{rrs:>8}  {verdict}{flag}")

    if args.update:
        with open(baseline_path, "w") as f:
            f.write("# spclust 性能基线：由 `python3 tests/bench.py --update` 生成\n")
            f.write("# 只在 tests/gen_bench.py 的参数与机器都相同时才可比\n")
            f.write(f"# 机器: {fp}\n")
            f.write("# 协议: 每用例重复 N 次，取中位；min/下四分位/中位差异 >15% 视为噪声大\n")
            f.write("# 列: 用例 二进制 中位墙钟秒 峰值RSS_MB\n")
            for case, which, reps, row, flag in rows:
                _, _, _, md, rss, _, _ = row
                f.write(f"{case:<14} {which:<10} {md:.4f} {rss:.1f}\n")
        print(f"\n已写入基线: {baseline_path}")
        if noisy:
            print(f"提醒: {len(noisy)} 个用例这次噪声大（{', '.join(noisy)}）——"
                  f"基线里这几个数字的误差就大，建议机器空闲时重记")
        return 0

    print("-" * 78)

    # 哨兵检查：机器状态和记基线时差太多的话，下面的比值都不作数
    canary_ratio = None
    for case, which, reps, row, flag in rows:
        if case == CANARY_CASE and row[5]:
            canary_ratio = row[5]
    if canary_ratio is not None and abs(canary_ratio - 1.0) > CANARY_TOL:
        # rt = 本次中位 / 基线中位；比值小于 1 说明耗时更少，也就是**更快**
        direction = "快" if canary_ratio < 1 else "慢"
        print(f"不作结论: 机器状态与记基线时明显不同 —— 哨兵用例 {CANARY_CASE} 本次{direction} "
              f"{abs(canary_ratio - 1) * 100:.0f}%（{canary_ratio:.2f}x）。"
              f"要么等机器空闲后重跑，要么用 `make bench-update` 重记基线。")
        return 0

    if noisy:
        print(f"警告: {len(noisy)} 个用例噪声大（{', '.join(noisy)}）—— 这些数字不作为结论，"
              f"建议机器空闲时重跑")
    if exits:
        print(f"失败: {len(exits)} 个用例超出容差（{', '.join(exits)}）")
        return 1
    print("结论: 与基线一致（超出容差的会标成『变慢』；只有超过 "
          f"{FATAL_TIME:.0f}x 时间 / {FATAL_RSS:.0f}x 内存才失败）")
    if canary_ratio:
        print(f"（哨兵 {CANARY_CASE} {canary_ratio:.2f}x，容差 ±{CANARY_TOL * 100:.0f}%）")
    return 0


if __name__ == "__main__":
    sys.exit(main())

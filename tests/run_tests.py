#!/usr/bin/env python3
"""
spclust 差分测试。

思路：spclust.cpp 是一个手写 malloc/free + 四向指针 + 十字链表的实现，
这类代码最容易在边界上出错，而它的正确性无法靠读代码确认。

所以这里做三方交叉验证：

  1. oracle(brute) —— 完全按定义、每步从原始稠密矩阵重算。最可信，但 O(n^4)。
  2. oracle(fast)  —— 用恒等式增量维护聚合量。先要求 brute/fast 一致，
                      一致后才有资格去测大数据。
  3. spclust       —— 被测对象。

三者独立实现。spclust 与 oracle 的逐次合并（成员集合 + 相似度）必须一致，
除了两种可解释的合法分歧：
  * 并列最大（合并顺序不确定）
  * 相似度被截断成整数后产生并列

已知/已确认的分歧放在 tests/known_divergences.txt 里，需要写明理由。

另有一重回归护栏：若 spclust 退回用 unsigned short 存相似度（每步截断、误差会累积），
差分测试会把偏差归类为 TRUNC 并**判定失败** —— 那个 bug 已修复，不该回来。

用法：
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

# 超过这个规模就不跑 O(n^4) 的 brute oracle 了
BRUTE_MAX_N = 24

# 随机用例：(对象数, 边密度, 种子)
RANDOM_CONFIGS = [
    (10, 0.30, 1), (10, 0.30, 2), (10, 0.30, 3),
    (18, 0.45, 1), (18, 0.45, 2), (18, 0.45, 3),
    (24, 0.65, 1), (24, 0.65, 2),
    (12, 1.00, 1), (12, 1.00, 2),
    (44, 0.25, 1), (44, 0.80, 1),
]

OUTCOME_ORDER = ("PASS", "WARN", "KNOWN", "TRUNC", "DRIFT", "XPASS", "FAIL", "ORACLE")

# 相对/绝对容差：oracle 用精确分数，spclust 用 double，链接值是浮点累加出来的。
# 1e-9 足够宽到容忍舍入误差，又远小于任何真实偏差
#（1 个单位在 ~100 的量级上就是 1e-2）。
#多阈值切树检查用的阈值（覆盖整数与小数）
CUT_LEVELS = (1, 10, 25, 50, 66.6, 75, 85, 90, 95, 99)

REL_TOL = 1e-9
ABS_TOL = 1e-9


def obj_names(member_ids, order):
    return [order[i] for i in sorted(member_ids)]


def gen_random_case(rng, n, density, tag):
    """生成一个随机稀疏图用例。相似度取不重复的值，尽量避开并列。"""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if rng.random() < density]
    if not pairs:
        return None
    k = len(pairs)
    if k <= 299:
        # 从 0 开始取：让 sim=0 的边（="隐式零" 路径）也有机会被随机用例覆盖。
        # 这个区间的改动是有理由的 —— tests/cases/22_zero_row_max.csv 那个 UAF 就是
        # 老的 1..299 区间永远生成不出来的配置
        vals = rng.sample(range(0, 300), k)
    else:
        vals = rng.sample(range(1, 65535), k)
    lines = [
        f"# 随机用例 n={n} density={density} tag={tag}（相似度互不相同）",
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
            raise SystemExit(f"known_divergences.txt 格式错误（应为 '用例 linkage 理由'）: {raw!r}")
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


def run_spclust(bin_path, case_path, cut):
    proc = subprocess.run(
        # --floor 0：算完整树，这样能与 oracle 的完整树逐次比对（比只比阈值以上更强）。
        # --cut：本次要验证的切树位置。
        # --matrix：走公开的 linkage 矩阵接口拿逐次合并（不再解析 DEBUG 日志）
        [str(bin_path), str(case_path), "--floor", "0", "--cut", str(cut), "--matrix"],
        cwd=str(REPO), capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"spclust 退出码 {proc.returncode}:\n{proc.stderr}")
    return proc.stdout


def values_equal(oracle_val, spclust_val):
    """oracle 是精确分数，spclust 是 double —— 按相对误差比较。"""
    a, b = float(oracle_val), float(spclust_val)
    return abs(a - b) <= ABS_TOL + REL_TOL * max(abs(a), abs(b))


def merge_multiset(res):
    """把合并序列看成“合并了哪些集合对”的多重集（忽略顺序）。"""
    out = []
    for (A, B), _v in res["merges"]:
        a, b = tuple(sorted(A)), tuple(sorted(B))
        out.append((a, b) if a < b else (b, a))
    return sorted(out)


def trunc_explains(o, o_trunc, got, o_trunc_cut):
    """spclust 与精确定义的偏差是否完全由 unsigned short 截断解释。

    包含两件事：
      a) 逐次值被截断（每步少最多 1，且会累积）
      b) 截断可能把两个不同的值（如 262.5 和 262.0）压成同一个值，产生并列，
         而并列的打破顺序是实现细节

    因此要求：
      1) 在第一次“截断造成并列”之前，截断镜像与 spclust 逐次完全一致（集合 + 值）
      2) 两边执行的合并（作为多重集）相同
      3) 最终分组相同
    """
    ts = o["trunc_steps"]
    om_list, gm_list = o_trunc["merges"], got["merges"]
    if ts:
        # 只要求“第一次截断造成并列”之前完全一致，之后允许并列导致的顺序不同
        limit = min(min(ts), len(om_list), len(gm_list))
    else:
        # 没有截断并列，那就得全程一致
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
    """返回 (outcome, detail)。

    o 是 oracle 的**完整合并树**（min_sim=0）。
    o_cut 是把 o 在本次阈值处切开得到的分区；spclust 打印的分区也切在同一阈值，
    所以两者应当一致。合并序列本身则做完整比对 —— 比只比"阈值以上那段"更强。
    """
    om, gm = o["merges"], got["merges"]

    # 找第一次分歧
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
            return "FAIL", "merge 全一致但在该阈值处切开的分组不同"
        if drift:
            return "DRIFT", f"{len(drift)} 步相似度不一致，第 {drift[0]} 步: " \
                           f"oracle={float(om[drift[0]][1]):.6f} spclust={float(gm[drift[0]][1]):.6f}"
        return "PASS", ""

    ties = sorted(set(o["tie_steps"]) | set(o["trunc_steps"]))
    tie_explains = any(t <= k for t in ties)

    detail = []
    if k < len(om) and k < len(gm):
        A, B = om[k][0]
        C, D = gm[k][0]
        detail.append(f"第 {k} 步: oracle 合并 {obj_names(A, order)}|{obj_names(B, order)}"
                      f"(sim={float(om[k][1]):.6f}) "
                      f"vs spclust {obj_names(C, order)}|{obj_names(D, order)}(sim={float(gm[k][1]):.6f})")
    else:
        detail.append(f"第 {k} 步起停止位置不同: oracle 共 {len(om)} 次合并，spclust {len(gm)} 次")
    detail.append("oracle 分组=" + str(sorted(sorted(obj_names(g, order)) for g in o_cut)))
    detail.append("spclust 分组=" + str(sorted(sorted(obj_names(g, order)) for g in got["partition"])))
    if ties:
        detail.append(f"并列/截断发生在第 {ties} 步")

    return ("WARN" if tie_explains else "FAIL"), "; ".join(detail)


def main(argv=None):
    ap = argparse.ArgumentParser(description="spclust 差分测试")
    ap.add_argument("--bin-prefix", default="./spclust-linkage-",
                    help="按 LINKAGE_MTD 编号的二进制前缀，默认 ./spclust-linkage-")
    ap.add_argument("--min-sim", type=int, default=50,
                    help="必须与二进制编译时的 MIN_SIM_CLUST 一致")
    ap.add_argument("--only", choices=sorted(LINKAGE_BY_MTD.values()),
                    help="只测某一种 linkage")
    ap.add_argument("--case", help="只跑文件名包含该子串的用例")
    ap.add_argument("--quick", action="store_true", help="跳过大规模用例")
    ap.add_argument("--no-random", action="store_true", help="只跑 tests/cases 下的手写用例")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)

    linkages = [args.only] if args.only else sorted(LINKAGE_BY_MTD.values())

    # 检查二进制
    for lk in linkages:
        p = REPO / f"{args.bin_prefix}{MTD_BY_LINKAGE[lk]}"
        if not p.exists():
            raise SystemExit(f"找不到 {p}，先跑 `make linkages`（或直接用 `make test`）")

    known = load_known_divergences()
    cases = collect_cases(args)
    if args.case:
        cases = [c for c in cases if args.case in c[0]]
    if not cases:
        raise SystemExit("没有匹配的用例")

    print("=" * 78)
    print(f"spclust 差分测试   用例 {len(cases)} 个   linkage {linkages}   min_sim={args.min_sim}")
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
            notes.append(f"{name}: 对象少于 2 个，跳过")
            continue

        for lk in linkages:
            # oracle 自检 1：小规模时同时跑 brute（每步按定义重算）和 fast（增量维护），
            # 两者必须逐次一致，否则是测试工具自身的问题，而不是 spclust 的问题。
            # 这里算的是**完整树**（min_sim=0），所以比对的是整棵树，比只比阈值以上那段更强。
            use_brute = len(order) <= BRUTE_MAX_N
            mode = "brute" if use_brute else "fast"
            o = oracle.hac(order, edges, lk, 0, mode)
            o_cut = oracle.cut(o, args.min_sim)
            outcome, detail = None, ""
            if use_brute:
                o_fast = oracle.hac(order, edges, lk, 0, "fast")
                if o["merges"] != o_fast["merges"] or o["partition"] != o_fast["partition"]:
                    outcome, detail = "ORACLE", "brute 与 fast 两个 oracle 实现不一致"
                else:
                    selfcheck += 1

            # oracle 自检 2：切树的单调性。
            # "把完整树切在 X 处" 必须等于 "直接用 X 当停止条件跑一遍" ——
            # 这依赖相似度沿合并序列单调不增（三种链接都是可归约的）。不成立就说明
            # 切树语义有问题，而 spclust 的 --cut 正是这么切的。
            o_stop = oracle.hac(order, edges, lk, args.min_sim, mode)
            if o_stop["partition"] != o_cut:
                outcome, detail = "ORACLE", "切树与提前停止得到的分组不一致（单调性不成立？）"
            else:
                selfcheck += 1

            # “截断镜像”：验证 spclust 与精确定义的偏差是否完全来自 unsigned short 截断
            o_trunc = oracle.hac(order, edges, lk, 0, "trunc")
            o_trunc_cut = oracle.cut(o_trunc, args.min_sim)

            if outcome is None:
                got = oracle.parse_spclust(run_spclust(REPO / f"{args.bin_prefix}{MTD_BY_LINKAGE[lk]}",
                                                       case_path, args.min_sim), order)
                outcome, detail = compare(o, got, order, o_cut)

                if outcome != "PASS" and o_trunc["merges"] != o["merges"] \
                        and trunc_explains(o, o_trunc, got, o_trunc_cut):
                    # 偏差能由“每步截断成整数”解释 —— 即旧的 unsigned short 问题回来了
                    outcome = "TRUNC"
                    detail = "【回归】偏差可由 unsigned short 整数截断解释；" \
                             "spclust 的相似度应当是 double: " + detail
                elif (name, lk) in known and outcome in ("FAIL", "DRIFT", "WARN"):
                    detail = f"[已知] {known[(name, lk)]} -- {detail}"
                    outcome = "KNOWN"
                elif (name, lk) in known and outcome == "PASS":
                    outcome = "XPASS"
                    detail = f"[已知分歧已消失] {known[(name, lk)]}"

                # 多阈值切树一致性：既然两边都拿到了**完整树**，就可以在若干阈值处
                # 各自切一刀再比 —— 这直接验证 --cut 这个功能。
                # 只在合并序列完全一致（PASS）时做：有并列最大值的用例合并顺序本就可能不同，
                # 切出来的分组自然也可能不同。
                if outcome == "PASS":
                    got_tree = {"objects": o["objects"], "merges": got["merges"]}
                    for t in CUT_LEVELS:
                        if oracle.cut(got_tree, t) != oracle.cut(o, t):
                            outcome, detail = "FAIL", f"在 cut={t} 处切开的分组不一致"
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
                notes.append(f"{name}/{lk}: 输入含 {o['duplicate_pairs']} 条重复 pair")

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
    print(f"oracle 自检（brute/fast 一致 + 切树=提前停止）: {selfcheck} 次")
    print(f"多阈值切树一致性检查: {cut_checks} 次")
    print(f"合计: PASS={total['PASS']} WARN={total['WARN']} KNOWN={total['KNOWN']} "
          f"TRUNC={total['TRUNC']} DRIFT={total['DRIFT']} XPASS={total['XPASS']} "
          f"FAIL={total['FAIL']} ORACLE={total['ORACLE']}")

    if notes:
        print()
        print("备注:")
        for n in notes:
            print(f"  - {n}")

    if failures:
        print()
        print("失败详情:")
        for lk, name, detail in failures:
            print(f"  [{lk}] {name}")
            for part in detail.split("; "):
                print(f"      {part}")

    # XPASS（已知分歧消失）不算失败：它依赖阈值 —— 例如某用例在某阈值下根本不合并，
    # 重复 pair 之类的分歧就自然消失了。但它会被打印出来，便于清理过期的记录。
    # TRUNC 则算失败：它意味着 unsigned short 截断那个 bug 又回来了。
    ok = not failures and total["DRIFT"] == 0
    print()
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

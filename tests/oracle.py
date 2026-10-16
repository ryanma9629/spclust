#!/usr/bin/env python3
"""
spclust 的独立参考实现（oracle），用于差分测试。

刻意按 *定义* 实现，不抄 spclust.cpp 里那套 Lance-Williams 递推 + 十字链表，
这样两边是独立算出来的，结果必须逐字一致。

模型：相似度矩阵 S，缺失的 pair 视为 0（隐式零）。两个簇 A、B 的相似度：

    average : sum_{a in A, b in B} S(a,b) / (|A| * |B|)     分母是全部 |A|*|B| 对
    max     : max_{a in A, b in B} S(a,b)                   缺失算 0
    min     : min_{a in A, b in B} S(a,b)                   仅当结果 > 0 时才是候选

对应 spclust 的 LINKAGE_MTD = 0 / 1 / 2。三个定义都与 sparseAHC 的文档一致
（average 用全对分母，min/max 对缺失边不贡献）。

每步取全局最大的簇对合并；若该值 < min_sim 则停止。

三种模式：
  brute —— 每步直接从原始稠密矩阵按定义重算。最可信，但 O(n^4)，只适合小数据。
  fast  —— 用恒等式增量维护：
             sum(A∪B, C) = sum(A,C) + sum(B,C)
             max(A∪B, C) = max(max(A,C), max(B,C))
             min(A∪B, C) = min(min(A,C), min(B,C))
           然后每步仍按定义（除以 |A|*|B|）算出链接值。适合大数据。
  trunc —— 镜像 spclust **修复前**的 unsigned short 整数语义。
           它不是独立 oracle，而是回归护栏：若 spclust 又退回整数截断，
           差分测试会用它把偏差归类为 TRUNC 并判定失败。
runner 会先用 brute 跑全部小用例，并交叉验证 fast 与 brute 一致。

全程用 Fraction 精确运算（也允许小数输入，如 87.5），不引入浮点误差。

命令行用法（调试用）：
    python3 tests/oracle.py sim.csv --linkage average --min-sim 50
"""

import json
import re
import sys
from fractions import Fraction

LINKAGES = ("average", "max", "min")

# spclust 的 LINKAGE_MTD -> 本文件的 linkage 名
LINKAGE_BY_MTD = {0: "average", 1: "max", 2: "min"}


def to_fraction(s):
    """把字符串解析成精确分数。不合法返回 None。"""
    try:
        return Fraction(s.strip())  # "90" / "87.5" / "175/2"
    except (ValueError, ZeroDivisionError):
        pass
    try:
        # 容忍 "1e-05" 之类的指数写法；nan / inf 会在这里被挡掉
        return Fraction(*float(s).as_integer_ratio())
    except (ValueError, OverflowError):
        return None


def load_edges(path):
    """读入边表，返回 (order, edges, warnings)。

    order 是按“首次出现顺序”排列的对象名 —— 与 spclust 的 init_grp() 一致，
    差分比较 merge 序列时依赖这个顺序。

    解析规则与 readdb.cpp 对齐：支持 '#' 注释、逗号或空白分隔、可选表头、
    丢弃自环。重复 pair 保留，但矩阵语义上“后出现的覆盖先出现的”。
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
            warnings.append(f"{path}:{lineno}: 字段不足 3 个，跳过")
            continue

        a, b, s = fields[0], fields[1], fields[2]
        if not a or not b:
            warnings.append(f"{path}:{lineno}: 对象名为空，跳过")
            continue
        v = to_fraction(s)
        if v is None:
            # 只容忍开头的表头行
            if not edges and [x.lower() for x in fields[:3]] == ["obj1", "obj2", "sim"]:
                continue
            warnings.append(f"{path}:{lineno}: 非法 sim 值 {s!r}，跳过")
            continue
        if not 0 <= v <= 65535:
            warnings.append(f"{path}:{lineno}: sim={v} 越界，跳过")
            continue
        if a == b:
            warnings.append(f"{path}:{lineno}: 自环 {a!r}，跳过")
            continue

        for x in (a, b):
            if x not in seen:
                seen.add(x)
                order.append(x)
        edges.append((a, b, v))

    return order, edges, warnings


def _build_dense(order, edges):
    """原始相似度矩阵（brute 模式用），缺失为 0。返回 (S, duplicate_count)。

    同一个 pair 出现多次且值不同时**取最大值** —— 这是 spclust 的实际语义
    （它的合并值来自 get_max()，即行内最大值），实测两组相反顺序的输入结果一致。
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
        if v > S[i][j]:                     # 重复 pair 取最大值（不是"后出现的覆盖"）
            S[i][j] = S[j][i] = v
    return S, dup


def _canon(a, b):
    """把一对簇规范化成固定顺序，便于逐字比较。"""
    return (a, b) if min(a) < min(b) else (b, a)


def hac(order, edges, linkage="average", min_sim=50, mode="brute"):
    """跑一遍参考 HAC，返回结果字典。

    mode='brute' —— 每步从原始矩阵按定义重算（最可信，O(n^4)）
    mode='fast'  —— 增量维护精确聚合量（O(n^3)）
    mode='trunc' —— 镜像 spclust 的整数语义：每步把链接值截断成整数再往下用。
                    这不是独立 oracle，而是个“假设检验”：如果 spclust 与精确定义的
                    偏差 *完全* 来自 unsigned short 截断，那它与本模式必须逐次一致。
    """
    if linkage not in LINKAGES:
        raise ValueError(f"unknown linkage {linkage!r}, expect one of {LINKAGES}")
    if mode not in ("brute", "fast", "trunc"):
        raise ValueError(f"unknown mode {mode!r}")

    n = len(order)
    S, duplicate_pairs = _build_dense(order, edges)

    members = {i: [i] for i in range(n)}
    live = list(range(n))

    # fast / trunc 模式的累积量：(i,j) -> 聚合值
    #   fast  : average 存 sum，max/min 存 max/min
    #   trunc : 直接存 spclust 那种（已截断的）整数链接值
    acc = {}
    if mode in ("fast", "trunc"):
        for i in range(n):
            for j in range(i + 1, n):
                acc[(i, j)] = S[i][j]

    def value(i, j):
        """返回簇 i 与 j 的链接相似度，精确分数。"""
        A, B = members[i], members[j]
        if mode == "brute":
            if linkage == "average":
                # 分母是全部 |A|*|B| 对（缺失的 pair 贡献 0）
                return Fraction(sum(S[a][b] for a in A for b in B), len(A) * len(B))
            if linkage == "max":
                # 缺失 pair 算 0，因此下界是 0
                return Fraction(max(S[a][b] for a in A for b in B), 1)
            # min：缺失 pair 算 0，因此结果可能为 0（= 非候选）
            return Fraction(min(S[a][b] for a in A for b in B), 1)

        key = (i, j) if i < j else (j, i)
        a = acc[key]
        if mode == "fast" and linkage == "average":
            # 每步重新按定义除以 |A|*|B|
            return Fraction(a, len(members[i]) * len(members[j]))
        return Fraction(a, 1)

    merges = []
    tie_steps, trunc_steps = set(), set()
    step = 0
    next_id = n

    while len(live) > 1:
        # 用 max() 对 (Fraction, i, j) 取 argmax：显然不会写错，
        # 也不要手工写交叉相乘的比较（那正是最容易把簇号和相似度搞混的地方）。
        cand = []
        for x in range(len(live)):
            for y in range(x + 1, len(live)):
                i, j = live[x], live[y]
                v = value(i, j)
                if linkage == "min" and not (v > 0):
                    continue  # min 只在结果 > 0 时才是候选
                cand.append((v, i, j))
        if not cand:
            break

        best_val, i, j = max(cand)
        #相似度 <= 0 就不再合并 —— 与 spclust 的 get_max()（相似度为 0 不算候选）对齐。
        #于是 min_sim=0 得到的就是**完整合并树**；min_sim>0 只是提前停（等价于事后切树）。
        if best_val < min_sim or best_val <= 0:
            break

        best_floor = best_val.numerator // best_val.denominator
        # 并列最大？（会导致合并顺序不确定，两边可能合法地走不同分支）
        if sum(1 for v, _, _ in cand if v == best_val) > 1:
            tie_steps.add(step)
        # 截断后并列？（spclust 存 unsigned short，会把 88.5 和 88.4 压成同一个值）
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
                        # 镜像 spclust:
                        #   sim = (sim1*|G1| + sim2*|G2|) / (|G1|+|G2|)   ← 整数除法
                        # 再赋给 unsigned short，于是每步都会截断一次
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
    # partition 用对象下标表示，与 parse_spclust() 的输出可直接比较
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
    """把完整合并树切在 min_sim 处，返回扁平分组（对象下标的 frozenset 的 frozenset）。

    相似度沿合并序列单调不增，所以"切"就是取 sim >= min_sim 的那段前缀；
    每取一次合并就把两个子簇并起来。这与 spclust 的做法一致（它用并查集，
    这里直接合并成员集合，更直白）。
    传进来 result 必须是 min_sim=0 算出的**完整树**。
    """
    live = set(frozenset([i]) for i in range(result["objects"]))
    for (A, B), v in result["merges"]:
        if v < min_sim:
            break                      # 前缀结束（序列单调不增）
        live.discard(A)
        live.discard(B)
        live.add(A | B)
    return frozenset(live)


# ---------------------------------------------------------------------------
# 解析 spclust 的输出
# ---------------------------------------------------------------------------

RE_GROUP = re.compile(r"^Group (\d+) \(size=(\d+), sim=([0-9.eE+-]+)\): (.*)$")
RE_SINGLE = re.compile(r"^Unclustered \((\d+)\): (.*)$")
RE_OBJECTS = re.compile(r"^Objects: (\d+), groups: (\d+)$")


def parse_spclust(text, order):
    """把 CLI 的输出解析成和 hac() 同构的结构。

    合并序列来自 **--matrix**（公开的 linkage 矩阵：left right height size），
    节点 id 沿用 scipy 约定：0..n-1 是对象、n+i 是第 i 次合并 —— 于是逐行重建
    出每个节点的成员集合即可。扁平分组来自末尾那几行。

    以前这里解析的是 DEBUG=2 的合并日志，那让测试与一条调试输出格式隐性耦合
    （改一行输出就可能让 runner 静默地少算），现在改成走公开接口。

    order 是对象的首次出现顺序 —— CLI 就是按这个顺序分配 id 的，
    所以叶子 id i 对应 order[i]。
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
        raise ValueError("没找到 'Objects: N, groups: K' 行，CLI 输出格式可能变了")

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

    ap = argparse.ArgumentParser(description="spclust 的独立参考实现")
    ap.add_argument("csv")
    ap.add_argument("--linkage", default="average", choices=LINKAGES)
    ap.add_argument("--min-sim", type=int, default=50)
    ap.add_argument("--mode", default="brute", choices=("brute", "fast", "trunc"))
    ap.add_argument("--cut", type=float, default=None,
                    help="把结果按这个阈值切开（默认与 --min-sim 相同）")
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
    print(f"oracle({args.linkage}/{args.mode}, 完整树 {r['merge_count']} 次合并, "
          f"切在 {cut_at}): {r['objects']} objects, {len(groups)} groups")
    for g in groups:
        if len(g) >= 2:
            print(f"  size={len(g):<4} {g}")
    singles = sorted(x for g in groups if len(g) == 1 for x in g)
    if singles:
        print(f"  Unclustered ({len(singles)}): {singles}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

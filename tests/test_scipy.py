#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
#
# scipy 互操作测试：把 `spclust_cli --scipy` 的导出交给**真的 scipy** 去读。
#
# 这一层存在的理由：前面所有测试都用我们自己写的 oracle 当尺子。这里换成 scipy ——
# 一个完全不同来源、被广泛使用的实现。它验证两件事：
#   1. 格式真的合规（scipy 的 is_valid_linkage / to_tree / fcluster 能直接用）；
#   2. 语义真的对得上（在**没有并列**的用例上，逐行与 scipy 自己的 linkage 完全相同）。
#
# 语义对应的关键：本库的"缺失边 = 相似度 0"正好对应"距离 = C - 相似度"（C = 最大相似度），
# 于是 average 链接等价于 scipy 的 method='average' 作用在同一个稠密距离矩阵上。
#
# 关于并列：相似度相等时合并顺序本就不作保证（契约里写了）。注意**输入值全不同也会出现并列**
# ——派生出来的中间相似度会撞在一起（实测 15_many_merges.csv 就是这样）。所以逐行比对
# 只在"该用例该 linkage 没有并列"时做（并列与否问 oracle），否则只查格式合规。
#
# 运行：python3 tests/test_scipy.py --bin ./spclust [--venv /tmp/bvenv/bin/python]
# 没有 scipy 时整体跳过（返回 0），这样 CI 上不装 scipy 也能过。

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
import oracle  # noqa: E402

CASES = sorted(str(p) for p in (HERE / "cases").glob("*.csv")) + [str(REPO / "sim.csv")]

# 本库的 linkage -> scipy 的 method。注意第三个是**故意不同的**：
#   Linkage::MinNonzero 只对"两边都存在的边"取最小；scipy 的 complete 在缺失对上取到
#   最大距离（= 相似度 0）。稠密输入下两者相同（已实测），稀疏下不同。
# (本库 linkage 编号, oracle 里的名字, scipy 的 method 名)
MATCHING = [(0, "average", "average"), (1, "max", "single")]
KNOWN_DIFFERENT = (2, "complete")


def parse_matrix(text):
    """解析 --scipy 的输出（# 开头是注释；允许 inf）。"""
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
        raise RuntimeError(f"spclust 退出码 {proc.returncode}:\n{proc.stderr}")
    return parse_matrix(proc.stdout), proc.stderr


def dense_from_case(case):
    """(对象名顺序, 稠密相似度矩阵, 最大相似度)。缺失 pair 为 0 —— 与本库语义一致。"""
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
    """树不是完整树（合并次数 < n-1）就不做逐行比对 —— 桥接行的顺序是任意的。"""
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
    ap = argparse.ArgumentParser(description="scipy 互操作测试")
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
        print(f"跳过 scipy 互操作测试：{e}")
        print("（装上 numpy/scipy 再跑，例如 python3 -m venv /tmp/v && /tmp/v/bin/pip install numpy scipy）")
        return 0

    print("=" * 78)
    print(f"spclust x scipy 互操作    {bin_path}    scipy {__import__('scipy').__version__}"
          f"    用例 {len(CASES)} 个")
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

            # ① 格式合规：scipy 自己要认
            if not is_valid_linkage(Z, throw=False, warning=False) or not is_monotonic(Z):
                bad.append(f"{name} linkage={mtd}: scipy 认为是非法 linkage 矩阵")
                continue
            try:
                to_tree(Z)
                fcluster(Z, 2, criterion='maxclust')
            except Exception as e:                                    # pragma: no cover
                bad.append(f"{name} linkage={mtd}: scipy 下游函数失败 {e}")
                continue
            n_valid += 1

            # ② 语义对得上：**树是完整树**（我们的合并次数 = n-1）且没有并列时，
            # 应当逐行与 scipy 完全相同。森林（有分量永远合不起来，桥接行是补的）
            # 不在这个范围内 —— 桥接的距离语义正确，但分量的合并顺序是任意的。
            Z_sc = scipy_linkage(cond, method=method)
            if len(Z) != len(Z_sc):
                bad.append(f"{name} linkage={mtd}: 行数 {len(Z)} != scipy {len(Z_sc)}")
                continue
            if ties:
                n_tie += 1
                if args.verbose:
                    print(f"  (并列，跳过逐行比对) {name} linkage={mtd}")
                continue
            if n_forest_skip(case, oracle_name, n):
                if args.verbose:
                    print(f"  (森林，跳过逐行比对) {name} linkage={mtd}")
                continue
            same = True
            for i in range(len(Z)):
                if {int(Z[i][0]), int(Z[i][1])} != {int(Z_sc[i][0]), int(Z_sc[i][1])} \
                   or abs(Z[i][2] - Z_sc[i][2]) > 1e-9 or int(Z[i][3]) != int(Z_sc[i][3]):
                    bad.append(f"{name} linkage={mtd}: 第 {i} 行不同 "
                               f"{Z[i].tolist()} vs {Z_sc[i].tolist()}")
                    same = False
                    break
            if same:
                n_exact += 1

    # ③ 森林：残缺矩阵对 scipy 来说只是"另一个较小的 n"（它没能力发现少了几行），
    #    所以这里查的是真正的保证：必须警告，且 --scipy-bridge 之后行数恰好 n-1、
    #    最后一行的 size 恰好是 n。
    forest_checked = 0
    for case in CASES:
        order, edges, _ = oracle.load_edges(case)
        n = len(order)
        if n < 2:
            continue
        res = oracle.hac(order, edges, "min", 0, "brute" if n <= 24 else "fast")
        if len(res["merges"]) == n - 1:
            continue                      # 树是完整的，这条检查不适用
        forest_checked += 1
        plain, stderr_plain = run_scipy(bin_path, case, 2, ())
        if "森林" not in stderr_plain:
            bad.append(f"{Path(case).name}: 残缺树没有在 stderr 上警告")
        bridged, _ = run_scipy(bin_path, case, 2, ("--scipy-bridge",))
        if len(bridged) != n - 1:
            bad.append(f"{Path(case).name}: 桥接后应有 n-1={n-1} 行，实际 {len(bridged)}")
        elif int(bridged[-1][3]) != n:
            bad.append(f"{Path(case).name}: 桥接后最后一行的 size 应是 {n}，"
                       f"实际 {int(bridged[-1][3])}")
        if len(plain) >= len(bridged):
            bad.append(f"{Path(case).name}: 桥接前后行数没有变化")
        if forest_checked >= 4:
            break
    if forest_checked == 0:
        print("  （没有森林用例，跳过桥接检查）")
    else:
        print(f"  森林用例 {forest_checked} 个：都警告了，且桥接后行数 = n-1")
    # ④ 已知不同：MinNonzero vs scipy complete。
    #    稠密矩阵上两者必须相同；稀疏矩阵上允许不同（并应当真的不同，否则这条就没意义）。
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
        # 稠密 = 所有 pair 都有边（上三角无 0）
        dense = bool((np.triu(S) > 0).sum() == n * (n - 1) // 2)
        if dense:
            dense_total += 1
            dense_same += 1 if same else 0
        else:
            sparse_total += 1
            sparse_same += 1 if same else 0
    if dense_total and dense_same != dense_total:
        bad.append(f"MinNonzero 与 scipy complete 在稠密输入上应当相同，"
                   f"实际 {dense_same}/{dense_total}")
    print(f"  稠密输入上 MinNonzero == scipy complete: {dense_same}/{dense_total}")
    print(f"  稀疏输入上 MinNonzero == scipy complete: {sparse_same}/{sparse_total}"
          f"（本来就该不同：scipy 在缺失对上取到最大距离）")

    print("-" * 78)
    print(f"格式合规 {n_valid} 组；逐行与 scipy 完全相同 {n_exact} 组；"
          f"因并列跳过逐行比对 {n_tie} 组")
    if bad:
        print("失败明细:")
        for b in bad[:12]:
            print("  ✗", b)
    print(f"RESULT: {'PASS' if not bad else 'FAIL'}")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())

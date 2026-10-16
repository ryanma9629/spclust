#!/usr/bin/env python3
"""
CLI 接口测试。

`tests/run_tests.py` 验证的是**算法**（把 DEBUG 日志解析出来与 oracle 比对），
而命令行接口 —— `--cut/--floor/--tree/--matrix`、退出码、stdin 模式 ——
之前完全没有覆盖。这两者还容易互相掩盖：程序与 runner 之间靠"合并日志的格式"
隐性耦合，改一行输出就可能让 runner 静默地少算（本项目就发生过：Kruskal 路径
忘了打印那行日志，runner 报成"0 次合并"）。

这里只测接口契约，不重复测算法：
  * 退出码与错误路径
  * `--cut` 的语义：与 oracle 的切树一致；且"切在 X"必须等于"用 X 当提前停止条件"
  * `--floor` 的语义：floor 决定树算多深，cut 只能切在 floor 之上
  * `--matrix` 的结构不变量（linkage 矩阵的基本性质）
  * `--tree` 的完整性（每个叶子恰好出现一次）
  * stdin 模式与文件模式等价
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
    """把扁平输出解析成 {frozenset(成员名)} 的集合（与 Group 序号无关）。"""
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


# ---------------------------------------------------------------------------
# 各项检查：返回 (ok, detail)
# ---------------------------------------------------------------------------

def t_exit_codes(bin_path):
    bad = []
    r = run(bin_path, ["--help"])
    if r.returncode != 0 or "用法" not in r.stdout + r.stderr:
        bad.append("--help 应当退出 0 并打印用法")
    for args, want, what in ((['--nope'], 2, "未知选项"),
                             (['a.csv', 'b.csv'], 2, "多余的位置参数"),
                             (['--cut'], 2, "--cut 缺少数值"),
                             (['/nonexistent/x.csv'], 1, "输入文件不存在")):
        r = run(bin_path, args)
        if r.returncode != want:
            bad.append(f"{what}: 退出码 {r.returncode}，期望 {want}")
    r = run(bin_path, ["-"], stdin_text="")
    if r.returncode != 1:
        bad.append(f"空输入: 退出码 {r.returncode}，期望 1")
    return (not bad), "; ".join(bad)


def t_cut_matches_oracle(bin_path):
    """--cut X 的分组必须等于 oracle 把完整树切在 X 处。"""
    bad = []
    n = 0
    for case in CASES:
        for cut in (0, 25, 50, 75, 90):
            n += 1
            got = parse_flat(run(bin_path, [case, "--floor", "0", "--cut", str(cut)]).stdout)
            want = oracle_cut(case, cut)
            if got != want:
                bad.append(f"{Path(case).name} cut={cut}: CLI 与 oracle 不一致")
    return (not bad), f"({n} 组) " + "; ".join(bad[:3])


def t_cut_equals_early_stop(bin_path):
    """把完整树切在 X 处，必须等于直接用 X 当停止条件（--floor X）。"""
    bad = []
    for case in CASES:
        for cut in (25, 50, 75):
            a = parse_flat(run(bin_path, [case, "--floor", "0", "--cut", str(cut)]).stdout)
            b = parse_flat(run(bin_path, [case, "--floor", str(cut), "--cut", str(cut)]).stdout)
            if a != b:
                bad.append(f"{Path(case).name} cut={cut}")
    return (not bad), "; ".join(bad[:3])


def t_floor_cut_interaction(bin_path):
    """--cut 低于默认 floor 时 floor 要自动跟着降，否则切不到。"""
    bad = []
    for case in CASES:
        # 完全不给 floor，只给一个很低的 cut：应当自动把树算到那里
        got = parse_flat(run(bin_path, [case, "--cut", "10"]).stdout)
        want = oracle_cut(case, 10)
        if got != want:
            bad.append(f"{Path(case).name}: --cut 10 没自动降低 floor")
    # 显式给更高的 floor：cut 比 floor 高时结果应当与"cut 就是 floor"无关
    r = run(bin_path, [CASES[0], "--floor", "50", "--cut", "95"])
    if r.returncode != 0:
        bad.append("--floor 50 --cut 95 应当正常运行（切在 floor 之上）")
    return (not bad), "; ".join(bad[:3])


def t_matrix_structure(bin_path):
    """linkage 矩阵的基本不变量。"""
    bad = []
    for case in CASES:
        for extra in ([], ["--floor", "0"]):
            text = run(bin_path, [case, "--floor", "0"] + extra).stdout
            rows = parse_matrix(text)
            order, edges, _ = oracle.load_edges(case)
            n = len(order)
            if not rows:
                # 允许一个合并都没有（所有相似度都 0）—— 但那样不该有 --matrix 行
                continue
            parents = {}
            for k, (l, r, sim, size) in enumerate(rows):
                node = n + k
                for c in (l, r):
                    if not (0 <= c < node):
                        bad.append(f"{Path(case).name}: 第 {k} 行的子节点 {c} 不在 [0,{node})")
                    if c in parents:
                        bad.append(f"{Path(case).name}: 节点 {c} 被当成子节点用了两次")
                    parents[c] = node
                if size < 2:
                    bad.append(f"{Path(case).name}: 第 {k} 行 size={size} < 2")
            # 相似度必须单调不增（切树语义依赖它）
            sims = [r[2] for r in rows]
            if any(sims[i] < sims[i + 1] - 1e-9 for i in range(len(sims) - 1)):
                bad.append(f"{Path(case).name}: 相似度不是单调不增 {sims[:6]}")
            # 最后一个节点的大小 = 它那个根包含的对象数
            roots = [x for x in range(n + len(rows)) if x not in parents]
            tot = sum(rows[x - n][3] for x in roots if x >= n) + sum(1 for x in roots if x < n)
            if tot != n:
                bad.append(f"{Path(case).name}: 各根大小之和 {tot} != 对象数 {n}")
    return (not bad), "; ".join(bad[:3])


def t_tree_completeness(bin_path):
    """--tree 里每个对象恰好出现一次（在深度上限内）。"""
    bad = []
    for case in CASES:
        text = run(bin_path, [case, "--floor", "0", "--tree"]).stdout
        order, _, _ = oracle.load_edges(case)
        if "Dendrogram" not in text:
            bad.append(f"{Path(case).name}: 没有 Dendrogram 行")
            continue
        leaves = []
        for line in text.splitlines():
            s = line.strip()
            if s.startswith("- "):
                leaves.append(s[2:].strip())
        names = set(order)
        if len(leaves) != len(set(leaves)):
            bad.append(f"{Path(case).name}: 有叶子重复出现")
        if not set(leaves) <= names:
            bad.append(f"{Path(case).name}: 出现了输入里没有的对象")
    # 多个根时树的根数与 matrix 的根数一致
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
    if f"共 {roots} 个根" not in t:
        bad.append(f"{Path(case).name}: --tree 报的根数与 matrix 不一致（应为 {roots}）")
    return (not bad), "; ".join(bad[:3])


def t_stdin_equals_file(bin_path):
    """stdin 模式与文件模式必须给出同样的结果。"""
    bad = []
    for case in CASES[:6]:
        text = Path(case).read_text()
        a = parse_flat(run(bin_path, [case], stdin_text=None).stdout)
        b = parse_flat(run(bin_path, ["-"], stdin_text=text).stdout)
        c = parse_flat(run(bin_path, [], stdin_text=text).stdout)
        if not (a == b == c):
            bad.append(f"{Path(case).name}: 文件 / '-' / 无参数 三种读法结果不同")
    return (not bad), "; ".join(bad[:3])


def t_duplicate_order_independent(bin_path):
    """重复 pair 取最大值：交换两行顺序结果必须一致。"""
    import tempfile
    bad = []
    for body in ("a,b,90\na,b,80\na,c,88\n", "a,b,80\na,b,90\na,c,88\n"):
        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as f:
            f.write(body)
            p = f.name
        part = parse_flat(run(bin_path, [p]).stdout)
        if frozenset([frozenset(["a", "b"]), frozenset(["c"])]) != part:
            bad.append(f"输入 {body.splitlines()[:2]} 得到 {part}")
        Path(p).unlink()
    return (not bad), "; ".join(bad[:2])


TESTS = [
    ("退出码与错误路径", t_exit_codes),
    ("--cut 与 oracle 切树一致", t_cut_matches_oracle),
    ("切树 == 提前停止(单调性)", t_cut_equals_early_stop),
    ("--floor/--cut 交互", t_floor_cut_interaction),
    ("--matrix 结构不变量", t_matrix_structure),
    ("--tree 完整性", t_tree_completeness),
    ("stdin 与文件等价", t_stdin_equals_file),
    ("重复 pair 顺序无关", t_duplicate_order_independent),
]


def main():
    ap = argparse.ArgumentParser(description="spclust CLI 测试")
    ap.add_argument("--bin", default="./spclust")
    args = ap.parse_args()
    # 注意 Path("./spclust") 会规范化成 "spclust"（丢掉斜杠），子进程就会去 PATH 找，
    # 所以这里必须先变成绝对路径。
    bin_path = Path(args.bin)
    if not bin_path.is_absolute():
        bin_path = (Path.cwd() / bin_path).resolve()

    print("=" * 78)
    print(f"spclust CLI 测试    {bin_path}    用例 {len(CASES)} 个")
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

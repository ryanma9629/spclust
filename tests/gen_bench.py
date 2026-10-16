#!/usr/bin/env python3
"""
基准输入的**固定生成器**。

为什么必须固定：运行时强烈依赖相似度分布 —— 同样的对象数和边数，
值域取 0..100（大部分边低于阈值、合并少）与取 1..65535（几乎全部参与合并）
实测能差 3 倍以上。换一个随机种子或值域，前后两次的数字就不可比了。

所以：参数写死在下面，`tests/bench_baseline.txt` 里的数字只在参数不变时可比。
要改参数就必须重新生成基线（`python3 tests/bench.py --update`）。

用法:
    python3 tests/gen_bench.py <case> > /tmp/case.csv     # 生成某个用例
    python3 tests/gen_bench.py --list                      # 列出全部用例
    python3 tests/gen_bench.py --all <dir>                 # 全部生成到目录
"""

import os
import random
import sys

# (名字, 对象数n, 平均度deg, sim 最小值, sim 最大值, 种子)
# 边数 ≈ n * deg / 2。值域上界受读取端限制：sim ∈ [0, 65535]。
CASES = [
    ("bench_small",  2000, 20,     1, 65535, 1),   # 快，用来做多次重复
    ("bench_mid",   20000, 20,     1, 65535, 2),   # 大多数边都参与合并（重活）
    ("bench_large", 50000, 20,     1, 65535, 3),   # 大一点，看规模趋势
    ("bench_lowval", 20000, 20,    0,   100, 4),   # 大多数边低于阈值（轻活）
    ("bench_single", 50000, 20,    1, 65535, 5),   # 单链接（Kruskal 路径）用
]


def case_params(name):
    for c in CASES:
        if c[0] == name:
            return c
    raise SystemExit(f"未知用例 {name!r}；用 --list 看可用的")


def generate(name):
    """生成边表文本（对象名 o<i>，pair 按字典序排序，保证逐位可复现）。"""
    _, n, deg, lo, hi, seed = case_params(name)
    rng = random.Random(seed)
    target = n * deg // 2
    pairs = set()
    # 用集合去重抽 pair；边数远小于 n^2 时这样最快也最稳
    while len(pairs) < target:
        a = rng.randrange(n)
        b = rng.randrange(n)
        if a == b:
            continue
        if a > b:
            a, b = b, a
        pairs.add((a, b))
    # 值在生成 pair 之后再抽，这样"哪条边拿到哪个值"也是确定的
    lines = []
    for i, j in sorted(pairs):
        lines.append(f"o{i},o{j},{rng.randint(lo, hi)}")
    return "\n".join(lines) + "\n"


def main(argv):
    if len(argv) < 2 or argv[1] in ("-h", "--help"):
        print(__doc__)
        return 0
    if argv[1] == "--list":
        for name, n, deg, lo, hi, seed in CASES:
            print(f"{name:<14} n={n:<7} deg={deg:<3} sim∈[{lo},{hi}] seed={seed:<3} 边数≈{n * deg // 2}")
        return 0
    if argv[1] == "--all":
        outdir = argv[2]
        os.makedirs(outdir, exist_ok=True)
        for name, *_ in CASES:
            path = os.path.join(outdir, name + ".csv")
            with open(path, "w") as f:
                f.write(generate(name))
            print(f"wrote {path} ({os.path.getsize(path)} bytes)")
        return 0
    sys.stdout.write(generate(argv[1]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

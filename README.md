# spclust

基于**稀疏相似度矩阵**的凝聚型层次聚类（agglomerative hierarchical clustering）。

反复取出全局相似度最大的一对簇并合并，直到最大相似度低于给定下限；用**十字链表**
存稀疏上三角矩阵，因此内存是 O(nnz) 而不是 O(n²) —— 这是它相对稠密实现的主要价值。
合并树被完整保留，阈值只是"切树位置"，换阈值不需要重新聚类。

```
库：      include/spclust/hac.h + src/hac.cpp      不依赖任何第三方库，只认整数 id
使用者：  tools/                                    所有 IO 与展示（读 CSV、id↔名字、打印）
```

## 快速开始

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure      # 差分测试 + CLI 契约测试
./build/spclust_cli sim.csv                     # 跑示例数据
```

只用标准库，没有第三方依赖：一个 C++17 编译器即可。也可以用 Makefile（开发时的快捷方式）：

```bash
make            # 构建 ./spclust
make lib        # 只把库编成 libspclust.a
make run        # 跑 sim.csv
make test       # 差分测试 + CLI 测试
make bench      # 性能回归护栏（固定输入 + 记录的基线）
make check      # 三种 linkage 各跑一遍做对比
make asan       # AddressSanitizer + UBSan
make asan-cases # 用 ASan/UBSan 跑遍 tests/cases 全部用例
make leaks      # 查内存泄漏
make clean
```

## 作为库使用

```cmake
find_package(spclust REQUIRED)
target_link_libraries(your_target PRIVATE spclust::spclust)
```

```cpp
#include <spclust/hac.h>

spclust::reserve(n_objects, n_edges);          // 可选：预留容量
for (auto& e : edges)
    spclust::add_edge(e.a, e.b, e.similarity); // id 必须是稠密的 0..n-1
spclust::build(0.0);                           // 0 = 算完整棵树

for (const spclust::Merge& m : spclust::tree())
    ...                                        // 就是 linkage 矩阵：left right similarity size

spclust::Cut c = spclust::cut(50.0);            // 把树切在相似度 50 处（不需要重新聚类）
for (const auto& members : c.members)
    ...                                        // 每个簇的对象 id
spclust::reset();                              // 同一进程里接着跑下一次
```

**用之前请读 `include/spclust/hac.h` 顶部的契约**（缺失边视为 0、id 必须稠密、
重复 pair 取最大值、自环抛异常、并列时合并顺序不作保证、以及"目前一个进程一份状态"）。

## 输入格式

每行 `obj1,obj2,sim`，也可以从标准输入读（`spclust_cli < data.csv` 或 `spclust_cli -`）。

- `sim`：CLI 接受 0..65535 的实数（允许小数，如 `87.5`）。**上界是 CLI 的输入契约**
  （历史遗留）；库本身接受任意有限非负 `double`。低于切树阈值的边只作为候选，不触发合并。
- 分隔符：含逗号按逗号拆，否则按空白（空格 / Tab）拆
- `#` 开头是注释；空行忽略；可选的表头行 `obj1,obj2,sim` 自动跳过
- 越界 / 非数字的 `sim`、字段数不足、自环（`obj1 == obj2`）都会打印警告并跳过该行
- **重复 pair 取最大值**（语义已定义并被测试钉住）

## 输出

默认打印"把合并树在 `--cut`（默认 50）处切开"的扁平分组：

```
Objects: 16, groups: 3
Group 0 (size=5, sim=90.5): a1 a2 a3 a4 a5
Group 1 (size=4, sim=87.6667): b1 b2 b3 b4
Group 2 (size=3, sim=82): c1 c2 c3
Unclustered (4): d1 d2 e1 e2
```

`sim` 是组内平均相似度（average 链接下按组大小加权累加），可以当作簇的紧致度；
只列出 size ≥ 2 的簇，其余单元素对象列在 `Unclustered` 里，成员按 id 顺序（输入中首次出现的顺序）。

因为合并树是完整算出来的，阈值只是"切在哪里"：

```
spclust_cli data.csv --cut 86         # 切在 86（不重新聚类）
spclust_cli data.csv --tree           # 树状图（--tree-depth N 控制深度，0 = 不限）
spclust_cli data.csv --matrix         # linkage 矩阵：left right height size
spclust_cli data.csv --floor 0        # 把树算到最底部（完整树；默认 = cut，代价更低）
```

节点 id 沿用 scipy 的约定：`0..n-1` 是对象、`n+i` 是第 i 次合并 —— 这里组号天然就是节点 id。
**注意 `height` 是相似度（越大越先合并），scipy 用的是距离。**
`--cut` / `--floor` 的语义与成本见 [docs/design.md](docs/design.md)。

## 文档

| 文件 | 内容 |
|---|---|
| [docs/design.md](docs/design.md) | 数据结构、算法流程、可调参数 |
| [docs/known-issues.md](docs/known-issues.md) | 修过的 bug（含发现过程）+ 仍存在的问题 |
| [docs/performance.md](docs/performance.md) | 性能实测、被实测否掉的优化、方法学教训 |
| [docs/testing.md](docs/testing.md) | 差分测试 / CLI 测试 / 性能护栏的设计 |
| [CHANGELOG.md](CHANGELOG.md) | 对使用者有意义的变化 |

## 测试与性能的量级

- 差分测试：独立参考实现（不含 Lance-Williams 递推）逐次合并比对 + 22 个手写边界用例；
  另有 900 次多阈值切树交叉验证。当前 `PASS=90 WARN=12 KNOWN=0 FAIL=0`
  （WARN 是并列导致的合并顺序不确定，两边都合法）。
- ASan/UBSan 干净、0 泄漏；节点池手工毒化，所以"读已释放节点"仍能被 ASan 抓到。
- 累计性能（对比最开始的基线，详见 [docs/performance.md](docs/performance.md)）：
  800k 边 **21.9x**、2M 边 **51.1x**，峰值内存 **−55~61%**；单链接 **300x**。

## 许可

MIT，见 [LICENSE](LICENSE)。

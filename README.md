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
make test       # 差分测试 + 切树属性测试 + CLI 测试
make test-cut   # 只跑切树属性测试
make test-scipy # 只跑 scipy 互操作（没装 scipy 会自动跳过）
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

spclust::HAC hac(spclust::Linkage::Average);   // 链接方式是**运行时**参数

hac.reserve(n_objects, n_edges);               // 可选：预留容量
for (auto& e : edges)
    hac.add_edge(e.a, e.b, e.similarity);      // id 必须是稠密的 0..n-1

hac.build(0.0);                                // 0 = 算完整棵树

for (const spclust::Merge& m : hac.tree())
    ...                                        // 就是 linkage 矩阵：left right similarity size

// 三种切法，任意切换，都不需要重新聚类
spclust::Cut c = hac.cut(50.0);                // 切在相似度 50 处
spclust::Cut k = hac.cut_n_clusters(12);       // 要正好 12 个簇
spclust::Cut s = hac.cut_min_size(5);          // 合并到没有小于 5 的簇为止

for (const auto& members : c.members)
    ...                                        // 每个簇的对象 id（Cut 里还有 cohesion）

hac.reset();                                   // 同一个对象接着跑下一次

// 一个进程里可以有任意多个独立聚类，互不影响，也不需要加锁。
// （只要别让两个线程用**同一个** HAC 对象。）
spclust::HAC other(spclust::Linkage::MaxSimilarity);
```

**用之前请读 `include/spclust/hac.h` 顶部的契约**（缺失边视为 0、id 必须稠密、
重复 pair 取最大值、自环抛异常、并列时合并顺序不作保证、一个对象持有独立状态且
不与读接口并发安全）。链接方式对簇形的影响见 [docs/design.md](docs/design.md#链接方式)。

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
spclust_cli data.csv --n-clusters 12  # 要正好 12 个簇（需要完整树）
spclust_cli data.csv --min-size 5     # 合并到没有小于 5 的簇为止（需要完整树）
spclust_cli data.csv --linkage 1      # 链接方式：0/average 1/max|single 2/min|min-nonzero
spclust_cli data.csv --tree           # 树状图（--tree-depth N 控制深度，0 = 不限）
spclust_cli data.csv --matrix         # linkage 矩阵：left right height size
spclust_cli data.csv --floor 0        # 把树算到最底部（完整树；默认 = cut，代价更低）
```

节点 id 沿用 scipy 的约定：`0..n-1` 是对象、`n+i` 是第 i 次合并 —— 这里组号天然就是节点 id。
**注意 `height` 是相似度（越大越先合并），scipy 用的是距离。**
`--cut` / `--floor` 的语义与成本见 [docs/design.md](docs/design.md)。

## C 接口（给 FFI 用）

`include/spclust/spclust_c.h` 是**纯 C**头文件（只用 `<stddef.h>` / `<stdint.h>`），
能被 C / Rust / Go / Python ctypes 直接 include 或绑定：

```c
#include "spclust/spclust_c.h"

spclust_hac* h = spclust_hac_new(SPCLUST_LINKAGE_AVERAGE);
spclust_hac_add_edge(h, 0, 1, 90.0);
spclust_hac_build(h, 0.0);

size_t n = spclust_hac_tree_size(h);
spclust_merge* rows = malloc(n * sizeof *rows);
spclust_hac_tree(h, rows, n, &n);          /* 两段式：先问长度，再填缓冲 */

spclust_cut* c = spclust_hac_cut(h, 50.0);
for (size_t g = 0; g < spclust_cut_n_groups(c); g++) { ... }
spclust_cut_free(c);
spclust_hac_free(h);
```

设计要点：

- **C++ 异常绝不穿过 C 边界**。每个入口都 `try/catch`，转成 `spclust_status`
  （`SPCLUST_ERR_INVALID_ARGUMENT` / `OUT_OF_MEMORY` / `BUFFER_TOO_SMALL` / ...）。
  细节用 `spclust_last_error()` 取，它是**线程局部**的。
  （`tests/test_c_api.py` 专门测这条：从 Python 传自环，必须拿到错误码，而不是把进程 abort 掉。）
- **不透明句柄**（`spclust_hac*` / `spclust_cut*`），一个句柄一份状态，多个句柄互不影响。
- **缓存的归属清楚**：`spclust_hac_tree` 写调用方给的缓冲；`spclust_cut` 的一堆扁平视图
  （`members` / `offsets` / `cohesion`）归库所有，用 `spclust_cut_free` 一次释放。
- 链接方式、三种切法、scipy 导出都有对应入口；`spclust_linkage_from_name` 认
  `average` / `max-similarity` / `min-nonzero`（以及 `max`/`single`/`min` 别名）。

```bash
make lib        # libspclust.a
make lib-shared # libspclust.so（FFI 用；macOS 上 -shared 出来的也是 Mach-O，ctypes 能直接 dlopen）
make test-c     # 用 C 编译器编测试 + 从 Python 用 ctypes 调一遍
```

## 与 scipy 互操作

```bash
spclust_cli data.csv --scipy --scipy-bridge > Z.txt
```
```python
from scipy.cluster.hierarchy import fcluster, dendrogram, to_tree
import numpy as np
Z = np.loadtxt("Z.txt")          # 直接就是 scipy 的 linkage 矩阵
fcluster(Z, 3, criterion='maxclust')
```

`height` 在 scipy 里是**距离**（越小越先合并），本库内部是**相似度**（越大越先合并），
所以导出必须选一个映射（`--scipy-distance`）：

- `complement`（默认）`dist = C - sim`，`C` = 最大相似度。**这是语义上最忠实的一个**：
  本库的"缺失边 = 相似度 0"正好对应"距离 = C"，于是 `average` 链接等价于 scipy 的
  `method='average'` 作用在同一个稠密距离矩阵上 —— 实测在**没有并列**的用例上与
  scipy 逐行完全相同（`tests/test_scipy.py`，23 组）。
- `negate` `dist = -sim`，保序、无参，但距离为负。
- `reciprocal` `dist = 1/sim`，只对"相似度像亲和度"的数据有意义。

**三点要知道的差异**（不是 bug，是语义不同，`tests/test_scipy.py` 里都钉住了）：

1. **并列**。相似度相等时合并顺序本就不作保证。注意这包括**派生值**：输入 28 个值全不同，
   中间算出来的相似度照样会撞在一起（`15_many_merges.csv` 就是这样）。并列时两边树不同，
   但都合法。
2. **`MinNonzero` ≠ scipy 的 `complete`**。scipy 在缺失对上取到最大距离（= 相似度 0），
   稀疏输入下几乎处处为 0；本库只对"两边都存在的边"取最小。**稠密输入下两者相同**（实测 6/6）。
3. **森林**。`MinNonzero` 下两簇之间没有边时这个 pair 连候选都不是，分量可能永远合不起来，
   于是树不是完整树（行数 < n-1）。scipy 要求恰好 n-1 行，所以要么用 `--scipy-bridge`
   补桥接行（`complement` 下桥接距离就是 `C`，正是 scipy 自己会算出的值），要么换链接方式。
   不加 `--scipy-bridge` 时会在 stderr 上警告。

## 文档

| 文件 | 内容 |
|---|---|
| [docs/design.md](docs/design.md) | 数据结构、算法流程、可调参数 |
| [docs/known-issues.md](docs/known-issues.md) | 修过的 bug（含发现过程）+ 仍存在的问题 |
| [docs/performance.md](docs/performance.md) | 性能实测、被实测否掉的优化、方法学教训 |
| [docs/testing.md](docs/testing.md) | 切树属性 / 差分 / CLI / scipy 互操作 / C API（含 FFI）/ 性能护栏 |
| [CHANGELOG.md](CHANGELOG.md) | 对使用者有意义的变化 |

## 测试与性能的量级

- 差分测试：独立参考实现（不含 Lance-Williams 递推）逐次合并比对 + 22 个手写边界用例；
  另有 900 次多阈值切树交叉验证。当前 `PASS=90 WARN=12 KNOWN=0 FAIL=0`
  （WARN 是并列导致的合并顺序不确定，两边都合法）。
- 切树属性测试（`tests/test_cut.cpp`，库级）：期望值不是把并查集再写一遍，而是
  **从合并树自底向上收集叶子集合**算出——另一条代码路径。覆盖三种链接 × 21 个随机图 ×
  并列阈值 / 正好 k 个簇 / 最小簇大小 / 紧致度归属 / 单调性，外加多实例与 reset 复用。
  当前 `PASS=18281 FAIL=0`。
- scipy 互操作（`tests/test_scipy.py`）：把 `--scipy` 的导出交给**真的 scipy** 读。
  23 个用例 × 2 种链接，格式全部合规，且在**没有并列**的用例上与 scipy 自己的
  `linkage` **逐行完全相同**；并列用例按契约跳过并标注。
- ASan/UBSan 干净、0 泄漏；节点池手工毒化，所以"读已释放节点"仍能被 ASan 抓到。
- 累计性能（对比最开始的基线，详见 [docs/performance.md](docs/performance.md)）：
  800k 边 **21.9x**、2M 边 **51.1x**，峰值内存 **−55~61%**；单链接 **300x**。

## 许可

MIT，见 [LICENSE](LICENSE)。

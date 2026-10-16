# 设计与算法

本文写内部结构、算法流程和几个语义约定。使用接口见 [README](../README.md)，
历史与已知问题见 [known-issues](known-issues.md)，性能数据见 [performance](performance.md)。

## 数据结构

输入是**边表**而不是矩阵：`(obj1, obj2, sim)`，缺失的 pair 视为相似度 0，不存储。

```
vec_dim1 (行向量)                 vec_dim2 (列向量)
 row i: { head -> [单元链] , max -> 本行最大值的单元 }
 col j: { head -> [单元链] }
```

- `sp_unit` 的 `left/right/up/down` 四个链接是**下标**（不是指针），行链和列链共用同一批节点；
  全部节点放在一块连续内存里，`sizeof(sp_unit)` 32 字节。
- **只存上三角**：`add_grp()` 强制 `group1 < group2`（小的当行、大的当列）。
- 每行维护一个 `max` 指针 —— 这是关键优化：`get_max()` 只扫"每行的 max"这一列指针，
  复杂度是 O(组数) 而不是 O(nnz)。

## 算法流程

```
build(floor)                        src/hac.cpp
  add_edge(a, b, sim)               每出现一个新 id 就登记一个单元素组；边插进结构
  while (组数 > 1):
      unit = get_max()              全局最大相似度单元（堆顶，摊还 O(1)）
      if (无单元 || sim < floor) break
      merge_grp(group1, group2, sim)
```

`merge_grp()` 四步：

1. 遍历 group1/group2 的行链和列链，把需要重算的**行/列下标收进一个复用的小数组**
   （不是扫整个索引空间）、并标记 `max` 指针失效的行 —— 稀疏性的价值所在。
2. 按 linkage 方式算出新组与各邻居的相似度，新组 id 追加在末尾（`grp.size()`）；
   顺带把这一对相似度也记进 `marks`，所以下一步不需要再做链查找。
3. 摘除并释放 group1/group2 的整行整列，做四向指针 splice，最后释放中心节点。
4. 重扫失效的行修复 `max` 指针；把这一次合并记进 `merges`（就是 linkage 矩阵）。

单链接（`Linkage::MaxSimilarity`）走的是另一条路：簇间相似度就是两个簇之间所有**原始边**的
最大值，所以整棵树等价于"边按相似度降序 + 并查集"（Kruskal），完全不建矩阵。
见 [performance](performance.md)。

## 链接方式

三种，构造 `HAC` 时给。**命名按相似度语义**，而且注意原来代码里的注释是反的
（把 max 注释成 complete、把 min 注释成 single），这里以行为为准：

| `Linkage` | 簇间相似度 | 按距离语义是 | 备注 |
|---|---|---|---|
| `Average` | 全部 `|A|×|B|` 对的平均（缺失边算 0） | average | 簇越大越容易被"稀释" |
| `MaxSimilarity` | 取最大 | **单链接** | 会产生 chaining；值只可能是原始边，所以走 Kruskal |
| `MinNonzero` | 取最小，但仅当两边都有边 | 全链接的保守变体 | 缺失边不参与，否则 min 恒为 0 |

`MinNonzero` **不是**标准 complete linkage：稀疏图上把缺失边当 0 会让 min 恒为 0、
失去意义，所以这里只对"两边都存在的边"取 min，名字也就显式写成 `MinNonzero`。
它的一个后果是**可能留下永远合不起来的分量**（两簇之间没有边 ⇒ 这对不是候选），
于是树可能是残缺的 —— `cut_n_clusters` / `cut_min_size` 在这种情况下返回树能给出的
最细切法。另外它的"组内紧致度"初值是空集的 +∞（`empty_cohesion`），
因为 0 是 min 的吸收元、会让报告恒为 0。

## 状态与实例

原来这里是 **34 个文件静态变量 + 一堆静态函数** —— 一个进程只能跑一次聚类，
换链接方式还得重新编译。现在这些状态都在 `spclust::detail::HACImpl` 里，
公开的头文件只有一个 `class HAC`（内部持一个 `unique_ptr`）：

- **一个对象 = 一份状态**：同一进程可以有任意多个聚类，交错喂边、交错 `build()`，
  结果与单独跑完全一致（`tests/test_cut.cpp` 里专门测了这条，包括 reset 复用）。
- **库里没有可变全局量**：`nm libspclust.a` 的 `.data`/`.bss` 里只有 C++ 运行时的
  常量（字符串、跳转表），没有算法状态。所以只要不共享同一个对象，多线程就是安全的。
- **链接方式运行时选**：`HAC(Linkage)`。分派只发生在公开入口（`add_edge` / `build` /
  `cut` / …），内层循环里链接方式是已知的同一个函数内联体，热路径上不留分支 ——
  实测代价 0.2~0.6%（交替 A/B，最小/四分位/中位一致）。
- **切树是纯函数**：并查集的父数组是 `cut()` 的局部量。以前它是文件静态量，
  于是 `cut()` 名义上 const 却在改全局状态，并发切树会互相踩。

## 切树

树是完整保留的（`build(floor)` 的 `floor` 只决定"算到多低"），所以阈值只是"切在哪"，
换切法不需要重新聚类。三个接口都归约到"取前 use 次合并"这一步：

| 接口 | 语义 | 需要的树 |
|---|---|---|
| `cut(t)` | 取所有 `similarity >= t` 的合并（并列会一起进来） | 算到 `t` 即可 |
| `cut_n_clusters(k)` | 正好 k 个簇 = 取前 `n-k` 次合并 | 要够深，拿不到 k 个就给最细的 |
| `cut_min_size(m)` | "所有簇都 >= m"的**最浅**切法 | 要够深 |

`cut_min_size` 不是"丢掉小簇"，而是"继续合并直到没有小簇"。注意有三种链接语义下
**两簇之间可以根本没有边**（`MinNonzero` 下这种 pair 连候选都不是），于是可能留下
互不相连的分量永远合不起来 —— 这时它返回树能给出的最细切法。

## 可调参数

**库本体没有编译期选项**（链接方式是构造参数）。剩下的都是命令行工具自己的，
在 `tools/cli_config.h` 里，用 `#ifndef` 包住所以命令行 `-D` 优先：

```bash
make MIN_SIM=60 DEBUG=2          # 默认切树阈值、调试级别
```

- `MIN_SIM_CLUST`：`--cut` 的默认值（默认 50）
- `MAX_SIM_READ`：输入相似度上界（默认 65535）。**这是历史遗留**（原来是
  `unsigned short`），库本身接受任意有限非负 `double`
- `DEBUG`：0=静默 1=阶段与进度 2=每次合并的细节
- `DEBUG_FREQ_READDB`：读入进度的打印频度

链接方式、切法、`--floor` 都是运行时参数，见 `spclust_cli --help`。

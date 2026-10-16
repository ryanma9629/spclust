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

单链接（`LINKAGE_MTD == 1`）走的是另一条路：簇间相似度就是两个簇之间所有**原始边**的最大值，
所以整棵树等价于"边按相似度降序 + 并查集"（Kruskal），完全不建矩阵。见 [performance](performance.md)。

## 可调参数

全部在 `src/hac.cpp` 里，都用 `#ifndef` 包住，所以可以在命令行覆盖、不用改源码：

```bash
make run LINKAGE=1 MIN_SIM=60 DEBUG=2        # 改默认切树阈值
make run MIN_SIM=0 LINKAGE=0                  # hmm，MIN_SIM=0 只改默认；树深度见 --floor
```

运行时参数（不改源码、不重编译）：

```bash
./spclust data.csv --cut 80 --tree            # 切在 80 并打印树
./spclust data.csv --floor 0 --matrix         # 完整树 + linkage 矩阵
```

| 变量 | 默认 | 含义 |
|---|---|---|
| `LINKAGE` | 0 | 0=average；1=max（原注释写 complete）；2=min，仅当两条边都存在时才更新（原注释写 single） |
| `MIN_SIM_CLUST` | 50 | **默认的切树阈值**（不再影响聚类本身）。可用 `--cut` 覆盖 |
| `DEBUG` | 1 | 0=静默，1=阶段/进度，2=每次合并的细节 |
| `DEBUG_FREQ_READDB` | 10000 | 读数据的进度打印频度 |
| `DEBUG_FREQ_CLUST` | 100 | 聚类进度打印频度（建议 >=100） |

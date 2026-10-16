# Changelog

格式参考 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)，版本号遵循
[语义化版本](https://semver.org/lang/zh-CN/)。

本文件记录对使用者有意义的变化。逐条的实现细节（含每次性能改动的实测数据与方法学修正）
在 git 历史里，每个提交都是独立可回滚的。

## [0.2.0] - 2026-10

库与使用者的分离（自由函数 API）。

### Added
- **库与使用者分离**：`include/spclust/hac.h` + `src/hac.cpp` 是库本体（不依赖任何东西、
  只认整数 id、不打印、不读文件）；`tools/` 是它的第一个使用者，负责所有 IO 与展示。
- **CMake**：`find_package(spclust)` + `spclust::spclust` + `install` + `ctest`；
  `SPCLUST_LINKAGE` / `SPCLUST_BUILD_TOOLS` / `SPCLUST_BUILD_TESTS` 选项。
  `Makefile` 保留作开发时的快捷方式。
- 库公开接口：`reserve` / `add_edge` / `build` / `num_objects` / `tree` / `cut` / `reset` /
  `linkage_name`；切树（`Cut`）现在也在库里（原来在 CLI 的匿名命名空间里，那属于核心功能）。
- `docs/`：设计、已修复与已知问题、性能实测、测试四份文档；README 只讲"怎么用"。
- GitHub Actions CI：两个平台 × (构建 + ctest + ASan/UBSan 全用例 + 只编库)。
- 完整的合并树输出：`--cut <X>`（把树切在相似度 X 处）、`--floor <X>`（树算到多低）、
  `--tree` / `--tree-depth <N>`（树状图）、`--matrix`（linkage 矩阵 `left right height size`）。
  阈值不再是停止条件，而是切树位置 —— 换阈值不需要重新聚类。
- 性能回归护栏：`make bench`（固定输入 + 记录基线 + 机器状态哨兵）、`make bench-update`。
- CLI 接口测试：`tests/test_cli.py`，覆盖退出码契约、切树与 oracle 的一致性、
  `--matrix` 结构不变量、`--tree` 完整性、stdin 等价性。
- 差分测试套件：`tests/oracle.py`（独立的按定义实现的参考实现，不含 Lance-Williams 递推）
  + `tests/run_tests.py`（逐次合并比对）+ 22 个手写边界用例。
- `LICENSE`（MIT）与每个源文件的 SPDX 头。

### Changed
- 库不再持有对象名：`add_edge(a, b, sim)` 收整数 id，`id ↔ 名字` 的映射归 CLI。
  库也不再有任何标准输出（原来每次合并的调试日志会打到 stdout）。
- 库内部不再 `exit()`：内存分配失败抛 `std::bad_alloc`，参数非法抛 `std::invalid_argument`。
- 相似度从 `unsigned short` 改为 `double`。原实现每步整数除法截断、误差向下累积，
  实测会把不同的相似度压成同一个值而产生并列，**改变最终簇数**。
- 单链接（`LINKAGE_MTD == 1`）改走 Kruskal：单链接的簇间相似度就是两个簇之间所有
  **原始边**的最大值，所以整棵树等价于"边按相似度降序 + 并查集"，不需要相似度矩阵。
  实测 63~300x（随规模），峰值内存从 15.8 GB 降到 50 MB（n=80000 / 800k 边）。
- 重复 pair 的语义定清为**取最大值**（此前是未定义行为）。
- 读取端改为**小批量流式**：不再把全部边缓存进内存，峰值内存下降 40~49%。
- 节点从 48 字节压到 32 字节（链接与组号改存 `uint32` 下标），峰值内存再降 ~18%，
  同时快 1.4~1.5x。
- `readdb` 不再内部 `exit(1)`：失败时返回 0（原因写 stderr），退出码交给调用方。
- 头文件不再 `using namespace std;`；`grp` 不再是跨文件的全局变量。
- `LINKAGE_MTD == 2` 的紧致度报告此前恒为 0（`min(...)` 的吸收元），现在报告簇内最小边。

### Fixed
- `get_max()` 的 heap-use-after-free：某行的相似度全为 0 时，max 指针不会被重新赋值，
  会保留一个已释放的节点。默认的 average 链接 + 一条 `sim=0` 的输入边即可触发。
  （差分测试对这个输入报 PASS —— 只有 ASan 抓到了它。）
- 未使用的 `MIN_ACCT_NUM` / `MIN_SIM_WG` 宏已删除。

### Performance
累计（对比最初的基线）：800k 边合计 **21.9x**、2M 边 **51.1x**；峰值内存 **-55~61%**。
单链接 **300x**。主要来源：去掉每次合并的 4 趟 O(n) 全量扫描、`get_max` 改最大堆、
`init` 改哈希表、去掉 `get_unit` 的线性查找、流式读入、32 字节节点、单链接走 Kruskal。

## [0.3.0] - 2026-10

这一版**改动了公开 API**（自由函数 → `spclust::HAC` 对象），因为之前的接口有一个
根本限制：状态是 34 个文件静态变量，**一个进程只能跑一次聚类**，而且换链接方式
要重新编译整个库。

### Changed

- **`class spclust::HAC` 取代原来的自由函数**（`reserve` / `add_edge` / `build` /
  `num_objects` / `tree` / `cut` / `reset` / `linkage_name` 都成了它的成员）。
  一个对象一份状态：同一进程可有任意多个聚类，交错使用互不影响；库里不再有可变全局量
  （`nm libspclust.a` 里 `.data`/`.bss` 无算法状态），所以只要不共享同一个对象，
  多线程安全。
- **链接方式从编译期 `#define LINKAGE_MTD` 变成构造参数 `spclust::Linkage`**：
  `Average` / `MaxSimilarity` / `MinNonzero`。以前要三个二进制才能跑三种链接，
  现在一个二进制就够（`spclust_cli --linkage 0|1|2`）。分派只在公开入口发生，
  热路径上不留分支 —— 实测代价 0.2~0.6%（交替 A/B，最小/四分位/中位一致）。
- **链接方式的名字按语义定稿**。原来代码里的注释是反的（`max` 注释成 complete、
  `min` 注释成 single），现在以行为为准：`max` = 单链接（走 Kruskal）、
  `min` = "仅对两边都存在的边取最小"，所以显式命名 `MinNonzero` 而不是 complete。
- `cut()` 变成纯函数：并查集的父数组改成局部量（原来是文件静态量，于是 `cut()`
  名义上 `const` 却在改全局状态，并发切树会互相踩）。
- 命令行工具：新增 `--linkage`、`--n-clusters <K>`、`--min-size <M>`；
  `--n-clusters` / `--min-size` 需要完整树，没给 `--floor` 时会自动用 `--floor 0`。
- CMake 的 `SPCLUST_LINKAGE` 选项删除（不再需要）；`SPCLUST_BUILD_TOOLS` /
  `SPCLUST_BUILD_TESTS` 现在用 `option()`，所以 `-D...=OFF` 真的生效
  （之前用 `set()`，会把用户传的 OFF 覆盖掉）。

### Added

- `Cut cut_n_clusters(std::size_t k)`：要正好 k 个簇（= 树的前 n-k 次合并）。
- `Cut cut_min_size(std::size_t m)`：合并到没有小于 m 的簇为止（满足条件的**最浅**切法，
  不是"丢掉小簇"）。
- `tests/test_cut.cpp`：库级切树属性测试。期望值用**从合并树自底向上收集叶子集合**
  独立算出（不是把并查集再写一遍），覆盖三种链接 × 21 个随机图 × 并列阈值 /
  k 个簇 / 最小簇大小 / 紧致度归属 / 单调性，外加**多实例与 reset 复用**。
  当前 `PASS=18281 FAIL=0`。
- `make test-cut`。

### Fixed

- `docs/design.md` 里的链接方式命名与流程图（原来写的是 `LINKAGE_MTD` 和旧函数名）。

## [0.1.0] - 2026-10

最初可用状态：把 2022 年的 Windows/ODBC 版本改成可移植（CSV / 标准输入）、建立差分测试套件。

# spclust

基于**稀疏相似度矩阵**的凝聚型层次聚类（agglomerative hierarchical clustering）。
反复取出全局相似度最大的一对簇并合并，直到最大相似度低于阈值或只剩一个簇。

核心实现是 `spclust.cpp`：用**十字链表**（orthogonal linked list）存稀疏上三角矩阵，
使每次合并只需重算"与这两个簇相邻的那几行几列"，而不是整个矩阵。

## 数据结构

输入是**边表**而不是矩阵：`(obj1, obj2, sim)`，缺失的 pair 视为相似度 0，不存储。

```
vec_dim1 (行向量)                 vec_dim2 (列向量)
 row i: { head -> [单元链] , max -> 本行最大值的单元 }
 col j: { head -> [单元链] }
```

- `sp_unit` 有 `left/right/up/down` 四个指针；行链和列链共用同一批节点。
- **只存上三角**：`add_grp()` 强制 `group1 < group2`（小的当行、大的当列）。
- 每行维护一个 `max` 指针 —— 这是关键优化：`get_max()` 只扫"每行的 max"这一列指针，
  复杂度是 O(组数) 而不是 O(nnz)。

## 算法流程

```
clust_grp()                       spclust.cpp
  init_grp()                      每个对象一个单元素组，把全部边插进结构
  while (组数 > 1):
      unit = get_max()            全局最大相似度单元
      if (无单元 || sim < MIN_SIM_CLUST) break
      merge_grp(group1, group2, sim)
```

`merge_grp()` 四步：

1. 遍历 group1/group2 的行链和列链，标记**受影响的边**（`update_dim1` / `update_dim2`）
   以及 `max` 指针失效的行（`update_max`）——稀疏性的价值所在。
2. 按 linkage 方式算出新组与各邻居的相似度，新组 id 追加在末尾（`grp.size()`）。
3. 摘除并释放 group1/group2 的整行整列，做四向指针 splice，最后释放中心节点。
4. 重扫受影响的行，修复 `max` 指针；把新组的成员列表与组内平均相似度写进 `grp`。

## 构建与运行（macOS）

不需要安装任何 C++ 运行时 —— Apple 的 libc++ 是操作系统自带的。
只需要 Xcode Command Line Tools（`xcode-select --install`），它提供 `clang++`、`make`、`lldb`。

```bash
make                # 构建 ./spclust
make run            # 构建并跑 sim.csv
make test           # 差分测试：独立 oracle 交叉验证正确性
make check          # 用三种 linkage 各跑一遍做对比
make asan           # AddressSanitizer + UBSan（强烈建议跑一次）
make asan-cases     # 用 ASan/UBSan 跑遍 tests/cases 全部用例
make leaks          # macOS 自带 leaks 工具查内存泄漏
make clean
make help
```

只用了标准库，没有第三方依赖，`g++`/`c++` 都能编，`-std=c++11` 也可以。

## 输入格式

`obj1,obj2,sim`，一行一条。也可以从标准输入读（`./spclust < data.csv` 或 `./spclust -`）。

- `sim`：0..65535 的实数（**允许小数**，如 87.5）；低于 `MIN_SIM_CLUST` 的边只作为候选，不会被合并
- 分隔符：含逗号按逗号拆，否则按空白（空格/Tab）拆
- `#` 开头是注释；空行忽略
- 可选的表头行 `obj1,obj2,sim` 会被自动跳过
- 越界/非数字的 `sim`、字段数不足、自环（`obj1 == obj2`）都会打印警告并跳过该行
  （自环会破坏"group1 < group2"的前提，必须丢弃）
- 重复的 pair 会被重复插入，但算法只会用到其中一个（结果不受影响，只是多占点内存）

## 输出

```
Linkage=0=average, min_sim=50
Objects: 16, groups: 3
Group 0 (size=5, sim=90.5): a4 a5 a3 a1 a2
Group 1 (size=4, sim=87.6667): b4 b3 b1 b2
Group 2 (size=3, sim=82): c3 c1 c2
Unclustered (4): d1 d2 e1 e2
```

`sim` 是**组内平均相似度**（average linkage 下按组大小加权累加得到），可以当作簇的紧致度。
只列出 size>=2 的簇，剩下的单元素对象列在 `Unclustered` 里。
注：修复前这里是截断后的整数（`sim=90` / `sim=87`），现在才是真实的平均值。

## 测试（差分验证）

`spclust.cpp` 是手写 `malloc/free` + 四向指针 + 十字链表的代码，这类代码的正确性
没法靠读代码确认。所以 `make test` 做**三方交叉验证**：

```
1. oracle(brute)  tests/oracle.py  每步从原始稠密矩阵按定义重算。最可信，O(n^4)。
2. oracle(fast)   tests/oracle.py  用恒等式增量维护聚合量。O(n^3)。
3. spclust                          被测对象。
```

`oracle.py` 里**没有** Lance-Williams 递推、没有链表，只有“按定义算 + 取全局 argmax”，
所以它是独立的。三种实现都用精确整数运算（`Fraction`），不引入浮点误差。

判定规则：spclust 与 oracle 的**逐次合并**（成员集合 + 相似度）以及最终分组必须一致。
不一致时分四类：

| 标记 | 含义 |
|---|---|
| `PASS` | 与精确定义完全一致 |
| `WARN` | 并列最大（或截断造成的并列）导致合并顺序不确定，两边都可能合法 |
| `TRUNC` | 偏差完全由 `unsigned short` 整数截断解释（见“已知问题 1”）|
| `KNOWN` | 记录在 `tests/known_divergences.txt` 里、已确认可由实现语义解释的分歧 |
| `DRIFT` / `FAIL` / `ORACLE` | 未解释的偏差 / 结果不同 / 测试工具自身不一致 —— 这三种都会让 `make test` 失败 |

oracle 自身的可信度靠 brute 与 fast 互相校验（小规模时两者必须逐次一致，报告里会打印次数）。

用例 = `tests/cases/*.csv`（22 个手写边界用例，每个文件的注释写明它在探测什么）
+ 12 个随机稀疏图用例（相似度互不相同，避免并列）。全套在 1 秒左右跑完。

```bash
make test                      # 默认 min_sim=50
make test MIN_SIM=70           # 换阈值重跑（oracle 会同步，被测二进制也会重编）
python3 tests/run_tests.py --only average --verbose
python3 tests/run_tests.py --quick --case 15_   # 只跑指定用例
python3 tests/oracle.py sim.csv --linkage average --min-sim 50   # 单独用 oracle 看结果
```

当前结果（`min_sim=50`，22 个手写 + 12 个随机 × 3 种 linkage = 102 次比对）：

```
PASS=88  WARN=8  KNOWN=6  TRUNC=0  DRIFT=0  FAIL=0  ORACLE=0
oracle 自检（brute 与 fast 逐次一致）: 96 次
RESULT: PASS
```

结论：**`spclust.cpp` 的链表维护本身是可靠的** —— 在 44 个对象、43 次连续合并、
七种输入形态（含稠密完全图、星型、长链、多连通分量、零相似度边、逆序输入、小数相似度）上，
`max`/`min` 两种链接的 64 次比对里只有 WARN（并列）和已登记的重复-pair 分歧。

另外 `average` 的相似度经修复后与精确定义的残差是纯浮点量级：

```
排除两个重复-pair 语义用例，146 次逐次相似度比对：
  最大绝对误差 = 3.3e-11
  最大相对误差 = 4.3e-13   （DEBUG 日志保留 12 位有效数字，这基本就是打印精度）
```

`TRUNC` 现在是**回归护栏**：如果谁把相似度改回 `unsigned short`，那 12 个 `TRUNC` 会
全部重新出现，而 `TRUNC` 算**失败**。

## 可调参数

全部在 `main.h` 里，都用 `#ifndef` 包住，所以可以在命令行覆盖、不用改源码：

```bash
make run LINKAGE=1 MIN_SIM=60 DEBUG=2
```

| 变量 | 默认 | 含义 |
|---|---|---|
| `LINKAGE` | 0 | 0=average；1=max（原注释写 complete）；2=min，仅当两条边都存在时才更新（原注释写 single） |
| `MIN_SIM_CLUST` | 50 | 相似度下限，低于此值停止合并 |
| `DEBUG` | 1 | 0=静默，1=阶段/进度，2=每次合并的细节 |
| `DEBUG_FREQ_READDB` | 10000 | 读数据的进度打印频度 |
| `DEBUG_FREQ_CLUST` | 100 | 聚类进度打印频度（建议 >=100） |

## 与原 Windows/ODBC 版本的差异

原版用 `<windows.h>` + `_tmain` + ODBC 从 MySQL 的 `test.sim100` 读数据，只能在 Windows/MSVC 下编译。
移植阶段 `spclust.cpp` **一行未改**；之后三次改动动了它：相似度精度、`get_max` 的 UAF、
以及把每次合并的 4 趟 O(n) 全量扫描改成 O(degree)（见「已修复」和「性能实测」）：

| 文件 | 改动 |
|---|---|
| `main.h` | Windows/ODBC 头文件包进 `#ifdef _WIN32`；补齐标准库头文件；可调参数加 `#ifndef` 以便命令行覆盖 |
| `main.cpp` | `_tmain(argc, _TCHAR* argv[])` → `main(argc, char* argv[])`；支持可选输入文件参数；打印结果（原来那段被注释掉了） |
| `readdb.h` | 签名由 `readdb(void)` 改为 `readdb(const char* path)`；删掉 ODBC 专用的 `odbc_show_err` |
| `readdb.cpp` | 整个 ODBC 读取换成 CSV / 标准输入读取；不再有 `char obj1[32]` 的截断限制 |
| `spclust.cpp` / `spclust.h` | 精度修复（`unsigned short` → `double`）、`get_max` 的 UAF、增量下标遍历替代全量扫描 |
| `Makefile` / `sim.csv` / `.gitignore` | 新增 |
| `tests/` | 新增：独立 oracle + 22 个边界用例 + 差分测试 runner |

原来的 `readdb.cpp` 还在 git 历史里：`git show 9c1afc7:readdb.cpp`。

## 已修复

### `average` 的相似度曾用 `unsigned short` 存，每步截断一次，误差会累积

修复前的 `spclust.cpp:196`：

```cpp
unsigned short sim = (sim1 * group1_size + sim2 * group2_size) / (group1_size + group2_size);
```

整数除法截断，结果又存回 `unsigned short`，于是**每次合并都向下丢一点**，并且丢掉的量会
参与下一步计算继续传播。`make test` 当时报了 12 个 `TRUNC`：

```
rand_n44d80s1.csv（44 对象）第 27 步: oracle=41869.6  spclust=41868   （15 次合并里 12 步偏低）
rand_n10d30s3.csv: 截断把 262.5 和 262.0 压成同一个值 262 → 产生并列
                   → 合并顺序变化 → oracle 得到 3 个簇，spclust 得到 4 个
```

关键点：这不只是显示问题 —— **判据本身**（取全局最大相似度）就在被截断后的值上做。

**修法**：`sp_sim::sim` / `sp_unit::sim` / `grp_s::sim` 与相关局部变量、
`add_grp`/`merge_grp` 的形参全部改成 `double`（`main.h`、`spclust.h`、`spclust.cpp`）；
`readdb.cpp` 的 `parse_sim` 改用 `strtod`，顺带开始接受小数（如 `87.5`）。
`max`/`min` 两种链接本来只比较不平均，不受影响。

**副作用（正向）**：以前 `if (sim > 0)` 会因为截断把真实存在但很小的平均值变成 0、
从而丢掉这个候选；现在不会了，所以候选集更贴合定义。

**验证**：修复后 `TRUNC=0`，12 个用例全部转为 `PASS`；
146 次逐次相似度比对的最大相对误差 `4.3e-13`（纯打印精度）。
并且把 `TRUNC` 改成了**失败**类别，所以这个 bug 再回来会直接让 `make test` 红。

### `get_max()` 里的 heap-use-after-free（行内相似度全为 0 时）

max 修复循环用的是 `if (u->sim > max)`，`max` 初值 0，于是**整行相似度都是 0 的行永远不会被重新赋值**。
如果这一行原来的 max 指向的单元正好在本次被删的列里，`vec_dim1[i].max` 就会保留一个已 free 的指针，
下一次 `get_max()` 直接读它。

触发条件（`tests/cases/22_zero_row_max.csv`）：某行的 max 是被删列里的一个 0 相似度单元，
同时该行还有另一个 0 单元**存活**（所以这一行没有变空 —— 变空的情况原本已经处理）。
默认的 average 链接 + 一条 `sim=0` 的输入边就能触发，不需要 `LINKAGE_MTD=1`。

修法：扫出该行的真实最大值并无条件赋值，指针永远要么是有效单元要么是 NULL。
并列行为不变（仍是"链上第一个取得最大值的单元"），而 0 值的 max 永远不会被 `get_max()` 选中
（那里是严格 `>`）。

**值得注意的是它是怎么被发现的**：差分测试对这个输入报 **PASS** —— 因为刚 free 的内存读回来恰好还是 0，
合并顺序也一样。**只有 ASan 抓到了它**。所以结果比对和内存安全是两个独立的检查，缺一不可。

### 每次合并的三趟 O(n) 全量扫描（性能）

旧 `merge_grp` 为每次合并分配 3 个 O(n) 的标记向量，然后

```cpp
for (size_t i = 0; i < dim1_size; i++)   // dim1_size ≈ 2n
    if (update_dim1[i]) { ... }          // 被标记的只有 O(degree) 个
```

即**用一个 O(n) 的全量扫描去筛出 O(degree) 个下标**。加上 `get_max()` 和 max 修复循环，
每次合并共 4 趟全量扫描。实测（n=40000、18293 次合并）：

```
update_dim1 循环   926,710,622 次迭代
update_dim2 循环   927,123,378 次迭代
max 修复循环       926,7xx,xxx 次迭代
get_max()          926,769,374 次迭代
------------------------------------------------
合计              ≈3.7e9 次迭代，等效约 0.8 ns/次
```

修法：只遍历被标记的下标，标记缓冲跨合并复用（只复位本次触及过的位置），于是每步是 O(degree) 而不是 O(n)。
顺带解决两件事：

- **两个更新循环的循环体逐字相同**，对同一个下标 x 都是 `get_unit(x, group1)` / `get_unit(x, group2)`。
  所以 x 同时落在两趟里时，旧实现会为同一个 `(x, grp_size)` 插入两个数值相同的**重复单元**。
  合并成一趟后这些重复单元消失（n=40000 时 malloc 次数从 1,121,925 降到 1,097,456，少 2.2%）。
- 每次合并不再分配 O(n) 内存（中间那 3 个标记向量）。

**实测 2.2 倍加速，且输出逐字节相同**（28 个输入 × 3 种 linkage，包含 DEBUG=2 的完整合并日志）。
但**复杂度类别没有变**，因为 `get_max()` 仍然是每次合并 O(n)，即 O(n²) —— 见下面的性能实测。

## 性能实测

环境：Apple M 系列，`-O2`，数据为随机稀疏图（平均度 20，相似度 0..100，`MIN_SIM_CLUST=50`）。

### 内存：线性，与稀疏性一致

| 边数 | n | 峰值 RSS | 字节/边 |
|---|---|---|---|
| 100k | 10k | 21.7 MB | 227 |
| 200k | 20k | 41.4 MB | 217 |
| 400k | 40k | 85.0 MB | 223 |
| 800k | 80k | 163.9 MB | 215 |

字节/边恒定 → 十字链表确实只存 nnz，这个设计目标达成了。

### 时间：修之前是 O(n²)，修之后仍是 O(n²)（但常数小了 2.2 倍）

聚类段耗时（`/usr/bin/time` + 代码内分段计时）：

| 边数 | 合并次数 | P0-1 修复前 | P0-1 修复后 | 加速 |
|---|---|---|---|---|
| 100k | 4,566 | 0.096 s | 0.048 s | 2.0x |
| 200k | 9,155 | 0.346 s | 0.154 s | 2.2x |
| 400k | 18,293 | 1.436 s | 0.669 s | 2.1x |
| 800k | 36,564 | 5.538 s | 2.406 s | 2.3x |

### 时间到底花在哪（P0-1 之后的实测构成）

| 边数 | 聚类 | 其中 `get_max()` | 占比 | 其余（合并主体 + malloc）|
|---|---|---|---|---|
| 100k | 0.052 s | 0.029 s | 57% | 0.022 s |
| 200k | 0.158 s | 0.109 s | 69% | 0.050 s |
| 400k | 0.660 s | 0.456 s | 69% | 0.204 s |
| 800k | 2.789 s | 2.187 s | **78%** | 0.602 s |

所以**下一个瓶颈已经很明确**：`get_max()` 每次合并扫一遍所有行（O(组数)），总共 O(n²)。
把它换成最大堆 + 懒删除（`sparseAHC` 就是这么做的），预计还能再快约 4 倍，
并把复杂度从 O(n²) 降到近似 O((n + nnz) log n)。

### 其它已排除的假设（都是实测排除的，不是推断）

| 怀疑对象 | 实测 |
|---|---|
| `grp[].objs` 的字符串复制 | 对象名 8 字符 vs 64 字符，耗时 2.915 vs 2.915 s → **不是问题** |
| 标记向量的 memset | `size_t` 改 `char`，2.91 → 2.80 s（4%）→ 不是主因 |
| `malloc`/`free` 开销 | n=40000 共 1.12M 次（2.8× 边数），推算 <10% |
| `get_max` 是最大头 | 单独计时：修前占 31%，修后占 57~78% |
| `LINKAGE_MTD=1` | 2.9 s → **24.0 s**（8.3x）→ 无条件 `add_grp` 导致稠密化是真的 |

## 已知问题

以下都是实测/静态分析得出的结论。
前两项**已由 `make test` 复现**（有对应标记），后面几项只是静态分析。

### 1. 重复 pair 是未定义输入，且 `get_max` 与 `get_unit` 会指向不同单元（`KNOWN`）

同一个 `(obj1, obj2)` 出现多次且 `sim` 不同时，读取端会把**每一条都插入结构**，于是：

- `get_max()` —— 取行内最大值 → 看到较大的那个
- `get_unit()` —— 沿行链从头找第一个匹配（`add_grp` 头插）→ 拿到最后插入的那个

两者指向**不同单元**。合并不崩溃（两个单元都会被恰好释放一次，`make asan-cases` 通过），
而且实测效果等价于“**重复 pair 取较大值，与写入顺序无关**”。但这是隐含行为，应该由读取端
明确策略（去重 / 报错 / 取最大）而不是留给聚类内部。
用例：`tests/cases/07_duplicate_edge.csv`、`19_duplicate_order.csv`。

### 2. 自环（`obj1 == obj2`）会破坏结构前提，读取端已丢弃

`add_grp` 的所有逻辑都建立在 `group1 < group2` 上，自环会让 `group1 == group2`，
行链和列链互相踩踏。当前的 `readdb.cpp` 会打印警告并跳过，这是**必要**的，不是可选项。

### 3. `max`/`min`（`LINKAGE_MTD` = 1/2）不是标准链接，且 `1` 会破坏稀疏性

`1`（注释写 complete）取 `max`，按距离语义的标准术语是**单链接**；
`2`（注释写 single）取 `min` 且要求两边都存在，是**全链接**的一个变体。
注意：这个命名与 `sparseAHC`（稀疏相似度矩阵 AHC 的参考实现）的文档完全一致，
在“相似度矩阵”这一支里是惯例，不一定是笔误。

真正的问题在 `1`：它的 `add_grp(i, grp_size, max(sim1, sim2))` 是**无条件插入**的
（哪怕值是 0），反复合并后矩阵会趋近稠密，内存和时间退化到 O(n²)。
`0` 和 `2` 都有 `if (sim > 0)` 保护，`1` 没有。

### 4. `LINKAGE_MTD=2` 时 `grp[].sim` 恒为 0

单元素组的 `sim` 初值是 0，而 `min(min(...), max_sim)` 一旦碰到 0 就永远是 0。
只影响输出的紧致度报告，不影响合并顺序。

### 5. “缺失边 = 0”这个语义本身会稀释簇间相似度

不是 bug，但值得写进文档：平均链接的分母是全部 `|A|·|B|` 对，所以簇越大，越容易被“稀释”。
用例 `05_chain_bridge.csv` 就是实例：`a3-b1 = 58` 是两组之间唯一的高于阈值的边，
但 `a1-a3` 缺边，一次合并后 `{a1,a2}` 对 `a3` 的相似度就只有 `(0+94)/2 = 47`，
桥在触发合并之前就沉到阈值以下了。`sparseAHC` 采用同样的分母，所以这不是实现差异，
而是这种建模方式的固有性质；如果你的簇形状应该由连通性（而不是平均数）决定，
要在归一化上做选择（比如分母改成实际存在的边数）。

### 6. 以下只是静态分析，测试未区分（因为结果不受影响）

- `grp` / `vec_dim1` / `vec_dim2` 单调增长到约 2n-1，被合并的槽位清空但不回收。
- `MIN_ACCT_NUM`、`MIN_SIM_WG` 定义了但从未使用。
- `using namespace std;` 放在头文件里，`grp` 靠全局变量 + `extern` 串起来。

### 验证记录

- `make asan` / `make asan-cases`：在 `sim.csv` 和全部 22 个用例上无内存错误、无 UB。
- `make leaks`：`0 leaks for 0 total leaked bytes`（`free_mem()` 释放了全部链表节点）。
- `make test`：`min_sim` = 20 / 50 / 70 / 200 四个阈值下均 `RESULT: PASS`，且 `TRUNC=0`。
- 等价性：22 个用例 + 3 个随机图 + `sim.csv` 共 28 个输入 × 3 种 linkage，
  完整 DEBUG 输出（含每次合并的相似度）**逐字节相同**，仅打印的整秒数不同。
- 性能：P0-1 后聚类段加速 2.0~2.3x，峰值内存 85→82 MB。

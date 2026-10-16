#ifndef MAIN_H
#define MAIN_H

#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <unordered_map>
#include <string>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------
// 可调参数。全部用 #ifndef 包住，因此可以在命令行覆盖而不改源码：
//     make LINKAGE=1 MIN_SIM=60 DEBUG=2
// 或在命令行直接给编译器加 -DLINKAGE_MTD=1
// ---------------------------------------------------------------------------

//读取数据过程中调试信息显示频度
#ifndef DEBUG_FREQ_READDB
#define DEBUG_FREQ_READDB 10000
#endif

//聚类过程中链接方法
//0:average linkage
//1:complete linkage
//2:single linkage
#ifndef LINKAGE_MTD
#define LINKAGE_MTD 0
#endif

//相似度下限，低于此值将停止聚类
#ifndef MIN_SIM_CLUST
#define MIN_SIM_CLUST 50
#endif

//聚类过程中调试信息出现的频度，实际应用中建议至少设定为100
#ifndef DEBUG_FREQ_CLUST
#define DEBUG_FREQ_CLUST 100
#endif

//调试信息级别：0=静默 1=阶段与进度 2=每次合并的细节
#ifndef DEBUG
#define DEBUG 1
#endif

//这里**不放** `using namespace std;`：那会把整个 std 命名空间注入到所有包含本头文件的
//翻译单元里。头文件里的类型一律写 std::，各 .cpp 自己决定要不要 using。

//单元素组的"组内平均相似度"初值。
//
//average / max 的合并式里 0 是合适的初值；但 min 链接的合并式是 min(min(g1,g2), max_sim)，
//0 是吸收元，会让紧致度报告恒为 0。单元素组本来就没有"组内边"，
//它的初值应该是**空集的 +∞**，这样 min 才有意义。
#if 2 == LINKAGE_MTD
static const double EMPTY_COHESION = std::numeric_limits<double>::infinity();
#else
static const double EMPTY_COHESION = 0.0;
#endif

//聚类过程中保存相似性的节点。
//
//四个链接和两个组号都存**下标**而不是指针，整个池是一块连续内存：
//   sizeof(sp_unit) = 4*4(链接) + 2*4(组号) + 8(sim) = 32 字节（原来 48）
//下标 0 保留为空哨兵（见 SP_NONE）。
//实测节点尺寸是一阶因素（把节点人为撑到 128 字节会慢 1.17x），所以这一压是值得的。
//
//组号范围是 [0, 2*对象数)，链接下标范围是 [0, 单元数)，两者 uint32 都够；
//代价是对象数上限约 21 亿 —— 到那个量级地址空间会先耗尽。
//空链接用 0。0 号节点是**保留的空哨兵**，不参与任何链表 ——
//这样 vec_dim1/vec_dim2 resize 时的值初始化（清零）天然就是"空"，
//与原来 head/max 是裸指针时零初始化即 NULL 的行为一致。
static const uint32_t SP_NONE = 0;

struct sp_unit {
	uint32_t right;
	uint32_t down;
	uint32_t left;
	uint32_t up;
	uint32_t group1;
	uint32_t group2;
	double sim;
};
//“左向量”中的元素
struct p_dim1 {
	uint32_t head;
	uint32_t max;
};
//"顶向量“中的元素
struct p_dim2 {
	uint32_t head;
};
//一次合并的记录 —— 也就是标准 linkage 矩阵的一行。
//
//节点编号沿用 scipy 的约定：0..n-1 是原始对象，n+i 是第 i 次合并产生的簇。
//这里**组号天然就是节点 id**：初始 n 个组的 id 是 0..n-1，之后每次合并取
//new_grp = grp.size()，所以第 i 次合并产生的组号恰好是 n+i。
struct merge_rec {
	uint32_t left;      // 左子节点
	uint32_t right;     // 右子节点
	double sim;         // 合并时的相似度（该簇的"高度"，越大越先合并）
	double cohesion;    // 合并后整个簇的组内平均相似度（沿用原来输出里 sim= 的含义）
	std::size_t size;   // 合并后的簇大小
};

//聚类过程中的活组：组内平均相似度 + 成员数。
//
//注意这里**不再保存成员列表**：成员名由合并树（merges）+ object_names 提供，
//输出时经 cut_tree/object_name 取得，所以原来那份 vector<string> 只被用来取 size()。
//去掉之后省掉了每次合并复制成员列表的开销，以及 2n 个 vector 对象的常驻内存
//（n=200k 时约 10 MB）。
struct grp_s {
	double sim;
	std::size_t n_members;
};
#endif

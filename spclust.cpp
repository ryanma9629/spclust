#include "spclust.h"
#include "readdb.h"
#include <utility>

//“左向量”，记录每个横向链的头、值最大节点
vector<p_dim1> vec_dim1;
//“顶向量”，记录每个纵向链的头
vector<p_dim2> vec_dim2;

// ---------------------------------------------------------------------------
// 节点池。
//
// 所有单元放在一块连续内存里（下标寻址），释放的节点用空闲链表回收。
// 换掉 malloc/free 本身几乎不省时间（实测 1.00x，见 README），**压小节点才省时间**：
//   sizeof(sp_unit) 48 -> 32 字节
//
// ASan 的手工毒化：池里的内存在整个生命周期内都可访问，ASan 原本就看不见
// "单元被释放之后又被读" —— 而那正是这个代码库里出现过两次的 bug 类型。
// 所以释放时毒化该节点、分配时解除毒化，UAF 依然会被 ASan 抓到。
// 池扩容换地址时，之前记在旧缓冲上的毒化标记会失效，按空闲链表重来一遍。
// ---------------------------------------------------------------------------
//池用 malloc/realloc 直接管，而不是 std::vector：vector 扩容要「新块 + 拷贝 + 释放旧块」，
//而 macOS 的分配器不会立刻把大块还给系统，于是每一次搬迁都留下一个常驻的中间块 ——
//实测 800k 边时按 1.25 倍扩要 150 MB（总共 26 MB 的节点！），按 2 倍扩也要 85 MB。
//realloc 对大块常常能原地扩张，既不拷贝也不产生中间块。
static sp_unit* arena = NULL;
static size_t arena_cap = 0;    // 容量（节点数）
static size_t arena_used = 0;   // 已经分配出去过的下标数（含 0 号哨兵）
//空闲下标单独存，不借用节点自己的字段 —— 否则遍历空闲表就要去读被毒化的节点，
//而毒化正是为了抓住"读已释放节点"这件事。
static vector<uint32_t> arena_free_list;
static size_t arena_live = 0;

#if defined(__SANITIZE_ADDRESS__)
# define SPCLUST_ASAN 1
#elif defined(__has_feature)
# if __has_feature(address_sanitizer)
#  define SPCLUST_ASAN 1
# endif
#endif

#ifdef SPCLUST_ASAN
# include <sanitizer/asan_interface.h>
# define SP_POISON(p)   __asan_poison_memory_region((p), sizeof(sp_unit))
# define SP_UNPOISON(p) __asan_unpoison_memory_region((p), sizeof(sp_unit))
# define SP_UNPOISON_RANGE(p, n) __asan_unpoison_memory_region((p), (n))
#else
# define SP_POISON(p)   ((void)0)
# define SP_UNPOISON(p) ((void)0)
# define SP_UNPOISON_RANGE(p, n) ((void)0)
#endif

static uint32_t node_alloc(void) {
	if (0 == arena_used) {
		arena = (sp_unit*)calloc(1, sizeof(sp_unit));   // 下标 0 保留为空哨兵，永不出借
		if (NULL == arena) {
			cerr << "Arena allocation error." << endl;
			exit(1);
		}
		arena_cap = 1;
		arena_used = 1;
	}
	if (!arena_free_list.empty()) {
		uint32_t i = arena_free_list.back();
		arena_free_list.pop_back();
		SP_UNPOISON(&arena[i]);
		arena_live++;
		return i;
	}
	if (arena_used == arena_cap) {
		size_t want = arena_cap * 2;
		//realloc 可能搬家，而搬家要读旧内容 —— 其中可能有被毒化的空闲节点，先整片解除
		SP_UNPOISON_RANGE(arena, arena_used * sizeof(sp_unit));
		sp_unit* p = (sp_unit*)realloc(arena, want * sizeof(sp_unit));
		if (NULL == p) {
			cerr << "Arena allocation error." << endl;
			exit(1);
		}
		arena = p;
		arena_cap = want;
		for (size_t k = 0; k < arena_free_list.size(); k++)
			SP_POISON(&arena[arena_free_list[k]]);
	}
	arena_live++;
	return (uint32_t)arena_used++;
}

static void node_free(uint32_t i) {
	arena_live--;
	arena_free_list.push_back(i);
	SP_POISON(&arena[i]);
}


// merge_grp 每次只需重算“与 group1/group2 相邻的那些行/列”，即 O(degree) 个下标，
// 但旧实现为此要扫描整个索引空间三次（两趟更新 + 一趟 max 修复）：实测这三次
// O(n) 全量扫描占了聚类段约 75% 的时间，而且每次合并还要新分配 3 个 O(n) 向量。
// 下面这几个缓冲跨合并复用，只复位本次触及过的位置，于是每步都是 O(degree)。
static vector<unsigned char> mark_flags;  // 位0=要重算相似度, 位1=max失效, 位2/位3=已记下对 group1/group2 的相似度
static vector<size_t> updated;            // 需要重算相似度的下标（行号或列号）
static vector<size_t> dirty_max;          // max 指针失效、需要重扫的行
static vector<size_t> touched;            // 本次置过位的位置，用于 O(k) 复位

//合并树。每次 merge_grp 成功合并一对就追加一条 —— 这就是层次聚类的完整输出，
//MIN_SIM_CLUST 从"停止条件"变成了"切树的位置"（见 main.cpp 的 cut_tree）。
vector<merge_rec> merges;
static size_t n_objects = 0;
static vector<string> object_names;   // id -> 对象名

//本次合并中，下标 x 到 group1 / group2 的相似度。相似度非负（读取端限制 [0,65535]），
//所以 -1 可以安全地表示"没有这条边"。
//
//为什么要缓存这两个值：更新每个下标原本要调两次 get_unit()，而 get_unit() 是沿链线性查找，
//实测占整个程序约一半时间。但这两条边本来就在我们已经遍历的那 4 条链上 ——
//对下标 x，S(x, group1) 这个单元在 x > group1 时位于"行 group1"链上，
//x < group1 时位于"列 group1"链上；对 group2 同理。收集阶段顺手记下来，更新阶段就
//完全不需要 get_unit() 了。代价是两张 O(组数) 的 double 表，且只复位本次触及过的位置。
static vector<double> sim_g1;
static vector<double> sim_g2;

enum { F_UPDATED = 1, F_DIRTY_MAX = 2, F_SIM_G1 = 4, F_SIM_G2 = 8 };

//给下标 x 置一位；返回 true 表示这一位本次第一次置上
static bool set_flag(size_t x, unsigned char bit) {
	if (x >= mark_flags.size()) {
		mark_flags.resize(x + 1, 0);
		sim_g1.resize(x + 1, -1.0);
		sim_g2.resize(x + 1, -1.0);
	}
	if (0 == mark_flags[x])
		touched.push_back(x);          // 登记以便本次结束时 O(k) 复位
	bool fresh = (0 == (mark_flags[x] & bit));
	mark_flags[x] |= bit;
	return fresh;
}

//标记下标 x “需要重算相似度”；返回 true 表示这是本次第一次标记（顺带去重）
static bool mark_updated(size_t x) { return set_flag(x, F_UPDATED); }

//标记下标 x 的 max 指针失效
static bool mark_dirty_max(size_t x) { return set_flag(x, F_DIRTY_MAX); }

//记录下标 x 到 group1 / group2 的相似度。**首次写入生效**：get_unit() 是沿链从头部
//（add_grp 头插，即最后插入的单元）找第一个匹配，首次写入对应同一个单元，
//所以重复 pair 的输入下行为也保持一致。
static void record_sim_g1(size_t x, double v) {
	set_flag(x, F_SIM_G1);
	if (sim_g1[x] < 0)
		sim_g1[x] = v;
}
static void record_sim_g2(size_t x, double v) {
	set_flag(x, F_SIM_G2);
	if (sim_g2[x] < 0)
		sim_g2[x] = v;
}

// ---------------------------------------------------------------------------
// 全局最大相似度：最大堆 + 版本号懒删除。
//
// P0-1 之后 get_max() 那次 O(组数) 扫描已经是聚类段 57~78% 的时间，换成堆后是 O(1)。
//
// 关键设计：**把 key 冻结在堆项里**（sim 存进项里，不实时去 vec_dim1[row].max 读）。
// 这一点是必须的：如果 key 实时读，一次合并里改掉 k 行的 key 之后，堆数组的次序就
// 与当前 key 不再对应，堆序被破坏；而 sift_up/sift_down 都假定“除待修的那一处之外
// 处处是有效堆”，所以逐个就地 sift 并不能把它修回来（这一点是实测踩出来的）。
//
// 行内 max 变化时就推进该行版本号并压入新项，旧项靠版本号作废。带来的好处：
//   * 不需要任何“项在堆里的位置”簿记，也不需要让更新幂等
//   * 不依赖调用顺序，可以在 add_grp / 摘除 / max 重扫处就地维护
//   * 堆操作完全不碰 vec_dim1，所以合并中途那些暂时悬垂的 max 指针与堆无关
// 作废项堆积超过活行数的若干倍时重建一次（每个活行只留一项），内存仍是 O(n)，
// 并且重建只按堆自身大小做，摊还 O(1)。
// ---------------------------------------------------------------------------
struct SimEntry {
	double sim;
	size_t row;
	size_t version;
};

static vector<SimEntry> heap;
static vector<size_t> row_version;   // 每行的当前版本号
static vector<char> row_live;        // 该行当前是否非空（用于控制重建阈值）
static size_t live_rows = 0;

//堆序：先比相似度，相等时行号小的算“更大”。
//这与旧实现“i 升序扫描、只在严格大于当前最大值时接管”的并列行为一致。
static bool entry_less(const SimEntry& a, const SimEntry& b) {
	if (a.sim != b.sim)
		return a.sim < b.sim;
	return a.row > b.row;
}

static void heap_swap(size_t i, size_t j) {
	SimEntry tmp = heap[i];
	heap[i] = heap[j];
	heap[j] = tmp;
}

static void heap_sift_up(size_t i) {
	while (i > 0) {
		size_t parent = (i - 1) / 2;
		if (!entry_less(heap[parent], heap[i]))
			break;
		heap_swap(parent, i);
		i = parent;
	}
}

static void heap_sift_down(size_t i) {
	size_t n = heap.size();
	while (true) {
		size_t left = 2 * i + 1, right = left + 1, top = i;
		if (left < n && entry_less(heap[top], heap[left]))
			top = left;
		if (right < n && entry_less(heap[top], heap[right]))
			top = right;
		if (top == i)
			break;
		heap_swap(i, top);
		i = top;
	}
}

//行 row 的 max 变了（也可能变成空）：推进版本号；非空则压入新项。
//**每一处写 vec_dim1[*].max 的地方都必须调用它**，否则那一行会被漏掉。
static void note_row_changed(size_t row) {
	if (row >= row_version.size()) {
		row_version.resize(row + 1, 0);
		row_live.resize(row + 1, 0);
	}
	row_version[row]++;

	bool live = (SP_NONE != vec_dim1[row].max);
	if (live != (0 != row_live[row])) {
		row_live[row] = live ? 1 : 0;
		live_rows += live ? 1 : -1;
	}
	if (!live)
		return;                      // 变空：只让旧项作废

	SimEntry e;
	e.sim = arena[vec_dim1[row].max].sim;
	e.row = row;
	e.version = row_version[row];
	heap.push_back(e);
	heap_sift_up(heap.size() - 1);

	//作废项堆得太多就重建：只保留每个活行的当前项，然后 Floyd 建堆。
	//重建代价按堆自身大小算，所以相对压入次数是摊还 O(1)，内存也稳定在 O(活行数)。
	if (heap.size() > 4 * live_rows + 64) {
		vector<SimEntry> fresh;
		fresh.reserve(heap.size());
		for (size_t i = 0; i < heap.size(); i++) {
			const SimEntry& e2 = heap[i];
			if (e2.row < row_version.size() && e2.version == row_version[e2.row])
				fresh.push_back(e2);
		}
		heap.swap(fresh);
		for (size_t i = heap.size() / 2; i-- > 0; )
			heap_sift_down(i);
	}
}

//将对象名称映射到初始组
//用 unordered_map 而不是 map：这一步是按名称找对象，不需要有序，只需要快。
//实测（800k 边 / 8 万对象的采样画像）字符串查找是初始化段的主要开销。
typedef unordered_map<string, size_t> obj_id_map;
obj_id_map map_obj_grp;
extern vector<grp_s> grp;

//将一个(group1, group2, sim)的节点加入结构，每次加入所在链表的头位置
uint32_t add_grp(size_t group1, size_t group2, double sim) {
	//始终保证group1小于group2
	if (group1 > group2)
		swap(group1, group2);

	//”左向量“和”顶向量“扩容
	if (vec_dim1.size() <= group1)
		vec_dim1.resize(group1 + 1);
	if (vec_dim2.size() <= group2)
		vec_dim2.resize(group2 + 1);

	//待加入的节点（从池里取下标）
	uint32_t unit = node_alloc();
	arena[unit].group1 = group1;
	arena[unit].group2 = group2;
	arena[unit].sim = sim;
	arena[unit].right = SP_NONE;
	arena[unit].left = SP_NONE;
	arena[unit].up = SP_NONE;
	arena[unit].down = SP_NONE;

	//如果该行的max不存在或者比待加入的节点小，更新”左向量“的max指针
	if (SP_NONE == vec_dim1[group1].max || arena[vec_dim1[group1].max].sim < sim) {
		vec_dim1[group1].max = unit;
		note_row_changed(group1);
	}

	//原来的横向链表头
	uint32_t prev_head = vec_dim1[group1].head;
	//新加入的节点成为链表头
	vec_dim1[group1].head = unit;
	//如果原链表头存在，右移一个位置，更新原链表头和新节点之间的链接关系
	if (SP_NONE != prev_head) {
		arena[unit].right = prev_head;
		arena[prev_head].left = unit;
	}

	//原来的纵向链表头，余同上
	prev_head = vec_dim2[group2].head;
	vec_dim2[group2].head = unit;
	if (SP_NONE != prev_head) {
		arena[unit].down = prev_head;
		arena[prev_head].up = unit;
	}
	return unit;
}

//获取全局结构中值最大的节点（堆顶，摊还 O(1)）
uint32_t get_max(void) {
	//丢掉被版本号作废的堆顶
	while (!heap.empty()) {
		const SimEntry& e = heap[0];
		if (e.row < row_version.size() && e.version == row_version[e.row])
			break;
		heap[0] = heap.back();
		heap.pop_back();
		if (!heap.empty())
			heap_sift_down(0);
	}
	if (heap.empty())
		return SP_NONE;
	uint32_t unit = vec_dim1[heap[0].row].max;
	//与原实现一致：相似度为 0 不算候选（原实现用 max 初值 0 + 严格 > 来筛）。
	//堆里放着所有非空行，所以堆顶 <= 0 就意味着不存在正相似度的候选。
	if (SP_NONE == unit || !(arena[unit].sim > 0))
		return SP_NONE;
	return unit;
}

//根据(group1, group2, sim)定位该节点
uint32_t get_unit(size_t group1, size_t group2) {
	//始终保证group1小于group2
	if (group1 > group2)
		swap(group1, group2);
	//如果越界
	if (group1 >= vec_dim1.size() || group2 >= vec_dim2.size())
		return SP_NONE;
	//沿着vec_dim1中group1那一行链表遍历节点
	uint32_t unit = vec_dim1[group1].head;
	while(SP_NONE != unit) {
		if (group2 == arena[unit].group2)
			break;
		else
			unit = arena[unit].right;
	}
	return unit;
}

// ---------------------------------------------------------------------------
// 单链接（LINKAGE_MTD == 1）走 Kruskal，完全不建稀疏结构。
//
// 单链接的簇间相似度 = 两个簇之间所有**原始边**的最大值（max 不稀释），
// 所以整棵合并树等价于"把边按相似度降序排序，再用并查集逐个合并"（Kruskal / MST）。
//
// 这避开了原本的稠密化：max 不稀释 => 簇级相似度图会变成原图的**传递闭包**，
// n=40000 / 400k 边实测要建 93,418,345 个单元（是边数的 234 倍），7.35 s。
// Kruskal 只用 m 条边加一次排序。
//
// 并列（相似度相同的多条边）下合并**顺序**可能与贪心实现不同，但单链接在两个实现里
// 都等于"只保留 sim >= X 的边之后各连通分量"，与顺序无关；所以任何阈值处切出来的
// 分组、以及每个簇的紧致度（簇内最大边）都完全一致。
// ---------------------------------------------------------------------------
struct krusk_edge {
	uint32_t a;
	uint32_t b;
	double sim;
};

static vector<krusk_edge> krusk_edges;
static vector<uint32_t> krusk_parent;

static bool krusk_order(const krusk_edge& x, const krusk_edge& y) {
	if (x.sim != y.sim)
		return x.sim > y.sim;      // 相似度降序
	if (x.a != y.a)
		return x.a < y.a;          // 并列时按 id 定序，只为了确定性
	return x.b < y.b;
}

static uint32_t krusk_find(uint32_t x) {
	while (krusk_parent[x] != x) {
		krusk_parent[x] = krusk_parent[krusk_parent[x]];   // 路径减半
		x = krusk_parent[x];
	}
	return x;
}

void clust_kruskal(double floor) {
	sort(krusk_edges.begin(), krusk_edges.end(), krusk_order);

	const size_t n = n_objects;
	//集合的代表就是"代表这个簇的节点 id"：初始是对象自己，每次合并产生 n+i
	krusk_parent.resize(n);
	for (size_t i = 0; i < n; i++)
		krusk_parent[i] = (uint32_t)i;
	vector<size_t> size_of(n, 1);
	vector<double> cohesion_of(n, 0.0);

	size_t next_id = n;
	for (size_t i = 0; i < krusk_edges.size(); i++) {
		const krusk_edge& e = krusk_edges[i];
		if (e.sim < floor || !(e.sim > 0))
			break;                  // 到 floor 了（或没有正相似度的边了）
		uint32_t ra = krusk_find(e.a);
		uint32_t rb = krusk_find(e.b);
		if (ra == rb)
			continue;               // 已经在同一个簇里（重复边或成环）

		uint32_t node = (uint32_t)next_id++;
		merge_rec rec;
		rec.left = ra;
		rec.right = rb;
		rec.sim = e.sim;
		//紧致度 = 簇内最大边 = 原来的 max(max(g1,g2), max_sim)
		rec.cohesion = max(max(cohesion_of[ra], cohesion_of[rb]), e.sim);
		rec.size = size_of[ra] + size_of[rb];
		merges.push_back(rec);

		//喂给测试 runner 的合并日志（格式与贪心路径一致）
		if (DEBUG >= 2)
			cout << "\tGroup " << ra << " and group " << rb <<
				" are merged into group " << node << " (sim=" << setprecision(12)
				<< e.sim << ")." << endl;

		//两个集合的代表都指向新节点
		krusk_parent[ra] = node;
		krusk_parent[rb] = node;
		krusk_parent.push_back(node);
		cohesion_of.push_back(rec.cohesion);
		size_of.push_back(rec.size);

		if (DEBUG >= 1 && 0 == merges.size() % DEBUG_FREQ_CLUST)
			cout << "\t" << merges.size() << " merges, sim=" << e.sim << "." << endl;
	}
}

//把一条边加入结构。原来叫 init_grp()：先把全部边读进 vector<sp_sim> sims，
//再在这里遍历。现在由 readdb 边读边回调，那份 O(边数) 的缓存直接不存在了
//（实测 1M 条边、长对象名时可省 169 MB 峰值内存）。
//obj1/obj2 是按值传进来的：readdb 把字符串 move 进参数，这里再 move 进 map，
//全程零字符串拷贝（原来要先拷进 sims，再拷进 map，再拷进成员列表）。
void init_add_edge(string obj1, string obj2, double sim) {
	//每条边 2 次哈希查找（原来 4 次树查找），字符串只在哈希冲突时才比较。
	//id 仍然是"首次出现顺序"，所以结果不变。
	pair<obj_id_map::iterator, bool> o1 =
		map_obj_grp.insert(pair<string, size_t>(std::move(obj1), map_obj_grp.size()));
	if (o1.second) {
		//映射为当前最大的组号+1（就是插入前的 size()），并让”分组向量“记录obj1自成的新组
		grp.resize(map_obj_grp.size());
		grp[map_obj_grp.size() - 1].sim = 0;
		grp[map_obj_grp.size() - 1].objs.push_back(o1.first->first);
		object_names.push_back(o1.first->first);   // id -> 名字，供树输出用
	}

	pair<obj_id_map::iterator, bool> o2 =
		map_obj_grp.insert(pair<string, size_t>(std::move(obj2), map_obj_grp.size()));
	if (o2.second) {
		grp.resize(map_obj_grp.size());
		grp[map_obj_grp.size() - 1].sim = 0;
		grp[map_obj_grp.size() - 1].objs.push_back(o2.first->first);
		object_names.push_back(o2.first->first);
	}

	if (1 == LINKAGE_MTD) {
		//单链接不需要稀疏结构，只攒边表（后面走 Kruskal）。这一支是编译期常量，
		//所以 add_grp 那条路根本不会被执行到。
		krusk_edge e;
		e.a = (uint32_t)o1.first->second;
		e.b = (uint32_t)o2.first->second;
		e.sim = sim;
		krusk_edges.push_back(e);
	} else
		//把组信息而不是对象信息加入全局结构
		add_grp(o1.first->second, o2.first->second, sim);
}

//将group1和group2两组聚合，形成新组grp_k，更新其它组到grp_k的相似度
size_t merge_grp(size_t group1, size_t group2, double max_sim) {
	//“左向量”的长度，由于把两个组合并过程中该向量长度会改变，用dim1_size记录合并前的长度，不能用vec_dim1.size()
	size_t dim1_size = vec_dim1.size();
	//聚合之前group1和group2两组的对象个数
	size_t group1_size = grp[group1].objs.size();
	size_t group2_size = grp[group2].objs.size();
	//”分组向量“的长度
	size_t grp_size = grp.size();

	//收集需要重算相似度的行/列，以及 max 指针失效的行。只遍历被标记的下标。

	//所有形如(group1, ..., ...)链上节点所在的列都需要更新
	uint32_t unit = vec_dim1[group1].head;
	while (SP_NONE != unit) {
		if (group2 != arena[unit].group2) {
			record_sim_g1(arena[unit].group2, arena[unit].sim);   // 单元(group1,x) 就是 S(x,group1)
			if (mark_updated(arena[unit].group2))
				updated.push_back(arena[unit].group2);
		}
		unit = arena[unit].right;
	}

	//所有形如(group2, ..., ...)链上节点所在的列都需要更新，如果该链存在
	if (dim1_size > group2) {
		unit = vec_dim1[group2].head;
		while (SP_NONE != unit) {
			record_sim_g2(arena[unit].group2, arena[unit].sim);   // 单元(group2,x) 就是 S(x,group2)
			if (mark_updated(arena[unit].group2))
				updated.push_back(arena[unit].group2);
			unit = arena[unit].right;
		}
	}

	//所有形如(..., group1, ...)链上节点所在的行都需要更新
	unit = vec_dim2[group1].head;
	while (SP_NONE != unit) {
		record_sim_g1(arena[unit].group1, arena[unit].sim);       // 单元(x,group1) 就是 S(x,group1)
		if (mark_updated(arena[unit].group1))
			updated.push_back(arena[unit].group1);
		if (unit == vec_dim1[arena[unit].group1].max && mark_dirty_max(arena[unit].group1))
			dirty_max.push_back(arena[unit].group1);
		unit = arena[unit].down;
	}

	//所有形如(..., group2, ...)链上节点所在的行都需要更新
	unit = vec_dim2[group2].head;
	while (SP_NONE != unit) {
		if (group1 != arena[unit].group1) {
			record_sim_g2(arena[unit].group1, arena[unit].sim);   // 单元(x,group2) 就是 S(x,group2)
			if (mark_updated(arena[unit].group1))
				updated.push_back(arena[unit].group1);
			if (unit == vec_dim1[arena[unit].group1].max && mark_dirty_max(arena[unit].group1))
				dirty_max.push_back(arena[unit].group1);
		}
		unit = arena[unit].down;
	}

	//排成升序，保持与旧实现（两趟 for i = 0..dim 扫描标记）一致的处理顺序
	sort(updated.begin(), updated.end());

	//重算这些下标到新组的相似度。
	//旧实现分“行”和“列”两趟遍历，两趟的循环体逐字相同：对同一个下标 x，
	//都是 get_unit(x, group1) 和 get_unit(x, group2)。所以 x 同时落在两趟里时，
	//旧实现会为同一个 (x, grp_size) 调用两次 add_grp、插入两个数值相同的重复单元。
	//合并成一趟后这些重复单元不再产生（用例 20_overlap_update.csv 就是这个配置）。
	for (size_t k = 0; k < updated.size(); k++) {
		size_t i = updated[k];
		//收集阶段已经把这 4 条链上的值顺手记下来了，这里直接读，不再做 get_unit() 的线性查找。
		//-1 表示没有这条边，等价于原来的 get_unit 返回 NULL -> 0。
		double sim1 = (sim_g1[i] < 0) ? 0 : sim_g1[i];
		double sim2 = (sim_g2[i] < 0) ? 0 : sim_g2[i];

		//Average Linkage
		if (0 == LINKAGE_MTD) {
			double sim = (sim1 * group1_size + sim2 * group2_size) / 
				(group1_size + group2_size);
			if (sim > 0)
				add_grp(i, grp_size, sim);
		//Complete Linkage
		} else if (1 == LINKAGE_MTD) {
			//和另外两个分支一样，只留下 > 0 的候选。值为 0 意味着两个簇之间根本没有边，
			//而 get_max() 本来就不会选相似度 <= 0 的单元，所以去掉这些单元不改变合并顺序，
			//却能把矩阵从"每次合并都加满"变成保持稀疏 —— 实测 n=40000 时差 116 倍
			//（8.71s vs 0.075s）。
			double sim = max(sim1, sim2);
			if (sim > 0)
				add_grp(i, grp_size, sim);
		}
		//Single Linkage
		else if (2 == LINKAGE_MTD)
			if (min(sim1, sim2) > 0)
				add_grp(i, grp_size, min(sim1, sim2));
	}

	size_t new_grp = grp.size();
	//”分组向量“扩容
	grp.resize(new_grp + 1);
	//扩容后的空间先存入group1内所有对象
	grp[new_grp].objs = grp[group1].objs;
	//再存入group2内所有对象（append）
	grp[new_grp].objs.insert(grp[new_grp].objs.end(), 
		grp[group2].objs.begin(), grp[group2].objs.end());
	if (0 == LINKAGE_MTD)
		grp[new_grp].sim = (grp[group1].sim * group1_size * (group1_size - 1)
			+ grp[group2].sim * group2_size * (group2_size - 1) + group1_size * group2_size * max_sim * 2) 
			/ ((group1_size + group2_size) * (group1_size + group2_size - 1));
	else if (1 == LINKAGE_MTD)
		grp[new_grp].sim = max(max(grp[group1].sim, grp[group2].sim), max_sim);
	else if (2 == LINKAGE_MTD)
		grp[new_grp].sim = min(min(grp[group1].sim, grp[group2].sim), max_sim);

	//清除group1和group2两组内的对象
	grp[group1].objs.clear();
	grp[group1].sim = 0;
	grp[group2].objs.clear();
	grp[group2].sim = 0;

	//定位(group1, group2, ...)这个点
	uint32_t unit_max = get_unit(group1, group2);
	uint32_t tmp = SP_NONE;

	//更新链表的链接关系
	//删除(group1, ..., ...)一整行
	unit = vec_dim1[group1].head;
	while (SP_NONE != unit) {
		//先不删除(group1,group2,sim)这个节点
		if (unit != unit_max) {
			//如果该节点上方没有节点（列链表头）
			if (SP_NONE == arena[unit].up) {
				//如果该节点下方也没有节点（列链表尾）
				if (SP_NONE == arena[unit].down)
					//该节点在列链表上是个孤立节点，直接把”顶向量“头记录清空
					vec_dim2[arena[unit].group2].head = SP_NONE;
				else {
					//”顶向量“的头记录更新
					vec_dim2[arena[unit].group2].head = arena[unit].down;
					//该节点的下方节点变成整个列链表的头节点
					arena[arena[unit].down].up = SP_NONE;
				}
			} else {
				//该节点的上方节点，直接向下链接该节点的下方节点
				arena[arena[unit].up].down = arena[unit].down;
				//如果该节点有下方节点
				if (SP_NONE != arena[unit].down)
					//下方节点直接链接上方节点
					arena[arena[unit].down].up = arena[unit].up;
			}
			//把沿着横向的更新路径推进一个节点
			tmp = arena[unit].right;
			//释放当前节点
			node_free(unit);
			unit = tmp;
		} else
			unit = arena[unit].right;
	}

	//删除(group2, ..., ...)一整行，这行不一定存在需要先判断
	if (dim1_size > group2) {
		unit = vec_dim1[group2].head;
		while (SP_NONE != unit) {
			if (SP_NONE == arena[unit].up) {
				if (SP_NONE == arena[unit].down)
					vec_dim2[arena[unit].group2].head = SP_NONE;
				else {
					vec_dim2[arena[unit].group2].head = arena[unit].down;
					arena[arena[unit].down].up = SP_NONE;
				}
			} else {
				arena[arena[unit].up].down = arena[unit].down;
				if (SP_NONE != arena[unit].down)
					arena[arena[unit].down].up = arena[unit].up;
			}
			tmp = arena[unit].right;
			node_free(unit);
			unit = tmp;
		}
	}

	//删除(..., group1, ...)一整列
	unit = vec_dim2[group1].head;
	while (SP_NONE != unit) {
		if (SP_NONE == arena[unit].left) {
			if (SP_NONE == arena[unit].right)
				vec_dim1[arena[unit].group1].head = SP_NONE;
			else {
				vec_dim1[arena[unit].group1].head = arena[unit].right;
				arena[arena[unit].right].left = SP_NONE;
			}
		} else {
			arena[arena[unit].left].right = arena[unit].right;
			if (SP_NONE != arena[unit].right)
				arena[arena[unit].right].left = arena[unit].left;
		}
		tmp = arena[unit].down;
		node_free(unit);
		unit = tmp;
	}

	//删除(..., group2, ...)一整列
	unit = vec_dim2[group2].head;
	while (SP_NONE != unit) {
		if (unit != unit_max) {
			if (SP_NONE == arena[unit].left) {
				if (SP_NONE == arena[unit].right)
					vec_dim1[arena[unit].group1].head = SP_NONE;
				else {
					vec_dim1[arena[unit].group1].head = arena[unit].right;
					arena[arena[unit].right].left = SP_NONE;
				}
			} else {
				arena[arena[unit].left].right = arena[unit].right;
				if (SP_NONE != arena[unit].right)
					arena[arena[unit].right].left = arena[unit].left;
			}
			tmp = arena[unit].down;
			node_free(unit);
			unit = tmp;
		} else
			unit = arena[unit].down;
	}

	//最后删除(group1,group2,sim)节点
	node_free(unit_max);
	unit_max = SP_NONE;
	
	//”左向量“中删除下标为group1,group2（如果存在）指针
	vec_dim1[group1].head = SP_NONE;
	vec_dim1[group1].max = SP_NONE;
	note_row_changed(group1);
	if (dim1_size > group2) {
		vec_dim1[group2].head = SP_NONE;
		vec_dim1[group2].max = SP_NONE;
		note_row_changed(group2);
	}
	//”顶向量“中删除下标为group1,group2的指针
	vec_dim2[group1].head = SP_NONE;
	vec_dim2[group2].head = SP_NONE;

	//遍历需要修复 max 指针的行（升序），只遍历本次失效的那几行
	sort(dirty_max.begin(), dirty_max.end());
	for (size_t k = 0; k < dirty_max.size(); k++) {
		size_t i = dirty_max[k];
		uint32_t u = vec_dim1[i].head;
		uint32_t best = SP_NONE;
		//从左到右遍历，找到最大值
		while (SP_NONE != u) {
			if (SP_NONE == best || arena[u].sim > arena[best].sim)
				best = u;
			u = arena[u].right;
		}
		//必须无条件赋值。旧写法是 `if (arena[u].sim > max)`，max 初值 0，
		//一旦整行相似度都是 0 就永远不会赋值，vec_dim1[i].max 会保留一个已经被 free 的指针，
		//下一次 get_max() 就是 heap-use-after-free（见 tests/cases/22_zero_row_max.csv）。
		//best 仍是“从左到右第一个取得最大值的单元”，并列行为与旧写法一致。
		vec_dim1[i].max = best;
		note_row_changed(i);
	}

	//复位本次触及过的标记与缓存值。缓冲跨合并复用，所以这里既不再分配、也不清零整个 O(n) 内存。
	for (size_t k = 0; k < touched.size(); k++) {
		mark_flags[touched[k]] = 0;
		sim_g1[touched[k]] = -1.0;
		sim_g2[touched[k]] = -1.0;
	}
	touched.clear();
	updated.clear();
	dirty_max.clear();

	return grp.size();
}

//释放结构占用内存
void free_mem(void) {
	//节点池是一整块连续内存，整个还掉即可，不需要逐节点遍历释放。
	//先把毒化标记整体解除，别把它们留给已经还给分配器的内存。
	if (NULL != arena)
		SP_UNPOISON_RANGE(arena, arena_used * sizeof(sp_unit));
	free(arena);
	arena = NULL;
	arena_cap = 0;
	arena_used = 0;
	arena_free_list.clear();
	vector<uint32_t>().swap(arena_free_list);
	arena_live = 0;

	//清除“左向量”和“顶向量”的指针
	for (size_t i = 0; i < vec_dim1.size(); i++) {
		vec_dim1[i].head = SP_NONE;
		vec_dim1[i].max = SP_NONE;
	}
	for (size_t i = 0; i < vec_dim2.size(); i++)
		vec_dim2[i].head = SP_NONE;

	vec_dim1.clear();
	vec_dim2.clear();

	//合并过程中复用的缓冲也一并释放
	mark_flags.clear();
	sim_g1.clear();
	sim_g2.clear();
	updated.clear();
	dirty_max.clear();
	touched.clear();
	krusk_edges.clear();
	vector<krusk_edge>().swap(krusk_edges);
	krusk_parent.clear();
	vector<uint32_t>().swap(krusk_parent);
	heap.clear();
	row_version.clear();
	row_live.clear();
	live_rows = 0;
}

size_t tree_object_count(void) {
	return n_objects;
}

//节点 id -> 对象名
const string& object_name(size_t id) {
	return object_names[id];
}

//聚类主过程
void clust_grp(double floor) {
	//结构在 readdb 边读边回调时就已经建好了（见 init_add_edge），这里只剩聚类
	clock_t s_time = clock();
	clock_t e_time = 0;
	if (DEBUG >= 1)
		cout << "Clustering..." << endl;

	//初始组（对象）数 = 树的叶子数
	n_objects = grp.size();

	if (1 == LINKAGE_MTD) {
		//单链接：边按相似度降序 + 并查集，O(m log m)
		clust_kruskal(floor);
		free_mem();
		e_time = clock();
		if (DEBUG >= 1)
			cout << "Done. (" << (double)(e_time - s_time) / CLOCKS_PER_SEC << "s)" << endl << endl;
		return;
	}

	size_t grp_size = grp.size();
	//合并到最大相似度 < floor 为止。floor 由调用方给：
	//  floor = 切树阈值  -> 成本与"直接在阈值处停下"完全一样（默认就是这样）
	//  floor <= 0        -> 完整树，代价是多做很多次合并（实测 2.2~4.2x，见 README）
	//记录下来的 merges 就是树；MIN_SIM_CLUST 不再是停止条件，而是切树的默认位置。
	while (grp_size > 1) {
		uint32_t unit = get_max();
		if (SP_NONE == unit)
			break;                      // 没有正相似度的候选了，树到此到底
		double max_sim = arena[unit].sim;
		if (max_sim < floor)
			break;                      // 到 floor 了；想要更深的树就把 floor 调低
		size_t group1 = arena[unit].group1;
		size_t group2 = arena[unit].group2;
		//每次合并，实际组数减少1
		if (DEBUG >= 1 && 0 == grp_size % DEBUG_FREQ_CLUST)
			cout << "\t" << grp_size << " groups left with max_sim=" << max_sim << "." << endl;
		grp_size--;
		size_t new_grp = merge_grp(group1, group2, max_sim) - 1;

		//记录这次合并（组号即节点 id，见 main.h 的 merge_rec）
		merge_rec rec;
		rec.left = (uint32_t)group1;
		rec.right = (uint32_t)group2;
		rec.sim = max_sim;
		rec.cohesion = grp[new_grp].sim;
		rec.size = grp[new_grp].objs.size();
		merges.push_back(rec);

		if (DEBUG >= 2)
			cout << "\tGroup " << group1 << " and group " << group2 << 
				" are merged into group " << new_grp << " (sim=" << setprecision(12)
				<< max_sim << ")." << endl;
	}
	//聚类结束，释放结构占用内存
	free_mem();

	e_time = clock();
	if (DEBUG >= 1)
		cout << "Done. (" << (double)(e_time - s_time) / CLOCKS_PER_SEC << "s)" << endl << endl;
}
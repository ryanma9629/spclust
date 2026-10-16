#include "spclust.h"
#include "readdb.h"
#include <utility>

//“左向量”，记录每个横向链的头、值最大节点
vector<p_dim1> vec_dim1;
//“顶向量”，记录每个纵向链的头
vector<p_dim2> vec_dim2;

// merge_grp 每次只需重算“与 group1/group2 相邻的那些行/列”，即 O(degree) 个下标，
// 但旧实现为此要扫描整个索引空间三次（两趟更新 + 一趟 max 修复）：实测这三次
// O(n) 全量扫描占了聚类段约 75% 的时间，而且每次合并还要新分配 3 个 O(n) 向量。
// 下面这几个缓冲跨合并复用，只复位本次触及过的位置，于是每步都是 O(degree)。
static vector<unsigned char> mark_flags;  // 位0=要重算相似度, 位1=max失效, 位2/位3=已记下对 group1/group2 的相似度
static vector<size_t> updated;            // 需要重算相似度的下标（行号或列号）
static vector<size_t> dirty_max;          // max 指针失效、需要重扫的行
static vector<size_t> touched;            // 本次置过位的位置，用于 O(k) 复位

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

	bool live = (NULL != vec_dim1[row].max);
	if (live != (0 != row_live[row])) {
		row_live[row] = live ? 1 : 0;
		live_rows += live ? 1 : -1;
	}
	if (!live)
		return;                      // 变空：只让旧项作废

	SimEntry e;
	e.sim = vec_dim1[row].max->sim;
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
sp_unit* add_grp(size_t group1, size_t group2, double sim) {
	//始终保证group1小于group2
	if (group1 > group2)
		swap(group1, group2);

	//”左向量“和”顶向量“扩容
	if (vec_dim1.size() <= group1)
		vec_dim1.resize(group1 + 1);
	if (vec_dim2.size() <= group2)
		vec_dim2.resize(group2 + 1);

	//待加入的节点
	sp_unit* unit = (sp_unit*)malloc(sizeof(sp_unit));
	if (NULL == unit) {
		cerr << "Malloc error." << endl;
		exit(1);
	}
	unit->group1 = group1;
	unit->group2 = group2;
	unit->sim = sim;
	unit->right = NULL;
	unit->left = NULL;
	unit->up = NULL;
	unit->down = NULL;

	//如果该行的max不存在或者比待加入的节点小，更新”左向量“的max指针
	if (NULL == vec_dim1[group1].max || vec_dim1[group1].max->sim < sim) {
		vec_dim1[group1].max = unit;
		note_row_changed(group1);
	}

	//原来的横向链表头
	sp_unit* prev_head = vec_dim1[group1].head;
	//新加入的节点成为链表头
	vec_dim1[group1].head = unit;
	//如果原链表头存在，右移一个位置，更新原链表头和新节点之间的链接关系
	if (NULL != prev_head) {
		unit->right = prev_head;
		prev_head->left = unit;
	}

	//原来的纵向链表头，余同上
	prev_head = vec_dim2[group2].head;
	vec_dim2[group2].head = unit;
	if (NULL != prev_head) {
		unit->down = prev_head;
		prev_head->up = unit;
	}
	return unit;
}

//获取全局结构中值最大的节点（堆顶，摊还 O(1)）
sp_unit* get_max(void) {
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
		return NULL;
	sp_unit* unit = vec_dim1[heap[0].row].max;
	//与原实现一致：相似度为 0 不算候选（原实现用 max 初值 0 + 严格 > 来筛）。
	//堆里放着所有非空行，所以堆顶 <= 0 就意味着不存在正相似度的候选。
	if (NULL == unit || !(unit->sim > 0))
		return NULL;
	return unit;
}

//根据(group1, group2, sim)定位该节点
sp_unit* get_unit(size_t group1, size_t group2) {
	//始终保证group1小于group2
	if (group1 > group2)
		swap(group1, group2);
	//如果越界
	if (group1 >= vec_dim1.size() || group2 >= vec_dim2.size())
		return NULL;
	//沿着vec_dim1中group1那一行链表遍历节点
	sp_unit* unit = vec_dim1[group1].head;
	while(NULL != unit) {
		if (group2 == unit->group2)
			break;
		else
			unit = unit->right;
	}
	return unit;
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
	}

	pair<obj_id_map::iterator, bool> o2 =
		map_obj_grp.insert(pair<string, size_t>(std::move(obj2), map_obj_grp.size()));
	if (o2.second) {
		grp.resize(map_obj_grp.size());
		grp[map_obj_grp.size() - 1].sim = 0;
		grp[map_obj_grp.size() - 1].objs.push_back(o2.first->first);
	}

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
	sp_unit* unit = vec_dim1[group1].head;
	while (NULL != unit) {
		if (group2 != unit->group2) {
			record_sim_g1(unit->group2, unit->sim);   // 单元(group1,x) 就是 S(x,group1)
			if (mark_updated(unit->group2))
				updated.push_back(unit->group2);
		}
		unit = unit->right;
	}

	//所有形如(group2, ..., ...)链上节点所在的列都需要更新，如果该链存在
	if (dim1_size > group2) {
		unit = vec_dim1[group2].head;
		while (NULL != unit) {
			record_sim_g2(unit->group2, unit->sim);   // 单元(group2,x) 就是 S(x,group2)
			if (mark_updated(unit->group2))
				updated.push_back(unit->group2);
			unit = unit->right;
		}
	}

	//所有形如(..., group1, ...)链上节点所在的行都需要更新
	unit = vec_dim2[group1].head;
	while (NULL != unit) {
		record_sim_g1(unit->group1, unit->sim);       // 单元(x,group1) 就是 S(x,group1)
		if (mark_updated(unit->group1))
			updated.push_back(unit->group1);
		if (unit == vec_dim1[unit->group1].max && mark_dirty_max(unit->group1))
			dirty_max.push_back(unit->group1);
		unit = unit->down;
	}

	//所有形如(..., group2, ...)链上节点所在的行都需要更新
	unit = vec_dim2[group2].head;
	while (NULL != unit) {
		if (group1 != unit->group1) {
			record_sim_g2(unit->group1, unit->sim);   // 单元(x,group2) 就是 S(x,group2)
			if (mark_updated(unit->group1))
				updated.push_back(unit->group1);
			if (unit == vec_dim1[unit->group1].max && mark_dirty_max(unit->group1))
				dirty_max.push_back(unit->group1);
		}
		unit = unit->down;
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
		} else if (1 == LINKAGE_MTD)
			add_grp(i, grp_size, max(sim1, sim2));
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
	sp_unit* unit_max = get_unit(group1, group2);
	sp_unit* tmp = NULL;

	//更新链表的链接关系
	//删除(group1, ..., ...)一整行
	unit = vec_dim1[group1].head;
	while (NULL != unit) {
		//先不删除(group1,group2,sim)这个节点
		if (unit != unit_max) {
			//如果该节点上方没有节点（列链表头）
			if (NULL == unit->up) {
				//如果该节点下方也没有节点（列链表尾）
				if (NULL == unit->down)
					//该节点在列链表上是个孤立节点，直接把”顶向量“头记录清空
					vec_dim2[unit->group2].head = NULL;
				else {
					//”顶向量“的头记录更新
					vec_dim2[unit->group2].head = unit->down;
					//该节点的下方节点变成整个列链表的头节点
					unit->down->up = NULL;
				}
			} else {
				//该节点的上方节点，直接向下链接该节点的下方节点
				unit->up->down = unit->down;
				//如果该节点有下方节点
				if (NULL != unit->down)
					//下方节点直接链接上方节点
					unit->down->up = unit->up;
			}
			//把沿着横向的更新路径推进一个节点
			tmp = unit->right;
			//释放当前节点
			free(unit);
			unit = tmp;
		} else
			unit = unit->right;
	}

	//删除(group2, ..., ...)一整行，这行不一定存在需要先判断
	if (dim1_size > group2) {
		unit = vec_dim1[group2].head;
		while (NULL != unit) {
			if (NULL == unit->up) {
				if (NULL == unit->down)
					vec_dim2[unit->group2].head = NULL;
				else {
					vec_dim2[unit->group2].head = unit->down;
					unit->down->up = NULL;
				}
			} else {
				unit->up->down = unit->down;
				if (NULL != unit->down)
					unit->down->up = unit->up;
			}
			tmp = unit->right;
			free(unit);
			unit = tmp;
		}
	}

	//删除(..., group1, ...)一整列
	unit = vec_dim2[group1].head;
	while (NULL != unit) {
		if (NULL == unit->left) {
			if (NULL == unit->right)
				vec_dim1[unit->group1].head = NULL;
			else {
				vec_dim1[unit->group1].head = unit->right;
				unit->right->left = NULL;
			}
		} else {
			unit->left->right = unit->right;
			if (NULL != unit->right)
				unit->right->left = unit->left;
		}
		tmp = unit->down;
		free(unit);
		unit = tmp;
	}

	//删除(..., group2, ...)一整列
	unit = vec_dim2[group2].head;
	while (NULL != unit) {
		if (unit != unit_max) {
			if (NULL == unit->left) {
				if (NULL == unit->right)
					vec_dim1[unit->group1].head = NULL;
				else {
					vec_dim1[unit->group1].head = unit->right;
					unit->right->left = NULL;
				}
			} else {
				unit->left->right = unit->right;
				if (NULL != unit->right)
					unit->right->left = unit->left;
			}
			tmp = unit->down;
			free(unit);
			unit = tmp;
		} else
			unit = unit->down;
	}

	//最后删除(group1,group2,sim)节点
	free(unit_max);
	unit_max = NULL;
	
	//”左向量“中删除下标为group1,group2（如果存在）指针
	vec_dim1[group1].head = NULL;
	vec_dim1[group1].max = NULL;
	note_row_changed(group1);
	if (dim1_size > group2) {
		vec_dim1[group2].head = NULL;
		vec_dim1[group2].max = NULL;
		note_row_changed(group2);
	}
	//”顶向量“中删除下标为group1,group2的指针
	vec_dim2[group1].head = NULL;
	vec_dim2[group2].head = NULL;

	//遍历需要修复 max 指针的行（升序），只遍历本次失效的那几行
	sort(dirty_max.begin(), dirty_max.end());
	for (size_t k = 0; k < dirty_max.size(); k++) {
		size_t i = dirty_max[k];
		sp_unit* u = vec_dim1[i].head;
		sp_unit* best = NULL;
		//从左到右遍历，找到最大值
		while (NULL != u) {
			if (NULL == best || u->sim > best->sim)
				best = u;
			u = u->right;
		}
		//必须无条件赋值。旧写法是 `if (u->sim > max)`，max 初值 0，
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
	sp_unit* tmp = NULL;
	sp_unit* unit = NULL;
	//沿着”左向量“逐行释放链表空间，并清除“左向量”的指针
	for (size_t i = 0; i < vec_dim1.size(); i++) {
		unit = vec_dim1[i].head;
		while (NULL != unit) {
			tmp = unit->right;
			free(unit);
			unit = tmp;
		}
		vec_dim1[i].head = NULL;
		vec_dim1[i].max = NULL;
	}
	//清除“顶向量”的指针
	for (size_t i = 0; i < vec_dim2.size(); i++)
		vec_dim2[i].head = NULL;

	vec_dim1.clear();
	vec_dim2.clear();

	//合并过程中复用的缓冲也一并释放
	mark_flags.clear();
	sim_g1.clear();
	sim_g2.clear();
	updated.clear();
	dirty_max.clear();
	touched.clear();
	heap.clear();
	row_version.clear();
	row_live.clear();
	live_rows = 0;
}

//聚类主过程
void clust_grp(void) {
	//结构在 readdb 边读边回调时就已经建好了（见 init_add_edge），这里只剩聚类
	clock_t s_time = clock();
	clock_t e_time = 0;
	if (DEBUG >= 1)
		cout << "Clustering..." << endl;

	//初始组（对象）数
	size_t grp_size = grp.size();
	//如果没有全部合并成1组，或剩余最大相似度小于阈值 MIN_SIM_CLUST
	while (grp_size > 1) {
		sp_unit* unit = get_max();
		double max_sim = (NULL != unit) ? unit->sim : 0;
		if (NULL == unit || max_sim < MIN_SIM_CLUST)
			break;
		else {
			size_t group1 = unit->group1;
			size_t group2 = unit->group2;
			//每次合并，实际组数减少1
			if (DEBUG >= 1 && 0 == grp_size % DEBUG_FREQ_CLUST)
				cout << "\t" << grp_size << " groups left with max_sim=" << max_sim << "." << endl;
			grp_size--;
			size_t new_grp = merge_grp(group1, group2, max_sim) - 1;
			if (DEBUG >= 2)
				cout << "\tGroup " << group1 << " and group " << group2 << 
					" are merged into group " << new_grp << " (sim=" << setprecision(12)
					<< max_sim << ")." << endl;
		}
	}
	//聚类结束，释放结构占用内存
	free_mem();

	e_time = clock();
	if (DEBUG >= 1)
		cout << "Done. (" << (double)(e_time - s_time) / CLOCKS_PER_SEC << "s)" << endl << endl;
}
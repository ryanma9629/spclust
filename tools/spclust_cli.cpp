// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

#include "spclust/hac.h"
#include "spclust/scipy.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "cli_config.h"
#include "readdb.h"

using namespace std;
using spclust::Cut;
using spclust::Merge;

//本文件是库的第一个使用者：它负责**所有** IO 与展示 —— 读文件、给对象分配 id、
//把 id 映射回名字、计时、打印。库里没有字符串、没有标准输出、没有全局名字表。

namespace {

//本程序只跑一个聚类，所以这里放一个指向它的全局指针，供下面几个打印函数用。
//（库那边是一对象一份状态，不再有全局量；这里只是 CLI 自己图个方便。）
spclust::HAC* g_hac = NULL;


//切法。三种互斥，命令行里后给的那个生效：
//   --cut <X>        切在相似度 X 处（默认，X 取编译期 MIN_SIM_CLUST）
//   --n-clusters <K> 要正好 K 个簇
//   --min-size <M>   合并到没有小于 M 的簇为止
//
//合并树本身是完整算出来的（spclust::build），所以这个阈值只决定"从哪里切开"，
//换一个阈值不需要重新聚类。
enum CutMode { CUT_SIM, CUT_N, CUT_SIZE };
CutMode cut_mode = CUT_SIM;
double cut_sim = (double)MIN_SIM_CLUST;
size_t cut_n_clusters = 0;
size_t cut_min_size = 0;

//按当前模式切树
Cut do_cut(void) {
	switch (cut_mode) {
	case CUT_N:    return g_hac->cut_n_clusters(cut_n_clusters);
	case CUT_SIZE: return g_hac->cut_min_size(cut_min_size);
	default:       return g_hac->cut(cut_sim);
	}
}

//一句话描述当前切法（用于提示信息）
string cut_desc(void) {
	char buf[128];
	switch (cut_mode) {
	case CUT_N:    snprintf(buf, sizeof(buf), "n-clusters=%zu", cut_n_clusters); break;
	case CUT_SIZE: snprintf(buf, sizeof(buf), "min-size=%zu", cut_min_size); break;
	default:       snprintf(buf, sizeof(buf), "cut=%g", cut_sim); break;
	}
	return string(buf);
}
//树算到多低。默认 = 切树阈值（于是成本与"直接在阈值处停下"完全一样）；
//用 --floor 显式指定（0 = 完整树）或者把 --cut 定得比默认更低时，会自动跟着降。
bool floor_given = false;
double tree_floor = (double)MIN_SIM_CLUST;
bool show_tree = false;
int tree_depth = 6;          // --tree 的深度上限；0 表示不限
bool show_matrix = false;
bool show_scipy = false;                 // --scipy：导出 scipy linkage 矩阵
spclust::DistanceMap scipy_map = spclust::DistanceMap::Complement;
bool scipy_bridge = false;               // 森林补成一棵树（scipy 要求 n-1 行）

//链接方式。以前这是编译期常量（LINKAGE_MTD），要换就得重编译；现在库里是
//构造参数，所以一个二进制就能跑三种，也能用 --linkage 在命令行选。
spclust::Linkage linkage = spclust::Linkage::Average;

//对象 id -> 名字。库只认 id；名字映射是 CLI 的事。
vector<string> names;
unordered_map<string, uint32_t> name_to_id;

//把一条边交给库。id 按"首次出现顺序"分配 —— 与节点 id 的约定一致，
//所以树里的叶子 id 就是这里的下标。
void edge_to_library(string obj1, string obj2, double sim) {
	uint32_t id1, id2;
	pair<unordered_map<string, uint32_t>::iterator, bool> r1 =
		name_to_id.insert(make_pair(std::move(obj1), (uint32_t)name_to_id.size()));
	if (r1.second)
		names.push_back(r1.first->first);
	id1 = r1.first->second;
	pair<unordered_map<string, uint32_t>::iterator, bool> r2 =
		name_to_id.insert(make_pair(std::move(obj2), (uint32_t)name_to_id.size()));
	if (r2.second)
		names.push_back(r2.first->first);
	id2 = r2.first->second;
	g_hac->add_edge(id1, id2, sim);
}

// ---------------------------------------------------------------------------
// 输出
// ---------------------------------------------------------------------------

//扁平分组（把树按 cut_sim 切开）
void print_result(void) {
	//spclust.cpp 的 DEBUG=2 日志把 cout 的精度调到了 12，这里显式改回 6 位有效数字
	cout << setprecision(6);

	Cut c = do_cut();

	//只列出 size>=2 的簇，按"簇大小降序、紧致度降序"
	vector<size_t> clusters;
	for (size_t i = 0; i < c.members.size(); i++)
		if (c.members[i].size() >= 2)
			clusters.push_back(i);
	sort(clusters.begin(), clusters.end(), [&c](size_t a, size_t b) {
		if (c.members[a].size() != c.members[b].size())
			return c.members[a].size() > c.members[b].size();
		return c.cohesion[a] > c.cohesion[b];
	});

	cout << "Objects: " << g_hac->num_objects() << ", groups: " << clusters.size() << endl;
	if (clusters.empty()) {
		cout << "No group reached 2 objects (" << cut_desc() << ")." << endl;
	}
	for (size_t k = 0; k < clusters.size(); k++) {
		const vector<uint32_t>& m = c.members[clusters[k]];
		cout << "Group " << k << " (size=" << m.size() << ", sim=" << c.cohesion[clusters[k]] << "):";
		for (size_t j = 0; j < m.size(); j++)
			cout << " " << names[m[j]];
		cout << endl;
	}

	//未成簇（仍是单元素）的对象
	vector<string> single;
	for (size_t i = 0; i < c.members.size(); i++) {
		if (c.members[i].size() == 1)
			single.push_back(names[c.members[i][0]]);
	}
	sort(single.begin(), single.end());
	if (!single.empty()) {
		cout << "Unclustered (" << single.size() << "):";
		for (size_t i = 0; i < single.size(); i++)
			cout << " " << single[i];
		cout << endl;
	}
}

//树状图（人类可读）。深度受限，太深的地方折叠成一行，完整数据用 --matrix。
void print_subtree(size_t node, int depth, int max_depth) {
	const size_t n = g_hac->num_objects();
	for (int i = 0; i < depth; i++)
		cout << "  ";
	if (node < n) {
		cout << "- " << names[node] << endl;
		return;
	}
	const Merge& r = g_hac->tree()[node - n];
	cout << "+ merged at sim=" << r.similarity << " (size=" << r.size << ")" << endl;
	if (max_depth > 0 && depth + 1 > max_depth) {
		for (int i = 0; i <= depth; i++)
			cout << "  ";
		cout << "... (" << r.size << " objects, 用 --matrix 取完整数据)" << endl;
		return;
	}
	print_subtree(r.left, depth + 1, max_depth);
	print_subtree(r.right, depth + 1, max_depth);
}

void print_tree(void) {
	const vector<Merge>& tr = g_hac->tree();
	const size_t n = g_hac->num_objects();
	//根 = 没有当过任何一次合并的子节点的节点
	vector<char> is_child(n + tr.size(), 0);
	for (size_t i = 0; i < tr.size(); i++) {
		is_child[tr[i].left] = 1;
		is_child[tr[i].right] = 1;
	}
	cout << "Dendrogram (" << n << " objects, " << tr.size() << " merges; "
		<< "阈值只是切树位置，树本身是完整的):" << endl;
	const int max_depth = tree_depth;
	size_t roots = 0;
	for (size_t node = 0; node < n + tr.size(); node++) {
		if (!is_child[node]) {
			if (roots > 0)
				cout << endl;      // 多个根（互不相连的分量）之间空一行
			print_subtree(node, 0, max_depth);
			roots++;
		}
	}
	cout << endl << "（共 " << roots << " 个根：相似度为 0 的簇不会再合并）" << endl;
}

//linkage 矩阵：每行 left right height size，节点 id 沿用 scipy 的约定
//（0..n-1 是对象，n+i 是第 i 次合并）。注意 height 这里是**相似度**（越大越先合并），
//scipy 用的是距离。
void print_matrix(void) {
	const vector<Merge>& tr = g_hac->tree();
	const size_t n = g_hac->num_objects();
	cout << "# linkage matrix: left right height size" << endl;
	cout << "# 节点 id：0.." << (n ? n - 1 : 0) << " 是对象（按输入中首次出现的顺序），"
		<< n << "+i 是第 i 次合并。height 是**相似度**（越大越先合并），不是距离。" << endl;
	cout << setprecision(12);
	for (size_t i = 0; i < tr.size(); i++) {
		cout << tr[i].left << " " << tr[i].right << " "
			<< tr[i].similarity << " " << tr[i].size << endl;
	}
}

//scipy 兼容导出：每行 left right distance size。
//注意 height 是**距离**（越小越先合并），与 --matrix 的相似度相反。
void print_scipy(void) {
	spclust::ScipyOptions o;
	o.map = scipy_map;
	o.bridge = scipy_bridge;
	spclust::ScipyExport e = spclust::to_scipy(g_hac->tree(), g_hac->num_objects(), o);

	if (e.n_components > 1 && !e.bridged)
		cerr << "注意: 合并树是森林（" << e.n_components
			<< " 个分量，树里没有把它们并起来）—— scipy 要求恰好 n-1 行。"
			<< "用 --scipy-bridge 补上桥接行，或者用 average/max 链接。" << endl;
	if (!e.monotonic)
		cerr << "注意: 距离沿合并序列不是单调递增的（MinNonzero 不保证），"
			<< "scipy 的 dendrogram/fcluster 假定它单调。" << endl;

	cout << "# scipy linkage matrix: left right distance size" << endl;
	cout << "# n=" << e.n_objects << " rows=" << e.rows.size()
		<< " components=" << e.n_components << (e.bridged ? " (bridged)" : "") << endl;
	cout << setprecision(12);
	for (size_t i = 0; i < e.rows.size(); i++)
		cout << e.rows[i].left << " " << e.rows[i].right << " "
			<< e.rows[i].distance << " " << e.rows[i].size << endl;
}

void usage(const char* prog) {
	cerr << "用法: " << prog << " [similarity.csv|-] [选项]" << endl
	     << "  similarity.csv   输入边表；省略或给 - 表示从标准输入读" << endl
	     << "  --cut <X>        把合并树切在相似度 X 处（默认 " << MIN_SIM_CLUST
	     << "，即编译期 MIN_SIM_CLUST）" << endl
	     << "  --n-clusters <K> 要正好 K 个簇（需要完整树，见下）" << endl
	     << "  --min-size <M>   合并到没有小于 M 的簇为止（需要完整树，见下）" << endl
	     << "  --floor <X>      树算到多低（默认 = cut；0 = 完整树，代价是慢 2~4 倍）" << endl
	     << "  --tree           打印树状图" << endl
	     << "  --tree-depth <N> 树状图的深度上限（默认 " << 6 << "，0 = 不限；隐含 --tree）" << endl
	     << "  --linkage <M>    链接方式：0/average（默认）、1/max|single、2/min|min-nonzero" << endl
	     << "  --matrix         打印 linkage 矩阵（left right height size，height 是相似度）" << endl
	     << "  --scipy          导出 scipy.cluster.hierarchy 的 linkage 矩阵（height 是距离）" << endl
	     << "  --scipy-distance <M>  complement（默认，C-sim）| negate（-sim）| reciprocal（1/sim）" << endl
	     << "  --scipy-bridge   森林补成一棵树（桥接距离 = inf），scipy 要求恰好 n-1 行" << endl
	     << "  -h, --help       显示这个" << endl;
}

} // namespace

int main(int argc, char* argv[]) {
	const char* path = NULL;
	for (int i = 1; i < argc; i++) {
		const string a = argv[i];
		if ("--cut" == a) {
			if (++i >= argc) {
				cerr << "--cut 需要一个数值" << endl;
				return 2;
			}
			cut_sim = atof(argv[i]);
		} else if ("--n-clusters" == a) {
			if (++i >= argc) {
				cerr << "--n-clusters 需要一个数值" << endl;
				return 2;
			}
			cut_n_clusters = (size_t)strtoul(argv[i], NULL, 10);
			cut_mode = CUT_N;
		} else if ("--min-size" == a) {
			if (++i >= argc) {
				cerr << "--min-size 需要一个数值" << endl;
				return 2;
			}
			cut_min_size = (size_t)strtoul(argv[i], NULL, 10);
			cut_mode = CUT_SIZE;
		} else if ("--floor" == a) {
			if (++i >= argc) {
				cerr << "--floor 需要一个数值" << endl;
				return 2;
			}
			tree_floor = atof(argv[i]);
			floor_given = true;
		} else if ("--tree" == a) {
			show_tree = true;
		} else if ("--tree-depth" == a) {
			if (++i >= argc) {
				cerr << "--tree-depth 需要一个数值" << endl;
				return 2;
			}
			tree_depth = atoi(argv[i]);
			show_tree = true;        // 给了深度显然是想看树

		} else if ("--linkage" == a) {
			if (++i >= argc) {
				cerr << "--linkage 需要一个值" << endl;
				return 2;
			}
			const string v = argv[i];
			if ("0" == v || "average" == v)
				linkage = spclust::Linkage::Average;
			else if ("1" == v || "max" == v || "single" == v || "max-similarity" == v)
				linkage = spclust::Linkage::MaxSimilarity;
			else if ("2" == v || "min" == v || "min-nonzero" == v)
				linkage = spclust::Linkage::MinNonzero;
			else {
				cerr << "未知的 linkage: " << v << "（用 0/average、1/max、2/min）" << endl;
				return 2;
			}

		} else if ("--scipy" == a) {
			show_scipy = true;
		} else if ("--scipy-distance" == a) {
			if (++i >= argc) {
				cerr << "--scipy-distance 需要一个值" << endl;
				return 2;
			}
			const string v = argv[i];
			if ("complement" == v)
				scipy_map = spclust::DistanceMap::Complement;
			else if ("negate" == v)
				scipy_map = spclust::DistanceMap::Negate;
			else if ("reciprocal" == v)
				scipy_map = spclust::DistanceMap::Reciprocal;
			else {
				cerr << "未知的距离映射: " << v << "（complement/negate/reciprocal）" << endl;
				return 2;
			}
		} else if ("--scipy-bridge" == a) {
			scipy_bridge = true;
		} else if ("--matrix" == a) {
			show_matrix = true;
		} else if ("-h" == a || "--help" == a) {
			usage(argv[0]);
			return 0;
		} else if (a.size() > 1 && '-' == a[0]) {
			//以 - 开头、但不是单独一个 "-"（那是 stdin）的，都是不认识的选项。
			//以前这里会把 --nope 当成文件名，于是报成"文件打不开"（退出码 1），
			//把"用法错误"和"输入错误"混在一起了。
			cerr << "未知选项: " << a << endl;
			usage(argv[0]);
			return 2;
		} else if (NULL == path) {
			path = argv[i];
		} else {
			cerr << "多余的参数: " << a << endl;
			usage(argv[0]);
			return 2;
		}
	}

	//一个进程一份状态这件事已经结束了：库那边是对象，这里自己造一个。
	spclust::HAC hac(linkage);
	g_hac = &hac;

	//无参数或 "-" 都表示从标准输入读取。
	//边一边读一边就进结构了（不缓存整个输入），所以 init 阶段并进了这里。
	//返回 0 表示失败（原因 readdb 已经写到 stderr）。以前是 readdb 内部直接 exit(1)，
	//那样没法当库用；现在由这里给退出码。
	if (0 == readdb(path, &edge_to_library))
		return 1;

	//只有 --cut 比默认更低时才需要把树算得更深：没给 --floor 就自动跟着降，
	//这样"只想在更低的阈值处切一刀"不需要额外参数，也不会白白多算。
	//按簇数 / 簇大小切则必须算完整树（不知道要合并到多低），代价见 usage。
	if (!floor_given) {
		if (CUT_SIM != cut_mode || show_scipy)
			tree_floor = 0.0;   // 要正好 k 个簇 / 最小簇大小 / scipy 矩阵，都得算完整树
		else if (cut_sim < tree_floor)
			tree_floor = cut_sim;
	}

	//这里打印的是**生效的** floor（可能在上面刚被调低）
	if (DEBUG >= 1)
		cout << "Linkage=" << hac.linkage_name() << ", " << cut_desc()
			<< ", floor=" << tree_floor << endl << endl;

	//算合并树。计时写在 CLI 这边 —— 库里不打印任何东西。
	if (DEBUG >= 1)
		cout << "Clustering..." << endl;
	clock_t s_time = clock();
	hac.build(tree_floor);
	if (DEBUG >= 1)
		cout << "Done. (" << (double)(clock() - s_time) / CLOCKS_PER_SEC << "s)" << endl << endl;

	//打印扁平结果（把树切在 cut_sim 处）
	print_result();

	if (show_tree)
		print_tree();
	if (show_matrix)
		print_matrix();
	if (show_scipy)
		print_scipy();

	return 0;
}

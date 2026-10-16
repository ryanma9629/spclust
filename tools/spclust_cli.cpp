// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

#include "spclust/hac.h"

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

//切树的相似度阈值：默认取编译期常量 MIN_SIM_CLUST，可用 --cut 覆盖。
//
//合并树本身是完整算出来的（spclust::build），所以这个阈值只决定"从哪里切开"，
//换一个阈值不需要重新聚类。
double cut_sim = (double)MIN_SIM_CLUST;
//树算到多低。默认 = 切树阈值（于是成本与"直接在阈值处停下"完全一样）；
//用 --floor 显式指定（0 = 完整树）或者把 --cut 定得比默认更低时，会自动跟着降。
bool floor_given = false;
double tree_floor = (double)MIN_SIM_CLUST;
bool show_tree = false;
int tree_depth = 6;          // --tree 的深度上限；0 表示不限
bool show_matrix = false;

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
	spclust::add_edge(id1, id2, sim);
}

// ---------------------------------------------------------------------------
// 输出
// ---------------------------------------------------------------------------

//扁平分组（把树按 cut_sim 切开）
void print_result(void) {
	//spclust.cpp 的 DEBUG=2 日志把 cout 的精度调到了 12，这里显式改回 6 位有效数字
	cout << setprecision(6);

	Cut c = spclust::cut(cut_sim);

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

	cout << "Objects: " << spclust::num_objects() << ", groups: " << clusters.size() << endl;
	if (clusters.empty()) {
		cout << "No group reached 2 objects (all remaining similarities below cut="
			<< cut_sim << ")." << endl;
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
	const size_t n = spclust::num_objects();
	for (int i = 0; i < depth; i++)
		cout << "  ";
	if (node < n) {
		cout << "- " << names[node] << endl;
		return;
	}
	const Merge& r = spclust::tree()[node - n];
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
	const vector<Merge>& tr = spclust::tree();
	const size_t n = spclust::num_objects();
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
	const vector<Merge>& tr = spclust::tree();
	const size_t n = spclust::num_objects();
	cout << "# linkage matrix: left right height size" << endl;
	cout << "# 节点 id：0.." << (n ? n - 1 : 0) << " 是对象（按输入中首次出现的顺序），"
		<< n << "+i 是第 i 次合并。height 是**相似度**（越大越先合并），不是距离。" << endl;
	cout << setprecision(12);
	for (size_t i = 0; i < tr.size(); i++) {
		cout << tr[i].left << " " << tr[i].right << " "
			<< tr[i].similarity << " " << tr[i].size << endl;
	}
}

void usage(const char* prog) {
	cerr << "用法: " << prog << " [similarity.csv|-] [选项]" << endl
	     << "  similarity.csv   输入边表；省略或给 - 表示从标准输入读" << endl
	     << "  --cut <X>        把合并树切在相似度 X 处（默认 " << MIN_SIM_CLUST
	     << "，即编译期 MIN_SIM_CLUST）" << endl
	     << "  --floor <X>      树算到多低（默认 = cut；0 = 完整树，代价是慢 2~4 倍）" << endl
	     << "  --tree           打印树状图" << endl
	     << "  --tree-depth <N> 树状图的深度上限（默认 " << 6 << "，0 = 不限；隐含 --tree）" << endl
	     << "  --matrix         打印 linkage 矩阵（left right height size）" << endl
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

	//无参数或 "-" 都表示从标准输入读取。
	//边一边读一边就进结构了（不缓存整个输入），所以 init 阶段并进了这里。
	//返回 0 表示失败（原因 readdb 已经写到 stderr）。以前是 readdb 内部直接 exit(1)，
	//那样没法当库用；现在由这里给退出码。
	if (0 == readdb(path, &edge_to_library))
		return 1;

	if (DEBUG >= 1)
		cout << "Linkage=" << spclust::linkage_name() << ", cut=" << cut_sim << endl << endl;

	//只有 --cut 比默认更低时才需要把树算得更深：没给 --floor 就自动跟着降，
	//这样"只想在更低的阈值处切一刀"不需要额外参数，也不会白白多算。
	if (!floor_given && cut_sim < tree_floor)
		tree_floor = cut_sim;

	//算合并树。计时写在 CLI 这边 —— 库里不打印任何东西。
	if (DEBUG >= 1)
		cout << "Clustering..." << endl;
	clock_t s_time = clock();
	spclust::build(tree_floor);
	if (DEBUG >= 1)
		cout << "Done. (" << (double)(clock() - s_time) / CLOCKS_PER_SEC << "s)" << endl << endl;

	//打印扁平结果（把树切在 cut_sim 处）
	print_result();

	if (show_tree)
		print_tree();
	if (show_matrix)
		print_matrix();

	return 0;
}

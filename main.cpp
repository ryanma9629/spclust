#include "main.h"
#include "spclust.h"
#include "readdb.h"

//聚类结果（spclust.cpp 以 extern 方式使用）
vector<grp_s> grp;

namespace {

//切树的相似度阈值：默认取编译期常量 MIN_SIM_CLUST，可用 --cut 覆盖。
//
//合并树本身是完整算出来的（见 clust_grp），所以这个阈值只决定"从哪里切开"，
//换一个阈值不需要重新聚类。
double cut_sim = (double)MIN_SIM_CLUST;
//树算到多低。默认 = 切树阈值（于是成本与"直接在阈值处停下"完全一样）；
//用 --floor 显式指定（0 = 完整树）或者把 --cut 定得比默认更低时，会自动跟着降。
bool floor_given = false;
double tree_floor = (double)MIN_SIM_CLUST;
bool show_tree = false;
bool show_matrix = false;

//注意：这里用中性名字。原代码把 1 标为 complete、2 标为 single，
//但实现里 1 取的是 max、2 取的是 min(仅当两条边都存在)，
//按标准术语正好相反（max=单链接，min=全链接）。见 README。
const char* linkage_name(void) {
#if 0 == LINKAGE_MTD
	return "0=average";
#elif 1 == LINKAGE_MTD
	return "1=max";
#elif 2 == LINKAGE_MTD
	return "2=min(nonzero)";
#else
	return "unknown";
#endif
}

// ---------------------------------------------------------------------------
// 并查集（切树用）
// ---------------------------------------------------------------------------
vector<size_t> uf_parent;

size_t uf_find(size_t x) {
	while (uf_parent[x] != x) {
		uf_parent[x] = uf_parent[uf_parent[x]];   // 路径减半
		x = uf_parent[x];
	}
	return x;
}

void uf_union(size_t a, size_t b) {
	a = uf_find(a);
	b = uf_find(b);
	if (a != b)
		uf_parent[b] = a;
}

//把合并树在 cut 处切开，得到扁平分组。
//
//相似度沿合并序列单调不增（average / max / min 三种链接都是可归约的），
//所以"切"就是取 sim >= cut 的那一段前缀，再用并查集把每对兄弟并起来。
//剩下的每个根就是一个簇；簇的紧致度取该簇在前缀里最后一次合并记下的值。
struct cut_result {
	vector<vector<size_t>> members;   // 每个簇的对象 id
	vector<double> cohesion;          // 对应的组内平均相似度
};

cut_result cut_tree(double cut) {
	const size_t n = tree_object_count();
	const size_t total = n + merges.size();

	uf_parent.resize(total);
	for (size_t i = 0; i < total; i++)
		uf_parent[i] = i;

	size_t k = 0;
	while (k < merges.size() && merges[k].sim >= cut)
		k++;
	//注意：除了把两个子节点并起来，还必须把"这次合并产生的新节点"（id = n+i）也指到
	//同一个集合里 —— 否则后面某次合并拿 n+i 当子节点时，它会被当成一个独立元素，
	//于是只有"两个原始对象直接合并"那种才会生效。
	for (size_t i = 0; i < k; i++) {
		uf_union(merges[i].left, merges[i].right);
		uf_parent[n + i] = uf_find(merges[i].left);
	}

	//每个簇的紧致度 = 该簇在前缀里最后一次合并写下的 cohesion（后写覆盖先写）
	vector<double> root_cohesion(total, 0.0);
	for (size_t i = 0; i < k; i++)
		root_cohesion[uf_find(merges[i].left)] = merges[i].cohesion;

	//按根把对象归拢
	const size_t NONE = (size_t)-1;
	vector<size_t> root_index(total, NONE);
	cut_result res;
	for (size_t id = 0; id < n; id++) {
		size_t r = uf_find(id);
		if (NONE == root_index[r]) {
			root_index[r] = res.members.size();
			res.members.push_back(vector<size_t>());
			res.cohesion.push_back(root_cohesion[r]);
		}
		res.members[root_index[r]].push_back(id);
	}
	return res;
}

// ---------------------------------------------------------------------------
// 输出
// ---------------------------------------------------------------------------

//扁平分组（把树按 cut_sim 切开）
void print_result(void) {
	//spclust.cpp 的 DEBUG=2 日志把 cout 的精度调到了 12，这里显式改回 6 位有效数字
	cout << setprecision(6);

	cut_result c = cut_tree(cut_sim);

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

	cout << "Objects: " << tree_object_count() << ", groups: " << clusters.size() << endl;
	if (clusters.empty()) {
		cout << "No group reached 2 objects (all remaining similarities below cut="
			<< cut_sim << ")." << endl;
	}
	for (size_t k = 0; k < clusters.size(); k++) {
		const vector<size_t>& m = c.members[clusters[k]];
		cout << "Group " << k << " (size=" << m.size() << ", sim=" << c.cohesion[clusters[k]] << "):";
		for (size_t j = 0; j < m.size(); j++)
			cout << " " << object_name(m[j]);
		cout << endl;
	}

	//未成簇（仍是单元素）的对象
	vector<string> single;
	for (size_t i = 0; i < c.members.size(); i++) {
		if (c.members[i].size() == 1)
			single.push_back(object_name(c.members[i][0]));
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
	const size_t n = tree_object_count();
	for (int i = 0; i < depth; i++)
		cout << "  ";
	if (node < n) {
		cout << "- " << object_name(node) << endl;
		return;
	}
	const merge_rec& r = merges[node - n];
	cout << "+ merged at sim=" << r.sim << " (size=" << r.size << ")" << endl;
	if (depth + 1 > max_depth) {
		for (int i = 0; i <= depth; i++)
			cout << "  ";
		cout << "... (" << r.size << " objects, 用 --matrix 取完整数据)" << endl;
		return;
	}
	print_subtree(r.left, depth + 1, max_depth);
	print_subtree(r.right, depth + 1, max_depth);
}

void print_tree(void) {
	const size_t n = tree_object_count();
	//根 = 没有当过任何一次合并的子节点的节点
	vector<char> is_child(n + merges.size(), 0);
	for (size_t i = 0; i < merges.size(); i++) {
		is_child[merges[i].left] = 1;
		is_child[merges[i].right] = 1;
	}
	cout << "Dendrogram (" << n << " objects, " << merges.size() << " merges; "
		<< "阈值只是切树位置，树本身是完整的):" << endl;
	const int max_depth = 6;
	size_t roots = 0;
	for (size_t node = 0; node < n + merges.size(); node++) {
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
	const size_t n = tree_object_count();
	cout << "# linkage matrix: left right height size" << endl;
	cout << "# 节点 id：0.." << (n ? n - 1 : 0) << " 是对象（按输入中首次出现的顺序），"
		<< n << "+i 是第 i 次合并。height 是**相似度**（越大越先合并），不是距离。" << endl;
	cout << setprecision(12);
	for (size_t i = 0; i < merges.size(); i++) {
		cout << merges[i].left << " " << merges[i].right << " "
			<< merges[i].sim << " " << merges[i].size << endl;
	}
}

void usage(const char* prog) {
	cerr << "用法: " << prog << " [similarity.csv|-] [选项]" << endl
	     << "  similarity.csv   输入边表；省略或给 - 表示从标准输入读" << endl
	     << "  --cut <X>        把合并树切在相似度 X 处（默认 " << MIN_SIM_CLUST
	     << "，即编译期 MIN_SIM_CLUST）" << endl
	     << "  --floor <X>      树算到多低（默认 = cut；0 = 完整树，代价是慢 2~4 倍）" << endl
	     << "  --tree           打印树状图" << endl
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
	readdb(path, &init_add_edge);

	if (DEBUG >= 1)
		cout << "Linkage=" << linkage_name() << ", cut=" << cut_sim << endl << endl;

	//只有 --cut 比默认更低时才需要把树算得更深：没给 --floor 就自动跟着降，
	//这样"只想在更低的阈值处切一刀"不需要额外参数，也不会白白多算。
	if (!floor_given && cut_sim < tree_floor)
		tree_floor = cut_sim;

	//聚类（并记下合并树）
	clust_grp(tree_floor);

	//打印扁平结果（把树切在 cut_sim 处）
	print_result();

	if (show_tree)
		print_tree();
	if (show_matrix)
		print_matrix();

	return 0;
}

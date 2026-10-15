#include "main.h"
#include "spclust.h"
#include "readdb.h"

//聚类结果（spclust.cpp 以 extern 方式使用）
vector<grp_s> grp;

namespace {

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

//按“簇大小降序、组内相似度降序”排序
bool grp_greater(const size_t& a, const size_t& b) {
	if (grp[a].objs.size() != grp[b].objs.size())
		return grp[a].objs.size() > grp[b].objs.size();
	return grp[a].sim > grp[b].sim;
}

void print_result(void) {
	vector<size_t> clusters;
	size_t n_objects = 0;
	for (size_t i = 0; i < grp.size(); i++) {
		n_objects += grp[i].objs.size();
		if (grp[i].objs.size() >= 2)
			clusters.push_back(i);
	}
	sort(clusters.begin(), clusters.end(), grp_greater);

	cout << "Objects: " << n_objects << ", groups: " << clusters.size() << endl;
	if (clusters.empty()) {
		cout << "No group reached 2 objects (all remaining similarities below min_sim="
			<< MIN_SIM_CLUST << ")." << endl;
	}
	for (size_t k = 0; k < clusters.size(); k++) {
		const grp_s& g = grp[clusters[k]];
		cout << "Group " << k << " (size=" << g.objs.size() << ", sim=" << g.sim << "):";
		for (size_t j = 0; j < g.objs.size(); j++)
			cout << " " << g.objs[j];
		cout << endl;
	}

	//未成簇（仍是单元素组）的对象
	vector<string> single;
	for (size_t i = 0; i < grp.size(); i++) {
		if (grp[i].objs.size() == 1)
			single.push_back(grp[i].objs[0]);
	}
	sort(single.begin(), single.end());
	if (!single.empty()) {
		cout << "Unclustered (" << single.size() << "):";
		for (size_t i = 0; i < single.size(); i++)
			cout << " " << single[i];
		cout << endl;
	}
}

} // namespace

int main(int argc, char* argv[]) {
	if (argc > 2) {
		cerr << "Usage: " << argv[0] << " [similarity.csv|-]" << endl;
		return 2;
	}

	//无参数或 "-" 都表示从标准输入读取
	readdb(argc > 1 ? argv[1] : NULL);

	if (DEBUG >= 1) {
		cout << "Linkage=" << linkage_name() << ", min_sim=" << MIN_SIM_CLUST << endl << endl;
	}

	//获得聚类结果
	clust_grp();

	//打印聚类结果
	print_result();

	return 0;
}

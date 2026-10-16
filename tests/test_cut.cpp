// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
//
// 库级属性测试：cut / cut_n_clusters / cut_min_size 三个切树接口。
//
// 差分测试（tests/run_tests.py）从 CLI 侧验证"合并序列"对不对，这个文件从库侧
// 验证"切法"对不对。期望值不是用同样的并查集再算一遍，而是**从合并树自底向上
// 收集叶子集合**（另一条代码路径）：树里每个节点天然对应一个叶子集合，于是
//   * 切在前 i 次合并之后 -> 该前缀里"没有父节点"的那些节点的叶子集合
//   * 正好 k 个簇          -> 切在 n-k 次之后
//   * 没有小于 m 的簇      -> 满足条件的最浅切法
// 三种期望值都能直接算出来，再和库里三个接口的返回逐一对比。
//
// 编译：见 Makefile 的 `make test-cut` / CMakeLists.txt 的 cut_property 测试。

#include "spclust/hac.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace {

using spclust::Cut;
using spclust::HAC;
using spclust::Linkage;
using spclust::Merge;

int g_pass = 0;
int g_fail = 0;
const char* g_ctx = "";

void check(bool ok, const std::string& what) {
	if (ok) {
		g_pass++;
	} else {
		g_fail++;
		std::printf("  ✗ [%s] %s\n", g_ctx, what.c_str());
	}
}

// 合并树自底向上算每个节点的叶子集合（独立于库内部那套并查集）。
struct TreeView {
	std::vector<std::vector<std::uint32_t> > leaves;   // 节点 id -> 叶子（升序）
	std::vector<std::size_t> size;
	std::vector<double> cohesion;

	explicit TreeView(std::size_t n, const std::vector<Merge>& tr) {
		leaves.resize(n + tr.size());
		size.assign(n + tr.size(), 0);
		cohesion.assign(n + tr.size(), 0.0);
		for (std::size_t i = 0; i < n; i++) {
			leaves[i].push_back((std::uint32_t)i);
			size[i] = 1;
		}
		for (std::size_t i = 0; i < tr.size(); i++) {
			std::size_t node = n + i;
			const Merge& m = tr[i];
			leaves[node] = leaves[m.left];
			leaves[node].insert(leaves[node].end(), leaves[m.right].begin(), leaves[m.right].end());
			std::sort(leaves[node].begin(), leaves[node].end());
			size[node] = size[m.left] + size[m.right];
			cohesion[node] = m.cohesion;
		}
	}

	// 树的前 use 次合并之后，"没有父节点"的节点就是各个簇。
	std::vector<std::size_t> roots_after(std::size_t n, const std::vector<Merge>& tr,
										 std::size_t use) {
		if (use > tr.size()) use = tr.size();      // 树可能是残缺的（见下）
		std::vector<char> has_parent(n + tr.size(), 0);
		for (std::size_t i = 0; i < use; i++) {
			has_parent[tr[i].left] = 1;
			has_parent[tr[i].right] = 1;
		}
		std::vector<std::size_t> out;
		const std::size_t limit = n + use;
		for (std::size_t node = 0; node < limit; node++)
			if (!has_parent[node])
				out.push_back(node);
		return out;
	}
};

// 把 Cut 规范化，便于比较：簇内升序、簇之间按"最小成员"升序。
std::vector<std::vector<std::uint32_t> > canon(const Cut& c) {
	std::vector<std::vector<std::uint32_t> > v = c.members;
	for (std::size_t i = 0; i < v.size(); i++)
		std::sort(v[i].begin(), v[i].end());
	std::sort(v.begin(), v.end());
	return v;
}

std::vector<std::vector<std::uint32_t> > canon_nodes(
	const TreeView& tv, const std::vector<std::size_t>& roots) {
	std::vector<std::vector<std::uint32_t> > v;
	for (std::size_t i = 0; i < roots.size(); i++)
		v.push_back(tv.leaves[roots[i]]);
	std::sort(v.begin(), v.end());
	return v;
}

std::string show(const std::vector<std::vector<std::uint32_t> >& v) {
	std::string s = "{";
	for (std::size_t i = 0; i < v.size(); i++) {
		s += "{";
		for (std::size_t j = 0; j < v[i].size(); j++) {
			if (j) s += ",";
			s += std::to_string(v[i][j]);
		}
		s += "}";
	}
	return s + "}";
}

// 一个连通图的随机用例：先随机生成一棵生成树（保证连通），再加若干随机边。
struct Case {
	std::vector<std::uint32_t> a, b;
	std::vector<double> sim;
	std::size_t n;
};

Case make_case(std::size_t n, std::size_t extra, std::mt19937& rng) {
	Case c;
	c.n = n;
	std::uniform_real_distribution<double> sim(1.0, 100.0);
	std::vector<std::size_t> order(n);
	for (std::size_t i = 0; i < n; i++) order[i] = i;
	std::shuffle(order.begin(), order.end(), rng);
	for (std::size_t i = 1; i < n; i++) {
		std::uniform_int_distribution<std::size_t> pick(0, i - 1);
		c.a.push_back((std::uint32_t)order[i]);
		c.b.push_back((std::uint32_t)order[pick(rng)]);
		c.sim.push_back(sim(rng));
	}
	std::uniform_int_distribution<std::size_t> pick(0, n - 1);
	for (std::size_t i = 0; i < extra; i++) {
		std::size_t x = pick(rng), y = pick(rng);
		if (x == y) continue;
		c.a.push_back((std::uint32_t)x);
		c.b.push_back((std::uint32_t)y);
		c.sim.push_back(sim(rng));
	}
	return c;
}

void run_one(const Case& c, Linkage mtd, const char* mtd_name, unsigned seed) {
	HAC hac(mtd);
	for (std::size_t i = 0; i < c.a.size(); i++)
		hac.add_edge(c.a[i], c.b[i], c.sim[i]);
	hac.build(0.0);                     // 完整树：这样三种切法都能到最深

	const std::vector<Merge>& tr = hac.tree();
	const std::size_t n = hac.num_objects();
	TreeView tv(n, tr);

	g_ctx = mtd_name;
	char buf[128];
	std::snprintf(buf, sizeof(buf), "%s n=%zu m=%zu seed=%u", mtd_name, n, tr.size(), seed);
	g_ctx = buf;

	// 树自洽：每个节点的叶子数 == 记录的大小 == 左右子节点大小之和
	for (std::size_t i = 0; i < tr.size(); i++) {
		std::size_t node = n + i;
		check(tv.size[node] == tr[i].size, "树里记的 size 与叶子数一致");
		check(tv.size[node] == tv.size[tr[i].left] + tv.size[tr[i].right], "size 可加");
		check(tr[i].left < node && tr[i].right < node, "子节点 id 小于父节点");
	}

	// ① cut(相似度)："切在 t"就是取所有 similarity >= t 的合并（并列会一起进来）。
	// 期望值在测试里独立数一遍前缀长度，再用叶子集合路线拼出分区。
	for (std::size_t i = 0; i < tr.size(); i++) {
		const double t = tr[i].similarity;
		std::size_t use = 0;
		while (use < tr.size() && tr[use].similarity >= t)
			use++;                                  // 并列整段
		Cut got = hac.cut(t);
		std::vector<std::vector<std::uint32_t> > want =
			canon_nodes(tv, tv.roots_after(n, tr, use));
		if (canon(got) != want)
			check(false, "cut(t_" + std::to_string(i) + "=" + std::to_string(t) + "): " +
				 show(canon(got)) + " vs " + show(want));
		else
			g_pass++;
	}
	// 阈值高于最大值 -> 一次都不合并；阈值 <= 0 -> 把树能给的都合并
	check(canon(hac.cut(tr.empty() ? 1.0 : tr[0].similarity + 1.0)) ==
		  canon_nodes(tv, tv.roots_after(n, tr, 0)), "阈值高于最大相似度 -> 全部分散");
	check(canon(hac.cut(0.0)) == canon_nodes(tv, tv.roots_after(n, tr, tr.size())),
		  "cut(0) -> 树能给的都合并了");

	// ② cut_n_clusters(k)：应当等于切在前 n-k 次之后；簇数正好 k（连通图 + 完整树）
	for (std::size_t k = 1; k <= n; k++) {
		Cut got = hac.cut_n_clusters(k);
		const std::size_t use = (n >= k) ? (n - k) : 0;
		std::vector<std::vector<std::uint32_t> > want =
			canon_nodes(tv, tv.roots_after(n, tr, use));
		// "正好 k 个簇"只有在树够深时才成立：min-nonzero 下，两簇之间根本没有边
		// 就不是候选，于是可能留下互不相连的分量（残缺树），这时只能给最细的切法。
		bool deep_enough = (use <= tr.size());
		bool ok = (canon(got) == want) && (!deep_enough || got.members.size() == k);
		if (!ok) {
			check(false, "cut_n_clusters(" + std::to_string(k) + "): " + show(canon(got)) +
				 " (" + std::to_string(got.members.size()) + " 个) vs " + show(want) +
				 " (" + std::to_string(want.size()) + " 个)");
			break;
		}
		g_pass++;
	}
	// 0 按 1 处理：与 k=1 完全等价（残缺树下"1"也只能给最细的切法）
	check(canon(hac.cut_n_clusters(0)) == canon(hac.cut_n_clusters(1)),
		  "cut_n_clusters(0) 按 1 处理");
	// 树完整（全部对象最终并成一个簇）时，k=1 必须正好是一个簇
	{
		std::vector<std::size_t> full = tv.roots_after(n, tr, tr.size());
		if (full.size() == 1 && tv.size[full[0]] == n)
			check(hac.cut_n_clusters(1).members.size() == 1, "完整树下 cut_n_clusters(1) 是一个簇");
	}
	check(canon(hac.cut_n_clusters(n + 5)) == canon(hac.cut_n_clusters(n)),
		  "k > n 时与 k == n 相同（全部分散）");

	// ③ cut_min_size(m)：应当是"所有簇都 >= m"的最浅切法
	for (std::size_t m = 1; m <= n + 1; m++) {
		std::size_t expect_use = tr.size();
		for (std::size_t use = 0; use <= tr.size(); use++) {
			std::vector<std::size_t> roots = tv.roots_after(n, tr, use);
			bool all_big = true;
			for (std::size_t i = 0; i < roots.size(); i++)
				if (tv.size[roots[i]] < m) { all_big = false; break; }
			if (all_big) { expect_use = use; break; }
		}
		Cut got = hac.cut_min_size(m);
		std::vector<std::vector<std::uint32_t> > want =
			canon_nodes(tv, tv.roots_after(n, tr, expect_use));
		if (canon(got) != want) {
			check(false, "cut_min_size(" + std::to_string(m) + "): " + show(canon(got)) +
				 " vs 期望 " + show(want));
			break;
		}
		// 只有在树里真的存在"所有簇都 >= m"的前缀时才断言；残缺树（无边的分量永远
		// 合不起来）达不到这个条件，这时返回的是最细切法。
		if (expect_use < tr.size() || tr.size() + 1 == n + (n - 1))
			for (std::size_t i = 0; i < got.members.size(); i++)
				check(got.members[i].size() >= m, "cut_min_size 的结果里没有小于 m 的簇");
		g_pass++;
	}

	// ④ 每个簇的紧致度 == 树里"造出这个簇的那次合并"记下的 cohesion
	for (std::size_t k = 1; k <= n; k++) {
		Cut got = hac.cut_n_clusters(k);
		std::vector<std::size_t> roots = tv.roots_after(n, tr, n - k);
		std::map<std::vector<std::uint32_t>, std::size_t> node_of;
		for (std::size_t i = 0; i < roots.size(); i++) node_of[tv.leaves[roots[i]]] = roots[i];
		for (std::size_t i = 0; i < got.members.size(); i++) {
			std::vector<std::uint32_t> mem = got.members[i];
			std::sort(mem.begin(), mem.end());
			std::size_t node = node_of[mem];
			check(got.cohesion[i] == tv.cohesion[node], "簇的紧致度来自造出它的那次合并");
		}
	}

	// ⑤ 单调性：k 越小应当越粗（每个簇都是更细切法里若干簇的并）
	for (std::size_t k = 2; k <= n; k++) {
		Cut fine = hac.cut_n_clusters(k), coarse = hac.cut_n_clusters(k - 1);
		std::vector<std::size_t> owner(n, (std::size_t)-1);
		for (std::size_t i = 0; i < fine.members.size(); i++)
			for (std::size_t j = 0; j < fine.members[i].size(); j++)
				owner[fine.members[i][j]] = i;
		bool ok = true;
		for (std::size_t i = 0; i < coarse.members.size() && ok; i++) {
			std::set<std::size_t> seen;
			for (std::size_t j = 0; j < coarse.members[i].size(); j++)
				seen.insert(owner[coarse.members[i][j]]);
			if (seen.size() != 1 && coarse.members[i].size() != 0) {
				// 一个粗簇可以由多个细簇组成，但必须**整块**是细簇的并：
				// 即每个细簇要么整个在里面，要么整个在外面
				for (std::size_t q = 0; q < fine.members.size() && ok; q++) {
					std::size_t in = 0;
					for (std::size_t j = 0; j < fine.members[q].size(); j++)
						for (std::size_t z = 0; z < coarse.members[i].size(); z++)
							if (coarse.members[i][z] == fine.members[q][j]) in++;
					if (in != 0 && in != fine.members[q].size()) ok = false;
				}
			}
		}
		check(ok, "cut_n_clusters(k-1) 的每个簇都是 cut_n_clusters(k) 里若干个簇的并");
	}
}

// 把一棵树拍成可比较的字符串
std::string tree_key(const std::vector<Merge>& tr) {
	std::string s;
	char buf[128];
	for (std::size_t i = 0; i < tr.size(); i++) {
		std::snprintf(buf, sizeof(buf), "%u,%u,%.17g,%zu;", tr[i].left, tr[i].right,
					  tr[i].similarity, tr[i].size);
		s += buf;
	}
	return s;
}

void feed(HAC& hac, const Case& c) {
	for (std::size_t i = 0; i < c.a.size(); i++)
		hac.add_edge(c.a[i], c.b[i], c.sim[i]);
}

// ⑥⑦ 多实例 / 可重入：这是把 34 个文件静态变量搬进对象之后最关键的一条性质 ——
// 以前一个进程只能跑一次聚类，现在两组状态交错使用必须互不影响。
void test_instances(const Case& c1, const Case& c2) {
	g_ctx = "多实例";
	// 单独跑一遍作为参照
	std::string solo1, solo1min, solo2;
	{
		HAC a(Linkage::Average);
		feed(a, c1);
		a.build(0.0);
		solo1 = tree_key(a.tree());
		HAC m(Linkage::MinNonzero);
		feed(m, c1);
		m.build(0.0);
		solo1min = tree_key(m.tree());
		HAC b(Linkage::Average);
		feed(b, c2);
		b.build(0.0);
		solo2 = tree_key(b.tree());
	}
	// 交错使用：边交错喂、树交错算
	{
		HAC a(Linkage::Average), m(Linkage::MinNonzero), b(Linkage::Average);
		std::size_t n = std::max(c1.a.size(), c2.a.size());
		for (std::size_t i = 0; i < n; i++) {
			if (i < c1.a.size()) {
				a.add_edge(c1.a[i], c1.b[i], c1.sim[i]);
				m.add_edge(c1.a[i], c1.b[i], c1.sim[i]);
			}
			if (i < c2.a.size())
				b.add_edge(c2.a[i], c2.b[i], c2.sim[i]);
		}
		a.build(0.0);
		m.build(0.0);
		b.build(0.0);
		check(tree_key(a.tree()) == solo1, "交错喂边后 average 的树不变");
		check(tree_key(m.tree()) == solo1min, "交错喂边后 min-nonzero 的树不变");
		check(tree_key(b.tree()) == solo2, "交错喂边后第二个数据集的树不变");
	}
	// reset 之后复用同一个对象，结果必须和新对象一样
	{
		HAC r(Linkage::Average);
		feed(r, c1);
		r.build(0.0);
		(void)r.cut(50.0);
		r.reset();
		check(r.num_objects() == 0 && r.tree().empty(), "reset 之后状态清空");
		feed(r, c2);
		r.build(0.0);
		check(tree_key(r.tree()) == solo2, "reset 复用同一个对象，结果与新对象一致");
	}
	// 对象之间互不影响：销毁其中一个，另一个的结果不变
	{
		HAC keep(Linkage::Average);
		feed(keep, c2);
		{
			HAC tmp(Linkage::Average);
			feed(tmp, c1);
			tmp.build(0.0);
		}
		keep.build(0.0);
		check(tree_key(keep.tree()) == solo2, "另一个对象析构后结果不变");
	}
}

} // namespace

int main() {
	const Linkage mtds[3] = {Linkage::Average, Linkage::MaxSimilarity, Linkage::MinNonzero};
	const char* names[3] = {"average", "max-similarity", "min-nonzero"};

	std::mt19937 rng(20261016);
	unsigned seed = 0;
	for (std::size_t n : {std::size_t(2), std::size_t(3), std::size_t(5), std::size_t(8),
						  std::size_t(13), std::size_t(21), std::size_t(40)}) {
		for (std::size_t extra : {std::size_t(0), std::size_t(3), std::size_t(25)}) {
			Case c = make_case(n, extra, rng);
			for (int m = 0; m < 3; m++)
				run_one(c, mtds[m], names[m], seed);
			seed++;
		}
	}

	// 多实例 / 可重入（用两组随机用例，交错喂边）
	{
		std::mt19937 rng2(777);
		Case c1 = make_case(11, 20, rng2);
		Case c2 = make_case(17, 5, rng2);
		test_instances(c1, c2);
	}

	std::printf("切树属性测试: PASS=%d FAIL=%d\n", g_pass, g_fail);
	return g_fail == 0 ? 0 : 1;
}

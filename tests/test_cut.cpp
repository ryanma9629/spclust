// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
//
// Library-level property test: the three cut entry points.
//
// The differential test (tests/run_tests.py) checks the merge sequence from the CLI side; this
// file checks the cut itself from the library side. The expectations are not a second copy of the
// union-find but come from **collecting leaf sets bottom-up from the merge tree** - a different code
// path. Every node of the tree is a leaf set, so:
//   * cut after i merges -> the leaf sets of the nodes with no parent in that prefix
//   * exactly k clusters -> cut after n-k merges
//   * no cluster below m -> the shallowest cut satisfying that
//
// Build: see `make test-cut` in the Makefile / the cut_property test in CMakeLists.txt.

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

// Compute every node's leaf set bottom-up (independent of the union-find inside the library).
struct TreeView {
	std::vector<std::vector<std::uint32_t> > leaves;   // node id -> leaves (ascending)
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

	// After the first `use` merges, the nodes with no parent are the clusters.
	std::vector<std::size_t> roots_after(std::size_t n, const std::vector<Merge>& tr,
										 std::size_t use) {
		if (use > tr.size()) use = tr.size();      // the tree may be incomplete (see below)
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

// Canonicalize a Cut for comparison: members ascending within a cluster, clusters sorted.
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

// A random connected case: generate a random spanning tree first (guaranteeing connectivity),
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
	hac.build(0.0);                     // full tree: so all three cut modes can reach the deepest level

	const std::vector<Merge>& tr = hac.tree();
	const std::size_t n = hac.num_objects();
	TreeView tv(n, tr);

	g_ctx = mtd_name;
	char buf[128];
	std::snprintf(buf, sizeof(buf), "%s n=%zu m=%zu seed=%u", mtd_name, n, tr.size(), seed);
	g_ctx = buf;

	// Tree self-consistency: a node's leaf count == its recorded size == the sum of its children
	for (std::size_t i = 0; i < tr.size(); i++) {
		std::size_t node = n + i;
		check(tv.size[node] == tr[i].size, "the recorded size matches the leaf count");
		check(tv.size[node] == tv.size[tr[i].left] + tv.size[tr[i].right], "sizes are additive");
		check(tr[i].left < node && tr[i].right < node, "child ids are below the parent id");
	}

	// (1) cut(similarity): "cut at t" takes every merge with similarity >= t (ties included).
	// The expected value counts the prefix length independently here, then builds the partition
	for (std::size_t i = 0; i < tr.size(); i++) {
		const double t = tr[i].similarity;
		std::size_t use = 0;
		while (use < tr.size() && tr[use].similarity >= t)
			use++;                                  // the whole tie block
		Cut got = hac.cut(t);
		std::vector<std::vector<std::uint32_t> > want =
			canon_nodes(tv, tv.roots_after(n, tr, use));
		if (canon(got) != want)
			check(false, "cut(t_" + std::to_string(i) + "=" + std::to_string(t) + "): " +
				 show(canon(got)) + " vs " + show(want));
		else
			g_pass++;
	}
	// A threshold above the maximum merges nothing; a threshold <= 0 merges everything the tree has
	check(canon(hac.cut(tr.empty() ? 1.0 : tr[0].similarity + 1.0)) ==
		  canon_nodes(tv, tv.roots_after(n, tr, 0)), "a threshold above the maximum similarity leaves everything separate");
	check(canon(hac.cut(0.0)) == canon_nodes(tv, tv.roots_after(n, tr, tr.size())),
		  "cut(0) merges everything the tree offers");

	// (2) cut_n_clusters(k): equal to cutting after n-k merges, with exactly k clusters
	for (std::size_t k = 1; k <= n; k++) {
		Cut got = hac.cut_n_clusters(k);
		const std::size_t use = (n >= k) ? (n - k) : 0;
		std::vector<std::vector<std::uint32_t> > want =
			canon_nodes(tv, tv.roots_after(n, tr, use));
		// "Exactly k clusters" only holds when the tree is deep enough: under MinNonzero two
		// clusters with no edge are not even candidates, so disconnected components may remain and
		bool deep_enough = (use <= tr.size());
		bool ok = (canon(got) == want) && (!deep_enough || got.members.size() == k);
		if (!ok) {
			check(false, "cut_n_clusters(" + std::to_string(k) + "): " + show(canon(got)) +
				 " (" + std::to_string(got.members.size()) + " clusters) vs " + show(want) +
				 " (" + std::to_string(want.size()) + " clusters)");
			break;
		}
		g_pass++;
	}
	// 0 is treated as 1: fully equivalent to k=1 (under an incomplete tree "1" also gives the
	check(canon(hac.cut_n_clusters(0)) == canon(hac.cut_n_clusters(1)),
		  "cut_n_clusters(0) is treated as 1");
	// With a complete tree (everything ends in one cluster), k=1 must be exactly one cluster
	{
		std::vector<std::size_t> full = tv.roots_after(n, tr, tr.size());
		if (full.size() == 1 && tv.size[full[0]] == n)
			check(hac.cut_n_clusters(1).members.size() == 1, "under a complete tree cut_n_clusters(1) is one cluster");
	}
	check(canon(hac.cut_n_clusters(n + 5)) == canon(hac.cut_n_clusters(n)),
		  "k > n behaves like k == n (everything separate)");

	// (3) cut_min_size(m): the shallowest cut where every cluster is >= m
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
				 " vs expected " + show(want));
			break;
		}
		// Only assert when such a prefix really exists; an incomplete tree (components with no edges
		// never merge) cannot reach it and returns the finest cut instead.
		if (expect_use < tr.size() || tr.size() + 1 == n + (n - 1))
			for (std::size_t i = 0; i < got.members.size(); i++)
				check(got.members[i].size() >= m, "cut_min_size produced no cluster below m");
		g_pass++;
	}

	// (4) Every cluster's cohesion equals the value recorded by the merge that created it
	for (std::size_t k = 1; k <= n; k++) {
		Cut got = hac.cut_n_clusters(k);
		std::vector<std::size_t> roots = tv.roots_after(n, tr, n - k);
		std::map<std::vector<std::uint32_t>, std::size_t> node_of;
		for (std::size_t i = 0; i < roots.size(); i++) node_of[tv.leaves[roots[i]]] = roots[i];
		for (std::size_t i = 0; i < got.members.size(); i++) {
			std::vector<std::uint32_t> mem = got.members[i];
			std::sort(mem.begin(), mem.end());
			std::size_t node = node_of[mem];
			check(got.cohesion[i] == tv.cohesion[node], "cluster cohesion comes from the merge that created it");
		}
	}

	// (5) Monotonicity: a smaller k must be coarser (each cluster is a union of finer ones)
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
				// A coarse cluster may consist of several fine ones, but each fine cluster must be
				// entirely inside or entirely outside it
				for (std::size_t q = 0; q < fine.members.size() && ok; q++) {
					std::size_t in = 0;
					for (std::size_t j = 0; j < fine.members[q].size(); j++)
						for (std::size_t z = 0; z < coarse.members[i].size(); z++)
							if (coarse.members[i][z] == fine.members[q][j]) in++;
					if (in != 0 && in != fine.members[q].size()) ok = false;
				}
			}
		}
		check(ok, "every cluster of cut_n_clusters(k-1) is a union of clusters from cut_n_clusters(k)");
	}
}

// Flatten a tree into a comparable string
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

// (6)(7) Multiple instances / reentrancy: the property that mattered most after moving 34 file
// statics into an object. One clustering per process used to be the limit; now two states used
void test_instances(const Case& c1, const Case& c2) {
	g_ctx = "multiple instances";
	// Run each on its own first, as a reference
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
	// Interleaved use: interleaved edges, interleaved builds
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
		check(tree_key(a.tree()) == solo1, "the average tree is unchanged after interleaved edges");
		check(tree_key(m.tree()) == solo1min, "the min-nonzero tree is unchanged after interleaved edges");
		check(tree_key(b.tree()) == solo2, "the second dataset's tree is unchanged after interleaved edges");
	}
	// Reusing one object after reset() must give the same result as a fresh one
	{
		HAC r(Linkage::Average);
		feed(r, c1);
		r.build(0.0);
		(void)r.cut(50.0);
		r.reset();
		check(r.num_objects() == 0 && r.tree().empty(), "reset clears the state");
		feed(r, c2);
		r.build(0.0);
		check(tree_key(r.tree()) == solo2, "after reset the reused object matches a fresh one");
	}
	// Objects do not affect each other: destroying one leaves the other unchanged
	{
		HAC keep(Linkage::Average);
		feed(keep, c2);
		{
			HAC tmp(Linkage::Average);
			feed(tmp, c1);
			tmp.build(0.0);
		}
		keep.build(0.0);
		check(tree_key(keep.tree()) == solo2, "the other object is unchanged after one is destroyed");
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

	// Multiple instances / reentrancy (two random cases, interleaved edges)
	{
		std::mt19937 rng2(777);
		Case c1 = make_case(11, 20, rng2);
		Case c2 = make_case(17, 5, rng2);
		test_instances(c1, c2);
	}

	std::printf("cut property test: PASS=%d FAIL=%d\n", g_pass, g_fail);
	return g_fail == 0 ? 0 : 1;
}

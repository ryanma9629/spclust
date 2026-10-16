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

//This file is the library's first user: it owns **all** IO and presentation -
//reading files, assigning ids to objects, mapping ids back to names, timing,
//printing. The library has no strings, no standard output and no global name table.

namespace {

//This program runs a single clustering, so a global pointer to it is enough for
//the print functions below. (The library holds one state per object and has no
//globals; this is just the CLI taking a shortcut for itself.)
spclust::HAC* g_hac = NULL;


//The cut mode. The three are mutually exclusive; the last one given on the
//command line wins:
//   --cut <X>        cut at similarity X (default; X is the compile-time MIN_SIM_CLUST)
//   --n-clusters <K> exactly K clusters
//   --min-size <M>   merge until no cluster is smaller than M
//
//The merge tree is computed in full (spclust::build), so this threshold only
//decides where to cut; changing it does not require re-clustering.
enum CutMode { CUT_SIM, CUT_N, CUT_SIZE };
CutMode cut_mode = CUT_SIM;
double cut_sim = (double)MIN_SIM_CLUST;
size_t cut_n_clusters = 0;
size_t cut_min_size = 0;

//Cut the tree according to the current mode
Cut do_cut(void) {
	switch (cut_mode) {
	case CUT_N:    return g_hac->cut_n_clusters(cut_n_clusters);
	case CUT_SIZE: return g_hac->cut_min_size(cut_min_size);
	default:       return g_hac->cut(cut_sim);
	}
}

//One-line description of the current cut mode (for messages)
string cut_desc(void) {
	char buf[128];
	switch (cut_mode) {
	case CUT_N:    snprintf(buf, sizeof(buf), "n-clusters=%zu", cut_n_clusters); break;
	case CUT_SIZE: snprintf(buf, sizeof(buf), "min-size=%zu", cut_min_size); break;
	default:       snprintf(buf, sizeof(buf), "cut=%g", cut_sim); break;
	}
	return string(buf);
}
//How far down to compute the tree. Default = the cut threshold (so the cost is
//identical to stopping right at that threshold). --floor overrides it (0 = full
//tree), and lowering --cut below the default lowers it automatically.
bool floor_given = false;
double tree_floor = (double)MIN_SIM_CLUST;
bool show_tree = false;
int tree_depth = 6;          // depth limit for --tree; 0 means unlimited
bool show_matrix = false;
bool show_scipy = false;                 // --scipy: export a scipy linkage matrix
spclust::DistanceMap scipy_map = spclust::DistanceMap::Complement;
bool scipy_bridge = false;               // bridge a forest into one tree (scipy wants n-1 rows)

//The linkage method. It used to be a compile-time constant (LINKAGE_MTD) and
//changing it meant recompiling; in the library it is a constructor argument, so
//one binary covers all three and --linkage selects it at run time.
spclust::Linkage linkage = spclust::Linkage::Average;

//object id -> name. The library only knows ids; the name mapping is the CLI's job.
vector<string> names;
unordered_map<string, uint32_t> name_to_id;

//Hand one edge to the library. Ids are assigned in order of first appearance,
//which matches the node-id convention, so a leaf id in the tree is an index here.
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
// Output
// ---------------------------------------------------------------------------

//Flat partition (the tree cut according to the current mode)
void print_result(void) {
	//Set 6 significant digits explicitly (a DEBUG=2 log used to leave cout at 12)
	cout << setprecision(6);

	Cut c = do_cut();

	//List only clusters of size >= 2, ordered by size descending then cohesion descending
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

	//Objects that did not join a cluster (still singletons)
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

//Human-readable dendrogram. Depth-limited; deeper parts collapse to one line, use --matrix for the full data.
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
		cout << "... (" << r.size << " objects; use --matrix for the full data)" << endl;
		return;
	}
	print_subtree(r.left, depth + 1, max_depth);
	print_subtree(r.right, depth + 1, max_depth);
}

void print_tree(void) {
	const vector<Merge>& tr = g_hac->tree();
	const size_t n = g_hac->num_objects();
	//A root is a node that never appeared as the child of any merge
	vector<char> is_child(n + tr.size(), 0);
	for (size_t i = 0; i < tr.size(); i++) {
		is_child[tr[i].left] = 1;
		is_child[tr[i].right] = 1;
	}
	cout << "Dendrogram (" << n << " objects, " << tr.size() << " merges; "
		<< "the threshold only decides where to cut; the tree itself is complete):" << endl;
	const int max_depth = tree_depth;
	size_t roots = 0;
	for (size_t node = 0; node < n + tr.size(); node++) {
		if (!is_child[node]) {
			if (roots > 0)
				cout << endl;      // blank line between roots (disconnected components)
			print_subtree(node, 0, max_depth);
			roots++;
		}
	}
	cout << endl << "(" << roots << " roots in total: clusters at similarity 0 never merge)" << endl;
}

//The linkage matrix: one row per merge as left right height size, with scipy's
//node-id convention (0..n-1 are the objects, n+i is merge i). Note that height
//here is a **similarity** (higher merges first); scipy uses a distance.
void print_matrix(void) {
	const vector<Merge>& tr = g_hac->tree();
	const size_t n = g_hac->num_objects();
	cout << "# linkage matrix: left right height size" << endl;
	cout << "# node ids: 0.." << (n ? n - 1 : 0) << " are objects (in order of first appearance in the input), "
		<< n << "+i is merge i. height is a **similarity** (higher merges first), not a distance." << endl;
	cout << setprecision(12);
	for (size_t i = 0; i < tr.size(); i++) {
		cout << tr[i].left << " " << tr[i].right << " "
			<< tr[i].similarity << " " << tr[i].size << endl;
	}
}

//scipy-compatible export: one row per merge as left right distance size.
//Note that height is a **distance** (smaller merges first), the opposite of the
//similarity printed by --matrix.
void print_scipy(void) {
	spclust::ScipyOptions o;
	o.map = scipy_map;
	o.bridge = scipy_bridge;
	spclust::ScipyExport e = spclust::to_scipy(g_hac->tree(), g_hac->num_objects(), o);

	if (e.n_components > 1 && !e.bridged)
		cerr << "note: the merge tree is a forest (" << e.n_components
			<< " components that the tree never joined) - scipy requires exactly n-1 rows. "
			<< "Use --scipy-bridge to append bridging rows, or use the average/max linkage." << endl;
	if (!e.monotonic)
		cerr << "note: distances are not non-decreasing along the merge sequence (MinNonzero "
			<< "does not guarantee it), while scipy's dendrogram/fcluster assume they are." << endl;

	cout << "# scipy linkage matrix: left right distance size" << endl;
	cout << "# n=" << e.n_objects << " rows=" << e.rows.size()
		<< " components=" << e.n_components << (e.bridged ? " (bridged)" : "") << endl;
	cout << setprecision(12);
	for (size_t i = 0; i < e.rows.size(); i++)
		cout << e.rows[i].left << " " << e.rows[i].right << " "
			<< e.rows[i].distance << " " << e.rows[i].size << endl;
}

void usage(const char* prog) {
	cerr << "usage: " << prog << " [similarity.csv|-] [options]" << endl
	     << "  similarity.csv   input edge list; omit it or pass - to read stdin" << endl
	     << "  --cut <X>        cut the merge tree at similarity X (default " << MIN_SIM_CLUST
	     << ", i.e. the compile-time MIN_SIM_CLUST)" << endl
	     << "  --n-clusters <K> exactly K clusters (needs the full tree, see below)" << endl
	     << "  --min-size <M>   merge until no cluster is smaller than M (needs the full tree)" << endl
	     << "  --floor <X>      how far down to compute the tree (default = cut; 0 = full tree, 2-4x slower)" << endl
	     << "  --tree           print a dendrogram" << endl
	     << "  --tree-depth <N> depth limit for the dendrogram (default " << 6 << ", 0 = unlimited; implies --tree)" << endl
	     << "  --linkage <M>    linkage: 0/average (default), 1/max|single, 2/min|min-nonzero" << endl
	     << "  --matrix         print the linkage matrix (left right height size; height is a similarity)" << endl
	     << "  --scipy          export a scipy.cluster.hierarchy linkage matrix (height is a distance)" << endl
	     << "  --scipy-distance <M>  complement (default, C-sim) | negate (-sim) | reciprocal (1/sim)" << endl
	     << "  --scipy-bridge   bridge a forest into one tree (bridging distance = inf); scipy needs exactly n-1 rows" << endl
	     << "  -h, --help       show this" << endl;
}

} // namespace

int main(int argc, char* argv[]) {
	const char* path = NULL;
	for (int i = 1; i < argc; i++) {
		const string a = argv[i];
		if ("--cut" == a) {
			if (++i >= argc) {
				cerr << "--cut needs a number" << endl;
				return 2;
			}
			cut_sim = atof(argv[i]);
		} else if ("--n-clusters" == a) {
			if (++i >= argc) {
				cerr << "--n-clusters needs a number" << endl;
				return 2;
			}
			cut_n_clusters = (size_t)strtoul(argv[i], NULL, 10);
			cut_mode = CUT_N;
		} else if ("--min-size" == a) {
			if (++i >= argc) {
				cerr << "--min-size needs a number" << endl;
				return 2;
			}
			cut_min_size = (size_t)strtoul(argv[i], NULL, 10);
			cut_mode = CUT_SIZE;
		} else if ("--floor" == a) {
			if (++i >= argc) {
				cerr << "--floor needs a number" << endl;
				return 2;
			}
			tree_floor = atof(argv[i]);
			floor_given = true;
		} else if ("--tree" == a) {
			show_tree = true;
		} else if ("--tree-depth" == a) {
			if (++i >= argc) {
				cerr << "--tree-depth needs a number" << endl;
				return 2;
			}
			tree_depth = atoi(argv[i]);
			show_tree = true;        // asking for a depth obviously means the tree is wanted

		} else if ("--linkage" == a) {
			if (++i >= argc) {
				cerr << "--linkage needs a value" << endl;
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
				cerr << "unknown linkage: " << v << " (use 0/average, 1/max, 2/min)" << endl;
				return 2;
			}

		} else if ("--scipy" == a) {
			show_scipy = true;
		} else if ("--scipy-distance" == a) {
			if (++i >= argc) {
				cerr << "--scipy-distance needs a value" << endl;
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
				cerr << "unknown distance mapping: " << v << " (complement/negate/reciprocal)" << endl;
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
			//Anything starting with - that is not a bare "-" (that one means stdin) is
			//an unknown option. This used to treat --nope as a filename, so it was
			//reported as "cannot open file" (exit code 1), conflating a usage error
			//with an input error.
			cerr << "unknown option: " << a << endl;
			usage(argv[0]);
			return 2;
		} else if (NULL == path) {
			path = argv[i];
		} else {
			cerr << "unexpected extra argument: " << a << endl;
			usage(argv[0]);
			return 2;
		}
	}

	//The days of one process holding one state are over: the library uses objects, so make one here.
	spclust::HAC hac(linkage);
	g_hac = &hac;

	//No argument or "-" both mean read from standard input.
	//Edges go into the structure as they are read (the whole input is never
	//buffered), so the init phase is folded in here.
	//A return of 0 means failure (readdb already wrote the reason to stderr).
	//readdb used to call exit(1) itself, which made it unusable as a library;
	//now the exit code is decided here.
	if (0 == readdb(path, &edge_to_library))
		return 1;

	//The tree only has to be computed deeper when --cut is below the default: with
	//no --floor given it follows automatically, so "just cut at a lower threshold"
	//needs no extra argument and nothing is computed needlessly.
	//Cutting by cluster count or minimum size requires the full tree (there is no
	//way to know how far down merging has to go); see usage for the cost.
	if (!floor_given) {
		if (CUT_SIM != cut_mode || show_scipy)
			tree_floor = 0.0;   // k clusters / a minimum size / a scipy matrix all need the full tree
		else if (cut_sim < tree_floor)
			tree_floor = cut_sim;
	}

	//This prints the **effective** floor, which may have just been lowered above
	if (DEBUG >= 1)
		cout << "Linkage=" << hac.linkage_name() << ", " << cut_desc()
			<< ", floor=" << tree_floor << endl << endl;

	//Compute the merge tree. Timing lives in the CLI - the library prints nothing.
	if (DEBUG >= 1)
		cout << "Clustering..." << endl;
	clock_t s_time = clock();
	hac.build(tree_floor);
	if (DEBUG >= 1)
		cout << "Done. (" << (double)(clock() - s_time) / CLOCKS_PER_SEC << "s)" << endl << endl;

	//Print the flat result (the tree cut according to the current mode)
	print_result();

	if (show_tree)
		print_tree();
	if (show_matrix)
		print_matrix();
	if (show_scipy)
		print_scipy();

	return 0;
}

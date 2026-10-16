// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

#include "spclust/hac.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------
// AddressSanitizer memory poisoning. The node pool is hand-written malloc, so freed
// slots are poisoned and "reading a released node" is still caught by ASan - one of
// the reasons the pool cannot simply become a std::vector.
// Note the macros take sizeof(sp_unit) at the point of use, so defining them at the
// top of the file is fine.
// ---------------------------------------------------------------------------
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


using namespace std;


namespace spclust {
namespace detail {

//Initial "mean intra-cluster similarity" of a singleton group.
//
//0 is a fine initial value for the average/max recurrences, but the min recurrence is
//min(min(g1,g2), max_sim) and 0 is an absorbing element there, which would make the
//cohesion report identically 0. A singleton has no intra-cluster edge at all, so its
//initial value should be the **+inf of the empty set**, and only then does min mean
//anything.

//The node that stores one similarity during clustering.
//
//The four links and both group numbers are **indices**, not pointers, so the whole
//pool is one contiguous block:
//   sizeof(sp_unit) = 4*4(links) + 2*4(groups) + 8(sim) = 32 bytes (was 48)
//Node size turned out to be a first-order factor (padding a node to 128 bytes costs
//1.17x), so this compression is worth it.
//
//An empty link is 0. Node 0 is a **reserved null sentinel** that belongs to no list,
//so value-initialization (zeroing) of vec_dim1/vec_dim2 is naturally "empty" - the
//same behaviour as when head/max were raw pointers and zero meant NULL.
const uint32_t SP_NONE = 0;

//-- Everything below is a stateless type definition shared by the three linkages --

struct SimEntry {
	double sim;
	size_t row;
	size_t version;
};

enum { F_UPDATED = 1, F_DIRTY_MAX = 2, F_SIM_G1 = 4, F_SIM_G2 = 8 };

struct krusk_edge {
	uint32_t a;
	uint32_t b;
	double sim;
};

struct sp_unit {
	uint32_t right;
	uint32_t down;
	uint32_t left;
	uint32_t up;
	uint32_t group1;
	uint32_t group2;
	double sim;
};
//An element of the "row vector"
struct p_dim1 {
	uint32_t head;
	uint32_t max;
};
//An element of the "column vector"
struct p_dim2 {
	uint32_t head;
};
//A live cluster during clustering: mean intra-cluster similarity + member count.
//The member list is **not** stored here: members come from the merge tree plus the
//caller's id-to-name mapping, and keeping a list would only duplicate memory.
struct grp_s {
	double sim;
	std::size_t n_members;
};

//The "row vector": for every row, the head of its linked list and its largest element
// ===========================================================================
// The implementation. This used to be 34 file-static variables plus static functions,
// which allowed exactly one clustering per process. All of it now lives in one struct,
// so:
//   * a process can run any number of independent clusterings;
//   * running several in parallel needs no locking (as long as no object is shared);
//   * the linkage is a member set at construction instead of a compile-time constant
//     (it used to be #define LINKAGE_MTD).
// Every internal symbol lives in spclust::detail; the header only forward-declares
// HACImpl.
// ===========================================================================
struct HACImpl {
	explicit HACImpl(Linkage mtd)
		: linkage(mtd),
		  empty_cohesion(mtd == Linkage::MinNonzero
						 ? numeric_limits<double>::infinity() : 0.0) {}

	const Linkage linkage;
	// Initial "mean intra-cluster similarity" of a singleton group.
	//
	// 0 is a fine initial value for the average/max recurrences, but the min recurrence is
	// min(min(g1,g2), max_sim) and 0 is absorbing there, which makes the cohesion report
	// identically 0. A singleton has no intra-cluster edge, so its init is the **+inf of the empty set**.
	const double empty_cohesion;

vector<p_dim1> vec_dim1;
//The "column vector": for every column, the head of its (vertical) linked list
vector<p_dim2> vec_dim2;

// ---------------------------------------------------------------------------
// The node pool.
//
// Every cell lives in one contiguous block (index-addressed); released nodes are recycled through a free list.
// Replacing malloc/free itself saves almost nothing (measured 1.00x); **shrinking the node does**:
//   sizeof(sp_unit) 48 -> 32 bytes
//
// Manual ASan poisoning: the pool stays accessible for its whole lifetime, so ASan cannot see
// "a cell was released and then read" on its own - the bug class that appeared twice here.
// So a node is poisoned on release and unpoisoned on allocation, and use-after-free is still caught.
// When the pool moves, poisoning recorded on the old buffer is invalidated and redone via the free list.
// ---------------------------------------------------------------------------
//The pool uses malloc/realloc rather than std::vector: growing a vector means new block + copy + free old,
//and the allocator does not return large blocks promptly, so every move leaves a resident intermediate block -
//measured at 800k edges, growth by 1.25x costs 150 MB (for 26 MB of nodes total!), 2x costs 85 MB.
//realloc can usually extend large blocks in place: no copy, no intermediate block.
sp_unit* arena = NULL;
size_t arena_cap = 0;    // capacity in nodes
size_t arena_used = 0;   // highest index ever handed out (including the 0 sentinel)
//Free indices are stored separately rather than borrowed from the node's own fields -
//otherwise walking the free list would read poisoned nodes, and the poisoning exists
//precisely to catch that.
vector<uint32_t> arena_free_list;
size_t arena_live = 0;


uint32_t node_alloc(void) {
	if (0 == arena_used) {
		arena = (sp_unit*)calloc(1, sizeof(sp_unit));   // index 0 is the null sentinel, never handed out
		if (NULL == arena)
			throw std::bad_alloc();
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
		//realloc may move the block, and moving reads the old contents - which may
		//include poisoned free nodes, so unpoison the whole block first
		SP_UNPOISON_RANGE(arena, arena_used * sizeof(sp_unit));
		sp_unit* p = (sp_unit*)realloc(arena, want * sizeof(sp_unit));
		if (NULL == p)
			throw std::bad_alloc();
		arena = p;
		arena_cap = want;
		for (size_t k = 0; k < arena_free_list.size(); k++)
			SP_POISON(&arena[arena_free_list[k]]);
	}
	arena_live++;
	return (uint32_t)arena_used++;
}

void node_free(uint32_t i) {
	arena_live--;
	arena_free_list.push_back(i);
	SP_POISON(&arena[i]);
}


// merge_grp only has to recompute the rows/columns adjacent to group1/group2, which is
// O(degree) indices, but the old implementation scanned the whole index space three
// times to find them (two update passes plus one max repair): those three O(n) sweeps
// were about 75% of the clustering phase, and each merge also allocated three fresh
// O(n) vectors. The buffers below are reused across merges and only reset at the
// positions actually touched, so every step is O(degree).
vector<unsigned char> mark_flags;  // bit0=recompute similarity, bit1=max invalid, bit2/bit3=similarity to group1/group2 recorded
vector<size_t> updated;            // indices (row or column) whose similarity must be recomputed
vector<size_t> dirty_max;          // rows whose max pointer is invalid and must be rescanned
vector<size_t> touched;            // positions whose bit was set, for O(k) reset

//The merge tree: one row appended per successful merge. This is the complete output of
//the hierarchical clustering; the threshold (--cut / cut()) is only where the tree is
//cut, no longer a stopping condition.
vector<Merge> merges;
size_t n_objects = 0;

//Similarity from index x to group1 / group2 in the current merge. Similarities are
//non-negative, so -1 safely means "there is no such edge".
//
//Why cache these two values: updating one index used to call get_unit() twice, and
//get_unit() is a linear walk along a chain that measured about half of the whole
//program's time. But those two edges are already on the four chains we traverse anyway:
//for index x, the cell S(x, group1) sits on the "row group1" chain when x > group1 and on
//the "column group1" chain when x < group1, and likewise for group2. Recording them
//during the collection phase removes get_unit() from the update phase entirely. The cost
//is two O(number of groups) double arrays, reset only where they were touched.
vector<double> sim_g1;
vector<double> sim_g2;


//Set one bit for index x; returns true when this bit was the first one set this round
bool set_flag(size_t x, unsigned char bit) {
	if (x >= mark_flags.size()) {
		mark_flags.resize(x + 1, 0);
		sim_g1.resize(x + 1, -1.0);
		sim_g2.resize(x + 1, -1.0);
	}
	if (0 == mark_flags[x])
		touched.push_back(x);          // record it so the round can reset in O(k)
	bool fresh = (0 == (mark_flags[x] & bit));
	mark_flags[x] |= bit;
	return fresh;
}

//Mark index x as needing a similarity recomputation; true means this is the first such mark (dedup)
bool mark_updated(size_t x) { return set_flag(x, F_UPDATED); }

//Mark index x's max pointer as invalid
bool mark_dirty_max(size_t x) { return set_flag(x, F_DIRTY_MAX); }

//Record the similarity from index x to group1 / group2. **The first write wins**:
//get_unit() walks the chain from the head (add_grp inserts at the head, i.e. the last
//cell inserted) and takes the first match, so the first write corresponds to the same
//cell - which keeps behaviour identical for inputs with duplicate pairs.
//Record the similarity from index x to group1 / group2. This is a **maximum**: the same
//pair can be on the structure more than once (duplicate input pairs), and the documented
//semantic is that the largest value wins - which must hold for derived similarities too,
//not just for the merge value itself. Taking the first cell the chain walk happens to
//reach instead would make the result depend on the order add_edge was called in.
//sim_g1/sim_g2 start at -1 (the "no edge" sentinel) and values are non-negative, so the
//first write still initialises the slot.
void record_sim_g1(size_t x, double v) {
	set_flag(x, F_SIM_G1);
	if (sim_g1[x] < v)
		sim_g1[x] = v;
}
void record_sim_g2(size_t x, double v) {
	set_flag(x, F_SIM_G2);
	if (sim_g2[x] < v)
		sim_g2[x] = v;
}

// ---------------------------------------------------------------------------
// The global maximum similarity: a max-heap with version-number lazy deletion.
//
// After the earlier optimisation, get_max()'s O(number of groups) sweep was 57-78% of
// the clustering phase; with a heap it is O(1).
//
// Key design point: **the key is frozen inside the heap entry** (the similarity is
// stored in the entry rather than read live from vec_dim1[row].max). This is
// mandatory: with a live key, changing k rows' keys during one merge leaves the heap
// array's order no longer matching the current keys, so the heap invariant breaks -
// and sift_up/sift_down both assume "the heap is valid everywhere except the one place
// being repaired", so sifting each position in place cannot restore it. That one was
// found by measurement.
//
// When a row's max changes, its version is bumped and a new entry is pushed; old entries
// are invalidated by their version. The benefits:
//   * no bookkeeping of "where an entry sits in the heap" is needed, and updates need not
//     be idempotent;
//   * it does not depend on call order, so it can be maintained in place in add_grp, in
//     removal, and in the max rescan;
//   * heap operations never touch vec_dim1, so the max pointers that dangle temporarily
//     mid-merge are irrelevant to the heap.
// Once invalidated entries outnumber live rows by a factor, the heap is rebuilt with one
// entry per live row: memory stays O(n) and the rebuild costs O(heap size), i.e. O(1)
// amortised.
// ---------------------------------------------------------------------------

vector<SimEntry> heap;
vector<size_t> row_version;   // current version number of each row
vector<char> row_live;        // whether the row is currently non-empty (drives the rebuild threshold)
size_t live_rows = 0;

//Heap order: compare similarity first; on a tie the smaller row index counts as
//"larger". Together with the column tie-break in the max repair this makes the merged
//pair the smallest (left, right) id pair among those tied on similarity - the rule
//documented in spclust/hac.h. That the ids decide, rather than the order edges were
//added in, is what makes the result reproducible.
static bool entry_less(const SimEntry& a, const SimEntry& b) {
	if (a.sim != b.sim)
		return a.sim < b.sim;
	return a.row > b.row;
}

void heap_swap(size_t i, size_t j) {
	SimEntry tmp = heap[i];
	heap[i] = heap[j];
	heap[j] = tmp;
}

void heap_sift_up(size_t i) {
	while (i > 0) {
		size_t parent = (i - 1) / 2;
		if (!entry_less(heap[parent], heap[i]))
			break;
		heap_swap(parent, i);
		i = parent;
	}
}

void heap_sift_down(size_t i) {
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

//Row `row`'s max changed (it may have become empty): bump the version and, if still
//non-empty, push a new entry. **Every place that writes vec_dim1[*].max must call
//this**, otherwise that row is silently skipped.
void note_row_changed(size_t row) {
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
		return;                      // became empty: just invalidate the old entry

	SimEntry e;
	e.sim = arena[vec_dim1[row].max].sim;
	e.row = row;
	e.version = row_version[row];
	heap.push_back(e);
	heap_sift_up(heap.size() - 1);

	//Too many invalidated entries: rebuild, keeping only each live row's current entry,
	//then Floyd heapify. The cost scales with the heap size, so it is O(1) amortised
	//against the number of pushes and memory stays O(live rows).
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

//The table of live clusters. It used to be defined in main.cpp and wired up here with
//extern; it now lives in this object (the CLI no longer reads it).
vector<grp_s> grp;

//Insert a (group1, group2, sim) cell into the structure, always at the head of its list
//Decide which of two candidate cells of the **same row** the row's max pointer should hold:
//higher similarity wins, and on a tie the smaller column id. This one rule, applied in
//add_grp (initial insertion) and in the max repair, is what makes the merged pair the
//smallest (left, right) id pair among those tied on similarity - the rule documented in
//spclust/hac.h. Leaving either site with a strict '>' would let insertion order decide.
static bool beats(double sim_a, uint32_t col_a, double sim_b, uint32_t col_b) {
	if (sim_a != sim_b)
		return sim_a > sim_b;
	return col_a < col_b;
}

uint32_t add_grp(size_t group1, size_t group2, double sim) {
	//group1 is always kept smaller than group2
	if (group1 > group2)
		swap(group1, group2);

	//grow the row vector and the column vector
	if (vec_dim1.size() <= group1)
		vec_dim1.resize(group1 + 1);
	if (vec_dim2.size() <= group2)
		vec_dim2.resize(group2 + 1);

	//the cell to insert (take an index from the pool)
	uint32_t unit = node_alloc();
	//Assign once (field order matches sp_unit): the compiler emits a few wide stores,
	//saving instructions over seven scattered field writes.
	arena[unit] = sp_unit{ SP_NONE, SP_NONE, SP_NONE, SP_NONE,
		(uint32_t)group1, (uint32_t)group2, sim };

	//If this row has no max yet, or the new cell beats the current one, update the row
	//vector's max pointer. "Beats" means higher similarity, and on a tie the **smaller
	//column id** - the same rule the max repair uses, so that the row's max is decided by
	//ids rather than by which edge was added first.
	if (SP_NONE == vec_dim1[group1].max ||
		beats(sim, (uint32_t)group2, arena[vec_dim1[group1].max].sim,
			  arena[vec_dim1[group1].max].group2)) {
		vec_dim1[group1].max = unit;
		note_row_changed(group1);
	}

	//the previous row-list head
	uint32_t prev_head = vec_dim1[group1].head;
	//the new cell becomes the list head
	vec_dim1[group1].head = unit;
	//If there was a previous head, shift it right and relink it behind the new cell
	if (SP_NONE != prev_head) {
		arena[unit].right = prev_head;
		arena[prev_head].left = unit;
	}

	//the previous column-list head; the rest is as above
	prev_head = vec_dim2[group2].head;
	vec_dim2[group2].head = unit;
	if (SP_NONE != prev_head) {
		arena[unit].down = prev_head;
		arena[prev_head].up = unit;
	}
	return unit;
}

//Get the largest cell in the whole structure (the heap top, amortised O(1))
uint32_t get_max(void) {
	//discard heap tops invalidated by their version
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
	//As in the original: similarity 0 is not a candidate (the original filtered with a
	//max initialised to 0 plus a strict >). The heap holds every non-empty row, so a top
	//<= 0 means no candidate with positive similarity exists.
	if (SP_NONE == unit || !(arena[unit].sim > 0))
		return SP_NONE;
	return unit;
}

//Locate the (group1, group2, sim) cell
uint32_t get_unit(size_t group1, size_t group2) {
	//group1 is always kept smaller than group2
	if (group1 > group2)
		swap(group1, group2);
	//Out of range
	if (group1 >= vec_dim1.size() || group2 >= vec_dim2.size())
		return SP_NONE;
	//Walk the row list of group1 in vec_dim1
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
// Single linkage (Linkage::MaxSimilarity) goes through Kruskal and never builds the sparse structure.
//
// Under single linkage the inter-cluster similarity is the maximum over all **original
// edges** between two clusters (the max is not diluted), so the whole merge tree is
// equivalent to "sort the edges by similarity descending and merge with a union-find"
// (Kruskal / MST).
//
// That avoids the densification the original suffered: the max is not diluted, so the
// cluster-level similarity graph becomes the **transitive closure** of the input graph -
// measured at n=40000 with 400k edges it built 93,418,345 cells (234x the edge count) in
// 7.35 s. Kruskal uses m edges plus one sort.
//
// With ties (several edges of equal similarity) the merge **order** may differ from the
// greedy implementation, but under single linkage both implementations amount to "the
// connected components of the graph restricted to edges with sim >= X", which does not
// depend on the order. So the groups at any threshold, and each cluster's cohesion (the
// maximum intra-cluster edge), are identical.
// ---------------------------------------------------------------------------

vector<krusk_edge> krusk_edges;
vector<uint32_t> krusk_parent;

static bool krusk_order(const krusk_edge& x, const krusk_edge& y) {
	if (x.sim != y.sim)
		return x.sim > y.sim;      // similarity descending
	if (x.a != y.a)
		return x.a < y.a;          // break ties by id, purely for determinism
	return x.b < y.b;
}

uint32_t krusk_find(uint32_t x) {
	while (krusk_parent[x] != x) {
		krusk_parent[x] = krusk_parent[krusk_parent[x]];   // path halving
		x = krusk_parent[x];
	}
	return x;
}

void clust_kruskal(double floor) {
	sort(krusk_edges.begin(), krusk_edges.end(), krusk_order);

	const size_t n = n_objects;
	//A set's representative is the node id representing that cluster: initially the object
	//itself, and every merge produces n+i
	krusk_parent.resize(n);
	for (size_t i = 0; i < n; i++)
		krusk_parent[i] = (uint32_t)i;
	vector<size_t> size_of(n, 1);
	vector<double> cohesion_of(n, 0.0);

	size_t next_id = n;
	for (size_t i = 0; i < krusk_edges.size(); i++) {
		const krusk_edge& e = krusk_edges[i];
		if (e.sim < floor || !(e.sim > 0))
			break;                  // reached floor (or no edge with positive similarity left)
		uint32_t ra = krusk_find(e.a);
		uint32_t rb = krusk_find(e.b);
		if (ra == rb)
			continue;               // already in the same cluster (duplicate edge or cycle)

		uint32_t node = (uint32_t)next_id++;
		Merge rec;
		rec.left = ra;
		rec.right = rb;
		rec.similarity = e.sim;
		//cohesion = the maximum intra-cluster edge = the old max(max(g1,g2), max_sim)
		rec.cohesion = max(max(cohesion_of[ra], cohesion_of[rb]), e.sim);
		rec.size = size_of[ra] + size_of[rb];
		merges.push_back(rec);

		//both set representatives point at the new node
		krusk_parent[ra] = node;
		krusk_parent[rb] = node;
		krusk_parent.push_back(node);
		cohesion_of.push_back(rec.cohesion);
		size_of.push_back(rec.size);

	}
}

//Feed one edge. The library **stores no object names**: ids come from the caller (they
//are expected to be dense 0..n-1), and a new id is registered as a new singleton group.
void add_edge(uint32_t a, uint32_t b, double sim) {
	//A self-loop would break the internal "row < column" premise that all of add_grp relies on
	if (a == b)
		throw std::invalid_argument("spclust::add_edge: self-loop is not supported");
	//Similarity must be finite and non-negative; !(sim >= 0.0) also rejects NaN
	if (!(sim >= 0.0))
		throw std::invalid_argument("spclust::add_edge: similarity must be finite and non-negative");
	if (a > b)
		swap(a, b);

	if (grp.size() <= b) {
		//A new id is registered as a singleton group. If the caller skips ids (sparse
		//ids), the slots in between also become objects with no edges - see the id
		//contract in the header.
		size_t old_size = grp.size();
		grp.resize((size_t)b + 1);
		for (size_t i = old_size; i <= b; i++) {
			grp[i].sim = empty_cohesion;
			grp[i].n_members = 1;
		}
	}

	if (Linkage::MaxSimilarity == linkage) {
		//Single linkage does not need the sparse structure, only an edge list (Kruskal
		//comes later). The branch is a constant per instance, so the add_grp path is
		//never reached.
		krusk_edge e;
		e.a = a;
		e.b = b;
		e.sim = sim;
		krusk_edges.push_back(e);
	} else
		//insert group information (not object information) into the structure
		add_grp(a, b, sim);
}

size_t merge_grp(size_t group1, size_t group2, double max_sim) {
	//The row-vector length changes while the two groups are being merged, so record the
	//pre-merge length in dim1_size rather than using vec_dim1.size()
	size_t dim1_size = vec_dim1.size();
	//member counts of group1 and group2 before merging
	size_t group1_size = grp[group1].n_members;
	size_t group2_size = grp[group2].n_members;
	//length of the group vector
	size_t grp_size = grp.size();

	//Collect the rows/columns whose similarity must be recomputed, plus the rows whose max
	//pointer is invalid. Only marked indices are visited.

	//Every column holding a cell of the (group1, ...) chain needs updating
	uint32_t unit = vec_dim1[group1].head;
	while (SP_NONE != unit) {
		if (group2 != arena[unit].group2) {
			record_sim_g1(arena[unit].group2, arena[unit].sim);   // cell (group1,x) is S(x,group1)
			if (mark_updated(arena[unit].group2))
				updated.push_back(arena[unit].group2);
		}
		unit = arena[unit].right;
	}

	//Every column holding a cell of the (group2, ...) chain needs updating, if that chain exists
	if (dim1_size > group2) {
		unit = vec_dim1[group2].head;
		while (SP_NONE != unit) {
			record_sim_g2(arena[unit].group2, arena[unit].sim);   // cell (group2,x) is S(x,group2)
			if (mark_updated(arena[unit].group2))
				updated.push_back(arena[unit].group2);
			unit = arena[unit].right;
		}
	}

	//Every row holding a cell of the (..., group1, ...) chain needs updating
	unit = vec_dim2[group1].head;
	while (SP_NONE != unit) {
		record_sim_g1(arena[unit].group1, arena[unit].sim);       // cell (x,group1) is S(x,group1)
		if (mark_updated(arena[unit].group1))
			updated.push_back(arena[unit].group1);
		if (unit == vec_dim1[arena[unit].group1].max && mark_dirty_max(arena[unit].group1))
			dirty_max.push_back(arena[unit].group1);
		unit = arena[unit].down;
	}

	//Every row holding a cell of the (..., group2, ...) chain needs updating
	unit = vec_dim2[group2].head;
	while (SP_NONE != unit) {
		if (group1 != arena[unit].group1) {
			record_sim_g2(arena[unit].group1, arena[unit].sim);   // cell (x,group2) is S(x,group2)
			if (mark_updated(arena[unit].group1))
				updated.push_back(arena[unit].group1);
			if (unit == vec_dim1[arena[unit].group1].max && mark_dirty_max(arena[unit].group1))
				dirty_max.push_back(arena[unit].group1);
		}
		unit = arena[unit].down;
	}

	//Sort ascending to keep the same processing order as the old implementation (two
	//passes of for i = 0..dim over the marks).
	//Note these two arrays are sorted once per merge and are not tiny - measured 38
	//elements on average (1,402,219 iterations over 36,564 merges at n=80000). Switching
	//to insertion sort made the whole program 21% slower (O(k^2) does not pay off at
	//k ~ 38). Do not replace the std::sort here with a small-array sort again.
	sort(updated.begin(), updated.end());

	//Recompute these indices' similarity to the new group.
	//The old implementation had two passes, one over rows and one over columns, with
	//byte-identical loop bodies: for the same index x both do get_unit(x, group1) and
	//get_unit(x, group2). So when x appeared in both passes, the old code called add_grp
	//twice for the same (x, grp_size) and inserted two duplicate cells with equal values.
	//Merging them into one pass removes those duplicates (case 20_overlap_update.csv is
	//exactly that configuration).
	for (size_t k = 0; k < updated.size(); k++) {
		size_t i = updated[k];
		//The collection phase already recorded the values on those four chains, so read them
		//directly instead of doing get_unit()'s linear walk. -1 means there is no such edge,
		//equivalent to the old get_unit returning NULL -> 0.
		double sim1 = (sim_g1[i] < 0) ? 0 : sim_g1[i];
		double sim2 = (sim_g2[i] < 0) ? 0 : sim_g2[i];

		//Average Linkage
		if (Linkage::Average == linkage) {
			double sim = (sim1 * group1_size + sim2 * group2_size) / 
				(group1_size + group2_size);
			if (sim > 0)
				add_grp(i, grp_size, sim);
		//MaxSimilarity (= single linkage in distance terms)
		} else if (Linkage::MaxSimilarity == linkage) {
			//As in the other two branches, only candidates > 0 are kept. A value of 0 means
			//there is no edge between the two clusters at all, and get_max() never picks a
			//cell with similarity <= 0 anyway, so dropping them does not change the merge
			//order - but it turns the matrix from "filled up on every merge" into a sparse
			//one: measured 116x faster at n=40000
			//（8.71s vs 0.075s）。
			double sim = max(sim1, sim2);
			if (sim > 0)
				add_grp(i, grp_size, sim);
		}
		//MinNonzero
		else if (Linkage::MinNonzero == linkage)
			if (min(sim1, sim2) > 0)
				add_grp(i, grp_size, min(sim1, sim2));
	}

	size_t new_grp = grp.size();
	//grow the group vector
	grp.resize(new_grp + 1);
	//the new slot first receives all objects of group1
	grp[new_grp].n_members = group1_size + group2_size;
	if (Linkage::Average == linkage)
		grp[new_grp].sim = (grp[group1].sim * group1_size * (group1_size - 1)
			+ grp[group2].sim * group2_size * (group2_size - 1) + group1_size * group2_size * max_sim * 2) 
			/ ((group1_size + group2_size) * (group1_size + group2_size - 1));
	else if (Linkage::MaxSimilarity == linkage)
		grp[new_grp].sim = max(max(grp[group1].sim, grp[group2].sim), max_sim);
	else if (Linkage::MinNonzero == linkage)
		grp[new_grp].sim = min(min(grp[group1].sim, grp[group2].sim), max_sim);

	//clear the objects of group1 and group2
	grp[group1].n_members = 0;
	grp[group1].sim = 0;
	grp[group2].n_members = 0;
	grp[group2].sim = 0;

	//locate the (group1, group2, ...) cell
	uint32_t unit_max = get_unit(group1, group2);
	uint32_t tmp = SP_NONE;

	//relink the lists
	//remove the whole (group1, ...) row
	unit = vec_dim1[group1].head;
	while (SP_NONE != unit) {
		//the (group1,group2,sim) cell itself is removed last
		if (unit != unit_max) {
			//no cell above it (it is the head of the column list)
			if (SP_NONE == arena[unit].up) {
				//and no cell below it either (it is the tail of the column list)
				if (SP_NONE == arena[unit].down)
					//the cell is alone in its column list: clear the column vector's head
					vec_dim2[arena[unit].group2].head = SP_NONE;
				else {
					//update the column vector's head
					vec_dim2[arena[unit].group2].head = arena[unit].down;
					//the cell below becomes the head of the column list
					arena[arena[unit].down].up = SP_NONE;
				}
			} else {
				//the cell above links directly to the cell below
				arena[arena[unit].up].down = arena[unit].down;
				//if there is a cell below
				if (SP_NONE != arena[unit].down)
					//the lower cell links directly to the upper one
					arena[arena[unit].down].up = arena[unit].up;
			}
			//advance one cell along the row path
			tmp = arena[unit].right;
			//release the current cell
			node_free(unit);
			unit = tmp;
		} else
			unit = arena[unit].right;
	}

	//remove the whole (group2, ...) row; it may not exist, so check first
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

	//remove the whole (..., group1, ...) column
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

	//remove the whole (..., group2, ...) column
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

	//finally remove the (group1,group2,sim) cell
	node_free(unit_max);
	unit_max = SP_NONE;
	
	//drop the group1/group2 entries from the row vector (if they exist)
	vec_dim1[group1].head = SP_NONE;
	vec_dim1[group1].max = SP_NONE;
	note_row_changed(group1);
	if (dim1_size > group2) {
		vec_dim1[group2].head = SP_NONE;
		vec_dim1[group2].max = SP_NONE;
		note_row_changed(group2);
	}
	//drop the group1/group2 entries from the column vector
	vec_dim2[group1].head = SP_NONE;
	vec_dim2[group2].head = SP_NONE;

	//Walk the rows whose max pointer must be repaired, ascending, only those invalidated now
	sort(dirty_max.begin(), dirty_max.end());
	for (size_t k = 0; k < dirty_max.size(); k++) {
		size_t i = dirty_max[k];
		uint32_t u = vec_dim1[i].head;
		uint32_t best = SP_NONE;
		//Walk the row and find the maximum. On a tie the **smaller column id** wins, not
		//the first cell encountered on the chain (see `beats`).
		while (SP_NONE != u) {
			if (SP_NONE == best || beats(arena[u].sim, arena[u].group2, arena[best].sim,
									   arena[best].group2))
				best = u;
			u = arena[u].right;
		}
		//This assignment must be unconditional. The old form was `if (arena[u].sim > max)`
		//with max initialised to 0: if a whole row had similarity 0 it never assigned, and
		//vec_dim1[i].max kept a pointer to freed memory, so the next get_max() was a
		//heap-use-after-free (see tests/cases/22_zero_row_max.csv).
		vec_dim1[i].max = best;
		note_row_changed(i);
	}

	//Reset the marks and cached values touched this round. The buffers are reused across
	//merges, so nothing is allocated here and the whole O(n) memory is not cleared.
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

//Release the group table. Safe to call twice (clear+swap leaves an empty vector).
void release_groups(void) { release(grp); }

//Release the memory used by the structure
//Shrink a vector to zero capacity. Plain clear() only drops the elements and keeps the
//capacity, which is exactly how ~220 MB per held result stayed resident at n = 1M: every
//one of the O(n) workspace buffers below is a vector, and after build() nothing reads any
//of them again. The cost of re-growing them in a later build() on the same object is a
//handful of large allocations, which the measurement below shows is negligible.
template <typename T>
static void release(std::vector<T>& v) {
	std::vector<T>().swap(v);
}

void free_mem(void) {
	//The node pool is one contiguous block: hand it back whole, no per-node traversal.
	//Unpoison everything first so no poisoning is left on memory returned to the allocator.
	if (NULL != arena)
		SP_UNPOISON_RANGE(arena, arena_used * sizeof(sp_unit));
	free(arena);
	arena = NULL;
	arena_cap = 0;
	arena_used = 0;
	arena_free_list.clear();
	vector<uint32_t>().swap(arena_free_list);
	arena_live = 0;

	//clear the row and column vectors
	for (size_t i = 0; i < vec_dim1.size(); i++) {
		vec_dim1[i].head = SP_NONE;
		vec_dim1[i].max = SP_NONE;
	}
	for (size_t i = 0; i < vec_dim2.size(); i++)
		vec_dim2[i].head = SP_NONE;

	release(vec_dim1);
	release(vec_dim2);
	release(mark_flags);
	release(sim_g1);
	release(sim_g2);
	release(updated);
	release(dirty_max);
	release(touched);
	release(krusk_edges);
	release(krusk_parent);
	release(heap);
	release(row_version);
	release(row_live);
	live_rows = 0;
}

//Compute the merge tree. floor is the similarity below which merging stops (<= 0 means
//compute the whole tree).
void build(double floor) {
	//the initial number of groups (objects) = the number of leaves
	n_objects = grp.size();

	//The id space is bounded by 2n-1 (n objects plus at most n-1 merges), so allocate once:
	//set_flag/record_sim_* in the hot path then never resize inside a loop (which would
	//trigger realloc and, under ASan, more poisoning work). The bounds check is kept, but
	//it is now always false and essentially free.
	mark_flags.assign(2 * n_objects, 0);
	sim_g1.assign(2 * n_objects, -1.0);
	sim_g2.assign(2 * n_objects, -1.0);

	if (Linkage::MaxSimilarity == linkage) {
		//single linkage: edges sorted by similarity descending + union-find, O(m log m)
		clust_kruskal(floor);
		free_mem();
		release_groups();
		return;
	}

	size_t grp_size = grp.size();
	//Merge until the largest similarity is below floor. floor comes from the caller:
	//  floor = the cut threshold -> costs exactly the same as stopping at that threshold
	//  floor <= 0                -> the full tree, at the cost of many more merges
	while (grp_size > 1) {
		uint32_t unit = get_max();
		if (SP_NONE == unit)
			break;                      // no candidate with positive similarity: the tree ends here
		double max_sim = arena[unit].sim;
		if (max_sim < floor)
			break;                      // reached floor; lower floor for a deeper tree
		size_t group1 = arena[unit].group1;
		size_t group2 = arena[unit].group2;
		grp_size--;
		size_t new_grp = merge_grp(group1, group2, max_sim) - 1;

		//record this merge (the group number is the node id; see spclust::Merge)
		Merge rec;
		rec.left = (uint32_t)group1;
		rec.right = (uint32_t)group2;
		rec.similarity = max_sim;
		rec.cohesion = grp[new_grp].sim;
		rec.size = grp[new_grp].n_members;
		merges.push_back(rec);
	}
	//Clustering finished: release the structure's memory, and the group table too.
	//`grp` has about 2n slots (16 bytes each) and nothing reads it after this point -
	//cut()/tree()/num_objects() all work off `merges` and `n_objects`. Keeping it would
	//leave ~32 bytes per object of dead resident memory in any process that holds the
	//result, which matters for a long-lived service or a Python session with several
	//HAC objects. This does not change peak memory: the peak is during the merge loop,
	//where the table is needed.
	free_mem();
	release_groups();
}

// ---------------------------------------------------------------------------
// Public API (see include/spclust/hac.h)
// ---------------------------------------------------------------------------

std::size_t num_objects() const {
	return n_objects;
}

const std::vector<Merge>& tree() const {
	return merges;
}

const char* linkage_name() const {
	switch (linkage) {
	case Linkage::Average:       return "average";
	case Linkage::MaxSimilarity: return "max-similarity";
	case Linkage::MinNonzero:    return "min-nonzero";
	}
	return "unknown";
}

void reserve(std::size_t n_objects_hint, std::size_t n_edges_hint) {
	//Reserve only what the caller explicitly knows. The node pool grows on demand rather
	//than taking one big block here - the pool's growth policy dominates peak memory.
	grp.reserve(n_objects_hint);
	merges.reserve(n_objects_hint);
	krusk_edges.reserve(n_edges_hint);
}

//Union-find used for cutting. The parent array belongs to the caller - it used to be a
//file static, so concurrent cuts in one process stepped on each other (and cut() was
//nominally const while mutating global state).
static std::size_t uf_find(std::vector<std::size_t>& parent, std::size_t x) {
	while (parent[x] != x) {
		parent[x] = parent[parent[x]];   // path halving
		x = parent[x];
	}
	return x;
}

// Treat the first `use` merges of the tree as the cut and collect every cluster's members
// and cohesion. cut / cut_n_clusters / cut_min_size all reduce to this step.
Cut partition_after(std::size_t use) const {
	const std::size_t n = n_objects;
	if (use > merges.size())
		use = merges.size();
	const std::size_t total = n + merges.size();

	std::vector<std::size_t> uf_parent(total);
	for (std::size_t i = 0; i < total; i++)
		uf_parent[i] = i;

	//Besides joining the two children, the node produced by this merge (id = n+i) must also
	//point into the same set - otherwise a later merge using n+i as a child would treat it
	//as a separate element.
	for (std::size_t i = 0; i < use; i++) {
		std::size_t a = uf_find(uf_parent, merges[i].left), b = uf_find(uf_parent, merges[i].right);
		if (a != b)
			uf_parent[b] = a;
		uf_parent[n + i] = uf_find(uf_parent, merges[i].left);
	}

	//A cluster's cohesion is the one written by the last merge of that cluster in the prefix
	//(later writes overwrite earlier ones)
	std::vector<double> root_cohesion(total, 0.0);
	for (std::size_t i = 0; i < use; i++)
		root_cohesion[uf_find(uf_parent, merges[i].left)] = merges[i].cohesion;

	const std::size_t NONE = (std::size_t)-1;
	std::vector<std::size_t> root_index(total, NONE);
	Cut res;
	for (std::size_t id = 0; id < n; id++) {
		std::size_t r = uf_find(uf_parent, id);
		if (NONE == root_index[r]) {
			root_index[r] = res.members.size();
			res.members.push_back(std::vector<std::uint32_t>());
			res.cohesion.push_back(root_cohesion[r]);
		}
		res.members[root_index[r]].push_back((std::uint32_t)id);
	}
	return res;
}

// Cut at similarity min_similarity: similarity is non-increasing along the merge sequence,
// so cutting is taking the prefix with similarity >= min_similarity.
Cut cut(double min_similarity) const {
	std::size_t k = 0;
	while (k < merges.size() && merges[k].similarity >= min_similarity)
		k++;
	return partition_after(k);
}

// Exactly k clusters: equivalent to taking the first n-k merges. If the tree was truncated
// with --floor it may not reach k clusters, in which case the finest cut the tree supports
// is returned.
Cut cut_n_clusters(std::size_t k) const {
	const std::size_t n = n_objects;
	if (k >= n)
		return partition_after(0);              // more than the object count: everything separate
	if (k < 1)
		k = 1;                                  // 0 clusters is meaningless; treat it as 1
	return partition_after(n - k);
}

// Merge until no cluster is smaller than min_size (the **shallowest** cut satisfying that).
// min_size <= 1 always holds, so nothing merges and every cluster is a singleton.
// Note this is **not** "discard small clusters" but "keep merging until there are none";
// if disconnected components remain (isolated points at similarity 0) they never join in,
// and then the finest cut the tree supports is returned.
Cut cut_min_size(std::size_t min_size) const {
	if (min_size <= 1)
		return partition_after(0);

	const std::size_t n = n_objects;
	const std::size_t total = n + merges.size();
	std::vector<std::size_t> parent(total), sz(total, 0);
	for (std::size_t i = 0; i < n; i++) {
		parent[i] = i;
		sz[i] = 1;
	}
	for (std::size_t i = n; i < total; i++)
		parent[i] = i;                          // merged nodes are aliases and count no elements of their own (sz=0)

	// small = how many clusters are still below min_size. Merging only ever grows clusters,
	// so this quantity is non-increasing and "the first position where none is below
	// min_size" is the answer - no need to rescan all clusters each time.
	std::size_t small = n;
	std::size_t use = merges.size();
	for (std::size_t i = 0; i < merges.size(); i++) {
		std::size_t a = uf_find(parent, merges[i].left);
		std::size_t b = uf_find(parent, merges[i].right);
		if (a != b) {
			if (sz[a] < min_size) small--;
			if (sz[b] < min_size) small--;
			parent[b] = a;
			sz[a] += sz[b];
			if (sz[a] < min_size) small++;
		}
		parent[n + i] = a;                      // the new node points at the merged cluster
		if (0 == small) {
			use = i + 1;
			break;
		}
	}
	return partition_after(use);
}

void reset() {
	free_mem();
	release_groups();
	release(merges);
	n_objects = 0;
}

}; // struct HACImpl

} // namespace detail

// ===========================================================================
// Public API: all thin wrappers, with no virtual function anywhere. Dispatch happens only
// at these entry points, while the inner loops (merge_grp / get_max / heap operations)
// stay inside one function where the linkage is known, so choosing the linkage at run time
// leaves no branch in the hot path.
// ===========================================================================

HAC::HAC(Linkage mtd) : impl_(new detail::HACImpl(mtd)) {}
HAC::~HAC() = default;
HAC::HAC(HAC&&) noexcept = default;
HAC& HAC::operator=(HAC&&) noexcept = default;

void HAC::reserve(std::size_t n, std::size_t m) { impl_->reserve(n, m); }
void HAC::add_edge(std::uint32_t a, std::uint32_t b, double similarity) {
	impl_->add_edge(a, b, similarity);
}
void HAC::build(double floor) { impl_->build(floor); }
std::size_t HAC::num_objects() const { return impl_->num_objects(); }
const std::vector<Merge>& HAC::tree() const { return impl_->tree(); }
Cut HAC::cut(double min_similarity) const { return impl_->cut(min_similarity); }
Cut HAC::cut_n_clusters(std::size_t k) const { return impl_->cut_n_clusters(k); }
Cut HAC::cut_min_size(std::size_t min_size) const { return impl_->cut_min_size(min_size); }
void HAC::reset() { impl_->reset(); }
Linkage HAC::linkage() const { return impl_->linkage; }
const char* HAC::linkage_name() const { return impl_->linkage_name(); }

} // namespace spclust

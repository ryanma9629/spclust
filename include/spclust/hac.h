// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

#ifndef SPCLUST_HAC_H
#define SPCLUST_HAC_H

// ===========================================================================
// Agglomerative hierarchical clustering over a sparse similarity matrix.
//
// Read this contract before using the library:
//
//   * The input is a **sparse** similarity matrix: any pair that is not given
//     has similarity 0 and costs no memory.
//   * Object ids are supplied by the caller and must be **dense 0..n-1** (the
//     internals allocate arrays indexed by the largest id). If you have
//     arbitrary or sparse ids, map them once on your side - that is exactly what
//     tools/spclust_cli.cpp does for string names, in about ten lines.
//   * If the same pair appears more than once, the **largest** value wins.
//   * A self-loop (a == b) would break the internal "row < column" invariant and
//     throws std::invalid_argument.
//   * Similarities must be finite and non-negative; negative values and NaN
//     throw std::invalid_argument. There is no upper bound - the historical
//     65535 limit was an artifact of an unsigned short representation.
//   * **The merge order of tied similarities is not specified.** Within one
//     linkage the behaviour is deterministic, but it may differ between linkages
//     and between versions. Under MaxSimilarity the resulting groups are
//     order-independent anyway (they are connected components); under Average and
//     MinNonzero they may differ.
//   * One HAC object holds one independent state; separate objects never affect
//     each other. A single object is **not thread-safe** - not even reading it
//     concurrently with writing.
//
// Relation to scipy: tree() is the linkage matrix, one row per merge:
//     left right similarity size
// Node numbering follows scipy's convention (0..n-1 are the objects, n+i is the
// cluster created by merge i). **Note that "height" here is a similarity (higher
// merges first), whereas scipy uses a distance.** See spclust/scipy.h for export.
// ===========================================================================

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace spclust {

// Linkage methods. Names follow the **similarity** semantics (as documented by
// sparseAHC):
//
//   Average       - the mean over **all |A|x|B| pairs** between two clusters
//                   (missing edges count as 0). The denominator is the full pair
//                   count, so larger clusters get diluted.
//   MaxSimilarity - take the maximum. In distance terms this is single linkage:
//                   one edge is enough to join two clusters, so it chains. Its
//                   inter-cluster similarity can only ever be an original edge
//                   value, so it runs through a Kruskal path and never builds the
//                   matrix - much faster on large sparse graphs.
//   MinNonzero    - take the minimum, but only over pairs where **both sides have
//                   an edge** (a missing pair does not participate). This is a
//                   conservative variant of complete linkage: naming it complete
//                   would be wrong, because treating missing edges as 0 makes the
//                   min identically 0 and meaningless on a sparse matrix.
enum class Linkage { Average, MaxSimilarity, MinNonzero };

// One merge - a row of the standard linkage matrix.
struct Merge {
	std::uint32_t left;      // left child
	std::uint32_t right;     // right child
	double similarity;       // similarity at which they merged (the cluster's height; higher merges first)
	double cohesion;         // mean intra-cluster similarity of the merged cluster
	std::size_t size;        // size of the merged cluster
};

// Result of cutting the tree: each cluster's members (object ids) and its mean
// intra-cluster similarity.
struct Cut {
	std::vector<std::vector<std::uint32_t>> members;
	std::vector<double> cohesion;
};

namespace detail {
struct HACImpl;   // defined in src/hac.cpp, so the header exposes no internals
}

// One hierarchical clustering. An object owns one independent state; any number
// of them can exist at the same time.
class HAC {
public:
	explicit HAC(Linkage linkage = Linkage::Average);
	~HAC();

	// Not copyable (it owns implementation state); movable.
	HAC(const HAC&) = delete;
	HAC& operator=(const HAC&) = delete;
	HAC(HAC&&) noexcept;
	HAC& operator=(HAC&&) noexcept;

	// Reserve capacity. Optional, but it avoids repeated growth when edges are
	// streamed in.
	void reserve(std::size_t n_objects, std::size_t n_edges);

	// Feed one edge. Missing pairs count as similarity 0; the largest value wins
	// when a pair repeats.
	void add_edge(std::uint32_t a, std::uint32_t b, double similarity);

	// Compute the merge tree.
	// floor is how far down the tree to compute: merging stops once the largest
	// remaining similarity is below floor. floor <= 0 computes the whole tree
	// (merging until no positive-similarity candidate is left).
	// Note this is the depth of the tree, not a clustering threshold; cut()
	// decides where to cut.
	void build(double floor = 0.0);

	// Number of objects (= number of leaves).
	std::size_t num_objects() const;

	// The merge tree, in merge order. Valid after build().
	const std::vector<Merge>& tree() const;

	// Cut the tree at similarity min_similarity, producing a flat partition.
	// Similarity is non-increasing along the merge sequence, so cutting is just
	// taking the prefix with similarity >= min_similarity.
	Cut cut(double min_similarity) const;

	// Exactly k clusters: equivalent to merging down to k clusters.
	// For k >= num_objects() everything stays separate; k == 0 is treated as 1.
	Cut cut_n_clusters(std::size_t k) const;

	// Only clusters of size >= min_size: take the shallowest prefix whose
	// clusters are all at least min_size (i.e. never cut out a cluster smaller
	// than min_size).
	Cut cut_min_size(std::size_t min_size) const;

	// Drop all state so the same object can be reused for another clustering.
	void reset();

	Linkage linkage() const;
	const char* linkage_name() const;

private:
	std::unique_ptr<detail::HACImpl> impl_;
};

} // namespace spclust

#endif

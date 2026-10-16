// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

#ifndef SPCLUST_SCIPY_H
#define SPCLUST_SCIPY_H

// ===========================================================================
// scipy interop: export this library's merge tree as a
// scipy.cluster.hierarchy linkage matrix.
//
// A scipy row is [idx1, idx2, height, count]:
//   * idx1/idx2 are the two children merged, **ascending**;
//   * height is a **distance** (smaller merges first), not a similarity - this is
//     the main difference between the two conventions;
//   * count is the number of samples in the merged cluster;
//   * node numbering matches ours: 0..n-1 are the leaves, n+i is merge i.
//
// The similarity -> distance mapping has to be chosen by the caller (a similarity
// is not a distance, and there is no single correct mapping):
//
//   Complement  dist = C - sim, with C defaulting to the largest similarity in the
//               tree. This is the **most faithful** one: this library's "a missing
//               edge is similarity 0" corresponds exactly to "distance C", so
//               Average linkage becomes scipy's method='average' applied to the
//               same dense distance matrix (measured: identical merge for merge,
//               differing only in tie order).
//   Negate      dist = -sim. Order-preserving and parameter-free, but distances
//               are negative.
//   Reciprocal  dist = 1/sim (+inf when sim is 0). Only meaningful when
//               similarities behave like affinities.
//
// **Two known differences** (semantics, not bugs):
//
//   1. `Linkage::MinNonzero` is not scipy's method='complete'. scipy reaches the
//      maximum distance (= similarity 0) on every missing pair, which on sparse
//      input makes almost everything look far apart; this library takes the
//      minimum only over pairs where both sides have an edge. **On dense input
//      the two agree** (measured).
//   2. The tree may be a **forest**: under `Linkage::MinNonzero` a pair of
//      clusters with no edge between them is not even a candidate, so components
//      may never merge. A scipy matrix must have exactly n-1 rows (scipy always
//      ends up with a single tree), so bridging rows can be requested explicitly:
//      the roots are merged in turn at a given distance. Nothing is bridged by
//      default; the caller decides based on `n_components`.
// ===========================================================================

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "spclust/hac.h"

namespace spclust {

enum class DistanceMap {
	Complement,   // C - sim
	Negate,       // -sim
	Reciprocal,   // 1/sim (+inf when sim <= 0)
};

// One row of the scipy linkage matrix.
struct ScipyRow {
	std::uint32_t left;      // the smaller child
	std::uint32_t right;     // the larger child
	double distance;
	std::uint32_t size;
};

struct ScipyOptions {
	DistanceMap map = DistanceMap::Complement;
	// The C used by Complement. 0 means "take the largest similarity in the tree".
	double reference = 0.0;
	// Whether to turn a forest into a single tree (scipy requires exactly n-1 rows)
	bool bridge = false;
	// Distance used for the bridging rows. NaN means automatic, which picks the
	// most faithful value for the mapping:
	//   Complement -> reference (C). With no edges between two components every
	//                 pair distance is C, which is exactly what scipy itself
	//                 would compute for them.
	//   Negate / Reciprocal -> +inf (nothing can be farther apart).
	double bridge_distance = std::numeric_limits<double>::quiet_NaN();
};

struct ScipyExport {
	std::vector<ScipyRow> rows;
	std::size_t n_objects = 0;
	std::size_t n_components = 0;   // number of roots in the tree; > 1 means a forest
	bool bridged = false;           // whether bridging rows were appended
	// Whether distances are non-decreasing along the merge sequence. scipy's
	// dendrogram and fcluster both assume this; Average and MaxSimilarity
	// guarantee it, MinNonzero does not.
	bool monotonic = true;
};

// Pure function: reads the tree, changes no state.
ScipyExport to_scipy(const std::vector<Merge>& tree, std::size_t n_objects,
					 const ScipyOptions& options = ScipyOptions());

// The distance function of a mapping under a given reference value (useful for
// converting thresholds yourself).
double map_distance(double similarity, DistanceMap map, double reference);

} // namespace spclust

#endif

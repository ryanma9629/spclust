// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

#include "spclust/scipy.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace spclust {

double map_distance(double similarity, DistanceMap map, double reference) {
	switch (map) {
	case DistanceMap::Negate:
		return -similarity;
	case DistanceMap::Reciprocal:
		return similarity > 0 ? 1.0 / similarity
							  : std::numeric_limits<double>::infinity();
	case DistanceMap::Complement:
		return reference - similarity;
	}
	return reference - similarity;
}

ScipyExport to_scipy(const std::vector<Merge>& tree, std::size_t n_objects,
					 const ScipyOptions& options) {
	ScipyExport out;
	out.n_objects = n_objects;

	// The Complement reference: the largest similarity in the tree by default. A monotone
	// tree would make it the first merge, but scanning all three is safer (MinNonzero need not be monotone).
	double reference = options.reference;
	if (DistanceMap::Complement == options.map && 0.0 == reference) {
		for (std::size_t i = 0; i < tree.size(); i++)
			if (tree[i].similarity > reference)
				reference = tree[i].similarity;
	}

	out.rows.reserve(tree.size() + n_objects);
	for (std::size_t i = 0; i < tree.size(); i++) {
		const Merge& m = tree[i];
		ScipyRow row;
		row.left = std::min(m.left, m.right);
		row.right = std::max(m.left, m.right);
		if (row.left == row.right)
			throw std::invalid_argument("to_scipy: a merge has identical children");
		row.distance = map_distance(m.similarity, options.map, reference);
		row.size = (std::uint32_t)m.size;
		out.rows.push_back(row);
	}

	// Monotonicity: scipy's dendrogram and fcluster both assume distances are non-decreasing.
	for (std::size_t i = 1; i < out.rows.size(); i++)
		if (out.rows[i].distance < out.rows[i - 1].distance)
			out.monotonic = false;

	// Number of components = nodes that never appeared as a child of any merge
	{
		const std::size_t total = n_objects + tree.size();
		std::vector<char> is_child(total, 0);
		for (std::size_t i = 0; i < tree.size(); i++) {
			is_child[tree[i].left] = 1;
			is_child[tree[i].right] = 1;
		}
		for (std::size_t node = 0; node < total; node++)
			if (!is_child[node])
				out.n_components++;
	}

	if (options.bridge && out.n_components > 1) {
		// Join the current roots one after another at the caller's bridging distance. scipy wants exactly n-1 rows,
		// and the last row's count must be n_objects.
		const std::size_t total = n_objects + tree.size();
		std::vector<char> is_child(total, 0);
		for (std::size_t i = 0; i < tree.size(); i++) {
			is_child[tree[i].left] = 1;
			is_child[tree[i].right] = 1;
		}
		std::vector<std::uint32_t> roots;
		for (std::size_t node = 0; node < total; node++)
			if (!is_child[node])
				roots.push_back((std::uint32_t)node);

		// Row i of the whole matrix represents node n+i, so the id of a newly merged node is
		// n + (current row count - 1). Leaves count 1; internal nodes are looked up in the rows already written.
		const std::size_t n = n_objects;
		std::uint32_t acc_id = roots[0];
		std::uint32_t acc_size = (roots[0] < n) ? 1u : out.rows[roots[0] - n].size;

		for (std::size_t i = 1; i < roots.size(); i++) {
			std::uint32_t sz = (roots[i] < n) ? 1u : out.rows[roots[i] - n].size;
			acc_size += sz;

			ScipyRow row;
			row.left = std::min(acc_id, roots[i]);
			row.right = std::max(acc_id, roots[i]);
			const double bd = (options.bridge_distance == options.bridge_distance)
								  ? options.bridge_distance            // not NaN
								  : (DistanceMap::Complement == options.map
										 ? reference
										 : std::numeric_limits<double>::infinity());
			row.distance = bd;
			row.size = acc_size;
			out.rows.push_back(row);

			// the id of the node this merge produces = n + the index of the row just written
			acc_id = (std::uint32_t)(n + out.rows.size() - 1);
		}
		out.bridged = true;
		// Bridging rows carry +inf (or the caller's value); if that is below the previous row, monotonicity is flagged too.
		for (std::size_t i = 1; i < out.rows.size(); i++)
			if (out.rows[i].distance < out.rows[i - 1].distance)
				out.monotonic = false;
	}

	return out;
}

} // namespace spclust

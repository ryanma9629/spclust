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

	// Complement 的参考值：默认取树里的最大相似度。树是单调的，所以就是第一次合并，
	// 但三种链接都扫一遍更稳妥（MinNonzero 不保证单调）。
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
			throw std::invalid_argument("to_scipy: 一次合并的两个子节点相同");
		row.distance = map_distance(m.similarity, options.map, reference);
		row.size = (std::uint32_t)m.size;
		out.rows.push_back(row);
	}

	// 单调性：scipy 的树状图与 fcluster 都假定距离沿合并序列非递减。
	for (std::size_t i = 1; i < out.rows.size(); i++)
		if (out.rows[i].distance < out.rows[i - 1].distance)
			out.monotonic = false;

	// 分量数 = 没当过任何一次合并子节点的节点数
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
		// 把当前的根依次并起来（距离用调用方给的桥接值）。scipy 要求恰好 n-1 行，
		// 且最后一行的 count 必须是 n_objects。
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

		// 第 i 行（全局下标）代表节点 n+i，所以合并出来的新节点 id 就是
		// n + (当前行数 - 1)。叶子大小算 1，内部节点查已经写好的行。
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
								  ? options.bridge_distance            // 非 NaN
								  : (DistanceMap::Complement == options.map
										 ? reference
										 : std::numeric_limits<double>::infinity());
			row.distance = bd;
			row.size = acc_size;
			out.rows.push_back(row);

			// 这次合并产生的节点 id = n + 刚写下的那一行的下标
			acc_id = (std::uint32_t)(n + out.rows.size() - 1);
		}
		out.bridged = true;
		// 桥接行是 +inf（或调用方给的值），如果它小于前一行，单调性也会被标出来。
		for (std::size_t i = 1; i < out.rows.size(); i++)
			if (out.rows[i].distance < out.rows[i - 1].distance)
				out.monotonic = false;
	}

	return out;
}

} // namespace spclust

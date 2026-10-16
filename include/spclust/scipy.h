// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

#ifndef SPCLUST_SCIPY_H
#define SPCLUST_SCIPY_H

// ===========================================================================
// scipy 互操作：把本库的合并树导出成 scipy.cluster.hierarchy 的 linkage 矩阵。
//
// scipy 的矩阵每行是 [idx1, idx2, height, count]：
//   * idx1/idx2 是这次合并的两个子节点，**升序**；
//   * height 是**距离**（越小越先合并），而非相似度 —— 这正是两边最大的差别；
//   * count 是合并后的样本数；
//   * 节点编号约定与我们相同：0..n-1 是叶子，n+i 是第 i 次合并。
//
// 相似度 -> 距离 的映射必须由调用方选（相似度不是距离，没有唯一正确映射）：
//
//   Complement  dist = C - sim，C 默认取树里的最大相似度。
//               这是**语义上最忠实**的一个：本库的"缺失边 = 相似度 0"正好对应
//               "距离 = C"，于是 average 链接等价于 scipy 的 method='average'
//               作用在同一个稠密距离矩阵上（已实测逐次一致，只差并列顺序）。
//   Negate      dist = -sim。保序、无参，但距离为负。
//   Reciprocal  dist = 1/sim（sim 为 0 时给 +inf）。只对"相似度像亲和度"的数据有意义。
//
// **两条已知差异**（不是 bug，是语义不同）：
//
//   1. `Linkage::MinNonzero` 与 scipy 的 method='complete' 不同。scipy 的 complete
//      在每个缺失对上取到最大距离（= 相似度 0），于是稀疏输入下几乎处处为 0；
//      本库只对"两边都存在的边"取最小。**稠密输入下两者相同**（已实测）。
//   2. 树可能是**森林**：`Linkage::MinNonzero` 下两簇之间根本没有边时，这个 pair
//      连候选都不是，于是分量可能永远合不起来。scipy 的矩阵要求恰好 n-1 行
//      （它总是把所有东西并成一棵树），所以这里可以显式要求补"桥接"行：
//      把这些根按给定距离依次合并。默认不补，`n_components > 1` 时由调用方决定。
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
	Reciprocal,   // 1/sim（sim <= 0 时 +inf）
};

// scipy linkage 矩阵的一行
struct ScipyRow {
	std::uint32_t left;      // 较小的子节点
	std::uint32_t right;     // 较大的子节点
	double distance;
	std::uint32_t size;
};

struct ScipyOptions {
	DistanceMap map = DistanceMap::Complement;
	// Complement 用的 C。0 表示"自动取树里的最大相似度"。
	double reference = 0.0;
	// 是否把森林补成一棵树（scipy 要求恰好 n-1 行）
	bool bridge = false;
	// 桥接用的距离。NaN = 自动，取当前映射下最忠实的值：
	//   Complement -> reference（C）。因为"两边没有边"在所有 pair 上的距离都是 C，
	//                 这也正是 scipy 自己对这两个分量算出来的值。
	//   Negate / Reciprocal -> +inf（没有比这更远的说法了）。
	double bridge_distance = std::numeric_limits<double>::quiet_NaN();
};

struct ScipyExport {
	std::vector<ScipyRow> rows;
	std::size_t n_objects = 0;
	std::size_t n_components = 0;   // 树里有几个根；>1 就是森林
	bool bridged = false;           // 是否补过桥接行
	// 距离是否沿合并序列非递减。scipy 的 dendrogram / fcluster 都假定这一点；
	// 三种链接里 average 与 max 保证单调，MinNonzero 不保证。
	bool monotonic = true;
};

// 纯函数：只读树，不改任何状态。
ScipyExport to_scipy(const std::vector<Merge>& tree, std::size_t n_objects,
					 const ScipyOptions& options = ScipyOptions());

// 该映射在给定参考值下的距离函数（导出时也可用于自己算阈值）
double map_distance(double similarity, DistanceMap map, double reference);

} // namespace spclust

#endif

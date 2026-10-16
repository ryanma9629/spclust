// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

#ifndef SPCLUST_HAC_H
#define SPCLUST_HAC_H

// ===========================================================================
// 基于稀疏相似度矩阵的凝聚型层次聚类。
//
// 用之前请先读完这段契约：
//
//   * 输入是**稀疏**相似度矩阵：没有给出的 pair 视为相似度 0，不占内存。
//   * 对象 id 由调用方给，必须是**稠密的 0..n-1**（内部按 max id 分配数组）。
//     如果手上是任意/稀疏的 id，在调用方做一次映射即可（tools/spclust_cli.cpp
//     里对字符串名做的就是这件事，约十行）。
//   * 同一个 pair 多次出现时**取最大值**。
//   * 自环（a == b）会破坏内部"行号 < 列号"的前提，抛 std::invalid_argument。
//   * 相似度必须是有限非负数；负数 / NaN 抛 std::invalid_argument。
//     （没有上界 —— 历史上限定 65535 是 unsigned short 时代的遗留。）
//   * **并列（相同相似度）时的合并顺序不作保证**：同一个链接方式内部是确定的，
//     但不同链接方式之间、不同版本之间可能不同。MaxSimilarity 下据此切出来的分组
//     与顺序无关（就是连通分量）；Average / MinNonzero 下可能不同。
//   * 一个 HAC 对象是一份独立状态，不同对象互不影响；但**同一个对象不是线程安全的**
//     （读接口也不与写并发安全）。
//
// 与 scipy 的对应关系：tree() 就是 linkage 矩阵，一行是
//     left right similarity size
// 节点编号与 scipy 一致（0..n-1 是对象，n+i 是第 i 次合并产生的簇）。
// **注意 height 在这里是相似度（越大越先合并），scipy 用的是距离。**
// ===========================================================================

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace spclust {

// 链接方式。命名按**相似度**语义（与 sparseAHC 的文档一致）：
//
//   Average       —— 两个簇之间**全部 |A|×|B| 对**相似度的平均（缺失边算 0）。
//                    分母是全部对数，所以簇越大越容易被"稀释"。
//   MaxSimilarity —— 取最大。按距离语义就是"单链接"：只要有任意一条边就把两簇连起来，
//                    会产生 chaining。它的簇间相似度只可能是**原始边**的值，
//                    所以内部走 Kruskal（不建矩阵），在大稀疏图上快得多。
//   MinNonzero    —— 取最小，但**仅当两边都有边**时才更新（缺边的那一对不参与）。
//                    即"全链接"的一个保守变体：按标准术语命名会叫 complete，
//                    但在稀疏图上把缺失边当 0 会让 min 恒为 0、失去意义，
//                    所以这里显式用 MinNonzero 这个名字。
enum class Linkage { Average, MaxSimilarity, MinNonzero };

// 一次合并 —— 标准 linkage 矩阵的一行。
struct Merge {
	std::uint32_t left;      // 左子节点
	std::uint32_t right;     // 右子节点
	double similarity;       // 合并时的相似度（该簇的"高度"，越大越先合并）
	double cohesion;         // 合并后整个簇的组内平均相似度
	std::size_t size;        // 合并后的簇大小
};

// 切树的结果：每个簇的成员（对象 id）与其组内平均相似度。
struct Cut {
	std::vector<std::vector<std::uint32_t>> members;
	std::vector<double> cohesion;
};

namespace detail {
struct HACImpl;   // 实现在 src/hac.cpp，这样头文件不暴露任何内部结构
}

// 一次层次聚类。一个对象持有一份独立状态，可以同时存在多个。
class HAC {
public:
	explicit HAC(Linkage linkage = Linkage::Average);
	~HAC();

	// 不可拷贝（持有实现状态）；可以移动
	HAC(const HAC&) = delete;
	HAC& operator=(const HAC&) = delete;
	HAC(HAC&&) noexcept;
	HAC& operator=(HAC&&) noexcept;

	// 预留容量。可选，但流式喂入大量边时能省掉中途扩容。
	void reserve(std::size_t n_objects, std::size_t n_edges);

	// 喂一条边。缺失的 pair 视为相似度 0；同一个 pair 多次出现取最大值。
	void add_edge(std::uint32_t a, std::uint32_t b, double similarity);

	// 算合并树。
	// floor 是"树算到多低"的相似度下限：合并到最大相似度 < floor 就停。
	// floor <= 0 表示算完整棵树（一直合到没有正相似度的候选为止）。
	// 注意阈值**不是**停止条件式的"聚类阈值"，而是树的深度：切在别处见 cut()。
	void build(double floor = 0.0);

	// 对象数（= 树的叶子数）。
	std::size_t num_objects() const;

	// 合并树。build() 之后有效，按合并先后排列。
	const std::vector<Merge>& tree() const;

	// 把树切在相似度 min_similarity 处，得到扁平分组。
	// 相似度沿合并序列单调不增，所以"切"就是取 similarity >= min_similarity 的那段前缀。
	Cut cut(double min_similarity) const;

	// 要正好 k 个簇：等价于按相似度从高到低合并到只剩 k 个。
	// k >= num_objects() 时返回全部分散；k == 0 当作 1。
	Cut cut_n_clusters(std::size_t k) const;

	// 只要"大小 >= min_size"的簇：取相似度最高的那段前缀，直到再合并会出现
	// 小于 min_size 的中间簇为止（即不切出比 min_size 更小的簇）。
	Cut cut_min_size(std::size_t min_size) const;

	// 清空全部状态，便于在同一对象上接着跑下一次聚类。
	void reset();

	Linkage linkage() const;
	const char* linkage_name() const;

private:
	std::unique_ptr<detail::HACImpl> impl_;
};

} // namespace spclust

#endif

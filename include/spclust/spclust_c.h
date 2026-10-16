/* SPDX-License-Identifier: MIT
 * Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
 *
 * spclust 的 C 接口 —— 给 FFI 用（Rust / Go / Python ctypes / C 项目）。
 *
 * 这个头文件是**纯 C**：只用 <stddef.h> / <stdint.h>，能被 C 编译器直接 include。
 * 用法与语义和 C++ 的 spclust::HAC 完全一致，契约见 include/spclust/hac.h，
 * 这里只重复要点：
 *   * 相似度矩阵是稀疏的，没给出的 pair 视为相似度 0；
 *   * 对象 id 必须是稠密的 0..n-1；
 *   * 同一个 pair 多次出现取最大值；
 *   * 自环、负数、NaN 会报错（返回非 0 错误码，不会崩）；
 *   * 并列时合并顺序不作保证；
 *   * 一个 handle 是一份独立状态，不同 handle 互不影响；同一个 handle 不是线程安全的。
 *
 * 错误处理：C++ 异常一律在库内被抓住，转成下面的错误码 —— 不会有异常穿过 C 边界。
 * 出错的细节可以用 spclust_last_error() 取（线程局部，下次调用前有效）。
 *
 * 典型用法：
 *
 *     spclust_hac* h = spclust_hac_new(SPCLUST_LINKAGE_AVERAGE);
 *     if (!h) { fprintf(stderr, "%s\n", spclust_last_error()); return 1; }
 *     if (spclust_hac_add_edge(h, 0, 1, 90.0) != SPCLUST_OK) { ... }
 *     if (spclust_hac_build(h, 0.0) != SPCLUST_OK) { ... }
 *
 *     size_t n = spclust_hac_tree_size(h);          // 先问长度
 *     spclust_merge* rows = malloc(n * sizeof *rows); // 自己分配缓冲
 *     if (spclust_hac_tree(h, rows, n, &n) != SPCLUST_OK) { ... }
 *     for (size_t i = 0; i < n; i++)
 *         printf("%u %u %g %zu\n", rows[i].left, rows[i].right,
 *                rows[i].similarity, rows[i].size);
 *     free(rows);
 *
 *     spclust_cut* c = spclust_hac_cut(h, 50.0);    // 返回值归调用方，用 spclust_cut_free 释放
 *     if (!c) { ... }
 *     for (size_t g = 0; g < spclust_cut_n_groups(c); g++) {
 *         const size_t* off = spclust_cut_offsets(c);
 *         const uint32_t* mem = spclust_cut_members(c);
 *         for (size_t k = off[g]; k < off[g + 1]; k++) printf(" %u", mem[k]);
 *     }
 *     spclust_cut_free(c);
 *     spclust_hac_free(h);
 */

#ifndef SPCLUST_C_H
#define SPCLUST_C_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ 版本 */

#define SPCLUST_VERSION_MAJOR 0
#define SPCLUST_VERSION_MINOR 3
#define SPCLUST_VERSION_PATCH 2

/* 运行时的版本字符串，形如 "0.3.2" */
const char* spclust_version(void);

/* ------------------------------------------------------------------ 错误码 */

typedef enum {
	SPCLUST_OK = 0,
	SPCLUST_ERR_INVALID_ARGUMENT = 1,  /* 自环 / 负数 / NaN / 非法 linkage 编号 ... */
	SPCLUST_ERR_OUT_OF_MEMORY = 2,     /* 分配失败（C++ 的 bad_alloc） */
	SPCLUST_ERR_STATE = 3,             /* 调用顺序不对（比如还没 build 就取 tree） */
	SPCLUST_ERR_BUFFER_TOO_SMALL = 4,  /* 给的输出缓冲不够 */
	SPCLUST_ERR_UNKNOWN = 5,           /* 其他 C++ 异常 */
} spclust_status;

/* 最近的错误描述。**线程局部**：每个线程各存一份，一直有效到该线程下次调用本库。
 * 没有错误时返回 NULL。 */
const char* spclust_last_error(void);

/* ------------------------------------------------------------------ 链接方式 */

typedef enum {
	SPCLUST_LINKAGE_AVERAGE = 0,         /* 全部 |A|x|B| 对的平均（缺失边算 0） */
	SPCLUST_LINKAGE_MAX_SIMILARITY = 1,  /* 取最大（距离语义下的单链接，走 Kruskal） */
	SPCLUST_LINKAGE_MIN_NONZERO = 2,     /* 取最小，但仅当两边都有边 */
} spclust_linkage;

/* 名字（"average" / "max-similarity" / "min-nonzero"）-> 编号，认不出来返回 -1 */
int spclust_linkage_from_name(const char* name);

/* ------------------------------------------------------------------ 句柄 */

typedef struct spclust_hac spclust_hac;

/* 建一个聚类器。linkage 见上面的枚举。失败返回 NULL（见 spclust_last_error）。 */
spclust_hac* spclust_hac_new(int linkage);

/* 释放。传 NULL 是安全的。 */
void spclust_hac_free(spclust_hac* h);

/* 预留容量（可选）。 */
spclust_status spclust_hac_reserve(spclust_hac* h, size_t n_objects, size_t n_edges);

/* 喂一条边。 */
spclust_status spclust_hac_add_edge(spclust_hac* h, uint32_t a, uint32_t b,
									double similarity);

/* 算合并树。floor <= 0 表示算完整棵树。 */
spclust_status spclust_hac_build(spclust_hac* h, double floor);

/* 清空状态，同一 handle 可以接着跑下一次。 */
spclust_status spclust_hac_reset(spclust_hac* h);

/* 对象数（= 树的叶子数）。 */
size_t spclust_hac_num_objects(const spclust_hac* h);

/* 合并树的行数。 */
size_t spclust_hac_tree_size(const spclust_hac* h);

/* 当前链接方式的编号 / 名字。 */
int spclust_hac_linkage(const spclust_hac* h);
const char* spclust_hac_linkage_name(const spclust_hac* h);

/* ------------------------------------------------------------------ linkage 矩阵 */

typedef struct {
	uint32_t left;
	uint32_t right;
	double similarity;
	double cohesion;
	size_t size;
} spclust_merge;

/* 取合并树。out 可以是 NULL（只问行数，*count 被写成实际行数）；
 * 否则 out 至少能放 cap 行，*count 被写成实际行数。
 * 缓冲不够时返回 SPCLUST_ERR_BUFFER_TOO_SMALL（*count 仍是需要的行数）。 */
spclust_status spclust_hac_tree(const spclust_hac* h, spclust_merge* out,
								size_t cap, size_t* count);

/* ------------------------------------------------------------------ 切树 */

/* 不透明结果。里面的内存归库所有，用 spclust_cut_free 释放整个结果。 */
typedef struct spclust_cut spclust_cut;

/* 切在相似度 min_similarity 处 */
spclust_cut* spclust_hac_cut(const spclust_hac* h, double min_similarity);
/* 要正好 k 个簇 */
spclust_cut* spclust_hac_cut_n_clusters(const spclust_hac* h, size_t k);
/* 合并到没有小于 min_size 的簇为止 */
spclust_cut* spclust_hac_cut_min_size(const spclust_hac* h, size_t min_size);

void spclust_cut_free(spclust_cut* c);

size_t spclust_cut_n_groups(const spclust_cut* c);
size_t spclust_cut_n_members(const spclust_cut* c);
/* 长度 n_members 的对象 id，按组连续存放 */
const uint32_t* spclust_cut_members(const spclust_cut* c);
/* 长度 n_groups+1；第 g 组的成员是 members[offsets[g] .. offsets[g+1]) */
const size_t* spclust_cut_offsets(const spclust_cut* c);
/* 长度 n_groups；每组的组内平均相似度 */
const double* spclust_cut_cohesion(const spclust_cut* c);

/* ------------------------------------------------------------------ scipy 导出 */

typedef enum {
	SPCLUST_DISTANCE_COMPLEMENT = 0,  /* C - sim（默认，语义最忠实） */
	SPCLUST_DISTANCE_NEGATE = 1,      /* -sim */
	SPCLUST_DISTANCE_RECIPROCAL = 2,  /* 1/sim */
} spclust_distance_map;

typedef struct {
	uint32_t left;
	uint32_t right;
	double distance;
	uint32_t size;
} spclust_scipy_row;

/* 导出成 scipy.cluster.hierarchy 的 linkage 矩阵（height 是**距离**）。
 * 语义与 include/spclust/scipy.h 一致：见那里的说明。
 * bridge != 0 时把森林补成一棵树（scipy 要求恰好 n-1 行），桥接距离为 NaN 时自动选。 */
spclust_status spclust_hac_to_scipy(const spclust_hac* h, int distance_map,
									double reference, int bridge, double bridge_distance,
									spclust_scipy_row* out, size_t cap,
									size_t* count, size_t* n_components);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SPCLUST_C_H */

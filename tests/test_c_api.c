/* SPDX-License-Identifier: MIT
 * Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
 *
 * C API 测试。这个文件是**真 C**（用 C 编译器编），所以它同时验证两件事：
 *   1. include/spclust/spclust_c.h 能被 C 编译器直接 include；
 *   2. extern "C" 的符号名与 ABI 在 C 侧可用。
 *
 * 编译：见 Makefile 的 `make test-c`（test_c_api.o 用 $(CC)，链接用 $(CXX)）。
 */

#include "spclust/spclust_c.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_pass = 0;
static int g_fail = 0;

static void check(int ok, const char* what) {
	if (ok) {
		g_pass++;
	} else {
		g_fail++;
		printf("  X %s\n", what);
	}
}

/* 6 个对象、3 个簇 {0,1} {2,3} {4,5}，全部值互不相同（不引入并列） */
static const struct { uint32_t a, b; double s; } EDGES[] = {
	{0, 1, 91.0}, {2, 3, 88.0}, {4, 5, 85.0},
	{0, 2, 31.0}, {1, 3, 22.0}, {0, 4, 13.0}, {2, 5, 7.0}, {1, 5, 3.0},
};
static const size_t N_EDGES = sizeof(EDGES) / sizeof(EDGES[0]);
#define N_OBJECTS 6

static spclust_hac* build_small(int linkage) {
	spclust_hac* h = spclust_hac_new(linkage);
	if (!h)
		return NULL;
	for (size_t i = 0; i < N_EDGES; i++)
		if (spclust_hac_add_edge(h, EDGES[i].a, EDGES[i].b, EDGES[i].s) != SPCLUST_OK)
			return NULL;
	if (spclust_hac_build(h, 0.0) != SPCLUST_OK)
		return NULL;
	return h;
}

static void test_version_and_names(void) {
	check(spclust_version() != NULL && strlen(spclust_version()) > 0, "版本字符串非空");
	check(strcmp(spclust_version(), "0.3.2") == 0, "版本号与头文件里的宏一致");
	check(SPCLUST_VERSION_MAJOR == 0 && SPCLUST_VERSION_MINOR == 3,
		  "头文件里的版本宏存在");
	check(spclust_linkage_from_name("average") == SPCLUST_LINKAGE_AVERAGE, "名字 average");
	check(spclust_linkage_from_name("max-similarity") == SPCLUST_LINKAGE_MAX_SIMILARITY,
		  "名字 max-similarity");
	check(spclust_linkage_from_name("single") == SPCLUST_LINKAGE_MAX_SIMILARITY, "别名 single");
	check(spclust_linkage_from_name("min") == SPCLUST_LINKAGE_MIN_NONZERO, "别名 min");
	check(spclust_linkage_from_name("nope") == -1, "认不出来的名字返回 -1");
	check(spclust_linkage_from_name(NULL) == -1, "NULL 名字返回 -1");
}

static void test_error_paths(void) {
	check(spclust_hac_new(7) == NULL, "非法 linkage -> NULL");
	check(spclust_last_error() != NULL, "非法 linkage 之后有错误描述");

	/* NULL 句柄必须安全（返回错误码而不是崩） */
	check(spclust_hac_add_edge(NULL, 0, 1, 1.0) == SPCLUST_ERR_INVALID_ARGUMENT,
		  "add_edge(NULL) 返回参数错误");
	check(spclust_hac_build(NULL, 0.0) == SPCLUST_ERR_INVALID_ARGUMENT, "build(NULL)");
	check(spclust_hac_cut(NULL, 1.0) == NULL, "cut(NULL) 返回 NULL");
	check(spclust_hac_num_objects(NULL) == 0, "num_objects(NULL) 返回 0");
	spclust_hac_free(NULL);        /* 不该崩 */
	spclust_cut_free(NULL);

	spclust_hac* h = spclust_hac_new(SPCLUST_LINKAGE_AVERAGE);
	check(h != NULL, "建 handle");
	check(spclust_hac_add_edge(h, 3, 3, 5.0) == SPCLUST_ERR_INVALID_ARGUMENT,
		  "自环 -> 参数错误");
	check(spclust_hac_add_edge(h, 0, 1, -1.0) == SPCLUST_ERR_INVALID_ARGUMENT,
		  "负相似度 -> 参数错误");
	check(spclust_hac_add_edge(h, 0, 1, NAN) == SPCLUST_ERR_INVALID_ARGUMENT,
		  "NaN -> 参数错误");
	check(strlen(spclust_last_error()) > 0, "错误描述非空");

	/* 出错之后再喂合法数据应当恢复 */
	check(spclust_hac_add_edge(h, 0, 1, 5.0) == SPCLUST_OK, "出错后仍可用");
	check(spclust_hac_linkage(h) == SPCLUST_LINKAGE_AVERAGE, "linkage 编号");
	check(strcmp(spclust_hac_linkage_name(h), "average") == 0, "linkage 名字");
	spclust_hac_free(h);
}

static void test_tree(void) {
	spclust_hac* h = build_small(SPCLUST_LINKAGE_AVERAGE);
	check(h != NULL, "建好小用例");
	if (!h)
		return;

	check(spclust_hac_num_objects(h) == N_OBJECTS, "对象数");
	check(spclust_hac_tree_size(h) == N_OBJECTS - 1, "完整树有 n-1 行");

	/* 先问长度（out = NULL） */
	size_t n = 0;
	check(spclust_hac_tree(h, NULL, 0, &n) == SPCLUST_OK, "只问行数");
	check(n == N_OBJECTS - 1, "问出来的行数对");

	/* 缓冲不够：返回 BUFFER_TOO_SMALL，*n 仍是需要的行数 */
	spclust_merge* rows = (spclust_merge*)malloc(n * sizeof(spclust_merge));
	check(rows != NULL, "分配输出缓冲");
	if (!rows) { spclust_hac_free(h); return; }
	check(spclust_hac_tree(h, rows, n - 1, &n) == SPCLUST_ERR_BUFFER_TOO_SMALL,
		  "缓冲差一行 -> BUFFER_TOO_SMALL");
	check(n == N_OBJECTS - 1, "缓冲不够时 *n 仍是需要的行数");

	check(spclust_hac_tree(h, rows, n, &n) == SPCLUST_OK, "取整棵树");
	/* 结构不变量：子节点 id 必须小于父节点 id（父节点 = n + 行号） */
	int ok_ids = 1, ok_height = 1;
	size_t last_size = 0;
	for (size_t i = 0; i < n; i++) {
		uint32_t parent = (uint32_t)(N_OBJECTS + i);
		if (rows[i].left >= parent || rows[i].right >= parent)
			ok_ids = 0;
		if (i + 1 < n && rows[i + 1].similarity > rows[i].similarity)
			ok_height = 0;
		last_size = rows[i].size;
	}
	check(ok_ids, "子节点 id < 父节点 id");
	check(rows[0].size == 2, "第一次合并的 size = 2");
	check(last_size == N_OBJECTS, "最后一次合并的 size = n");
	check(ok_height, "相似度沿合并序列单调不增");
	check(rows[0].similarity < 92.0 && rows[0].similarity > 90.0, "第一次合并的相似度");

	/* 三个簇的先合并：91 / 88 / 85 */
	check(fabs(rows[0].similarity - 91.0) < 1e-9 &&
		  fabs(rows[1].similarity - 88.0) < 1e-9 &&
		  fabs(rows[2].similarity - 85.0) < 1e-9, "前三次合并是组内那三对");

	free(rows);
	spclust_hac_free(h);
}

static void test_cuts(void) {
	spclust_hac* h = build_small(SPCLUST_LINKAGE_AVERAGE);
	if (!h)
		return;

	/* cut(50)：组内三对都 > 50，跨组的都 < 50 -> 正好 3 个簇 */
	spclust_cut* c = spclust_hac_cut(h, 50.0);
	check(c != NULL, "cut(50)");
	if (c) {
		check(spclust_cut_n_groups(c) == 3, "cut(50) 得到 3 个簇");
		check(spclust_cut_n_members(c) == N_OBJECTS, "成员总数 = n");
		const size_t* off = spclust_cut_offsets(c);
		const uint32_t* mem = spclust_cut_members(c);
		const double* coh = spclust_cut_cohesion(c);
		check(off != NULL && mem != NULL && coh != NULL, "三个视图非 NULL");
		char seen[N_OBJECTS];
		memset(seen, 0, sizeof seen);
		int ok = (off[0] == 0);
		for (size_t g = 0; g < spclust_cut_n_groups(c); g++) {
			if (off[g + 1] - off[g] != 2)
				ok = 0;                          /* 每组 2 个 */
			if (fabs(coh[g] - (g == 0 ? 91.0 : (g == 1 ? 88.0 : 85.0))) > 1e-9)
				ok = 0;
			for (size_t k = off[g]; k < off[g + 1]; k++) {
				if (mem[k] >= N_OBJECTS || seen[mem[k]])
					ok = 0;
				seen[mem[k]] = 1;
			}
		}
		check(ok, "每组 2 个成员、紧致度正确、成员是 0..n-1 的一个排列");
		check(off[spclust_cut_n_groups(c)] == N_OBJECTS, "offsets 最后一项 = 成员总数");
		spclust_cut_free(c);
	}

	/* cut_n_clusters */
	spclust_cut* k2 = spclust_hac_cut_n_clusters(h, 2);
	check(k2 && spclust_cut_n_groups(k2) == 2, "cut_n_clusters(2) 得到 2 个簇");
	spclust_cut_free(k2);
	spclust_cut* k1 = spclust_hac_cut_n_clusters(h, 1);
	check(k1 && spclust_cut_n_groups(k1) == 1 &&
		  spclust_cut_n_members(k1) == N_OBJECTS, "cut_n_clusters(1) 只有一个簇");
	spclust_cut_free(k1);

	/* cut_min_size */
	spclust_cut* m2 = spclust_hac_cut_min_size(h, 2);
	check(m2 && spclust_cut_n_groups(m2) == 3, "cut_min_size(2) 得到 3 个 size>=2 的簇");
	spclust_cut_free(m2);
	spclust_cut* m1 = spclust_hac_cut_min_size(h, 1);
	check(m1 && spclust_cut_n_groups(m1) == N_OBJECTS, "cut_min_size(1) = 全部分散");
	spclust_cut_free(m1);

	/* 极高阈值 -> 全部分散 */
	spclust_cut* hi = spclust_hac_cut(h, 1000.0);
	check(hi && spclust_cut_n_groups(hi) == N_OBJECTS, "cut(1000) = 全部分散");
	spclust_cut_free(hi);

	spclust_hac_free(h);
}

static void test_reset_and_instances(void) {
	spclust_hac* h = build_small(SPCLUST_LINKAGE_AVERAGE);
	if (!h)
		return;
	check(spclust_hac_reset(h) == SPCLUST_OK, "reset");
	check(spclust_hac_num_objects(h) == 0, "reset 之后对象数归零");
	check(spclust_hac_tree_size(h) == 0, "reset 之后树为空");

	/* 两个 handle 同时存在、互不影响 */
	spclust_hac* a = build_small(SPCLUST_LINKAGE_AVERAGE);
	spclust_hac* b = build_small(SPCLUST_LINKAGE_MIN_NONZERO);
	check(a && b, "两个 handle 同时建好");
	if (a && b) {
		check(spclust_hac_tree_size(a) == N_OBJECTS - 1, "a 的树不受 b 影响");
		check(spclust_hac_linkage(a) == SPCLUST_LINKAGE_AVERAGE &&
			  spclust_hac_linkage(b) == SPCLUST_LINKAGE_MIN_NONZERO, "两个 handle 的 linkage 独立");
		spclust_hac_free(b);
		check(spclust_hac_tree_size(a) == N_OBJECTS - 1, "b 释放后 a 仍然可用");
	}
	spclust_hac_free(a);
	spclust_hac_free(h);
}

static void test_scipy_export(void) {
	spclust_hac* h = build_small(SPCLUST_LINKAGE_AVERAGE);
	if (!h)
		return;
	size_t n = 0, comps = 0;
	check(spclust_hac_to_scipy(h, SPCLUST_DISTANCE_COMPLEMENT, 0.0, 0, NAN,
							   NULL, 0, &n, &comps) == SPCLUST_OK, "scipy 导出：问行数");
	check(n == N_OBJECTS - 1, "scipy 行数 = n-1");
	check(comps == 1, "完整树的连通分量数 = 1");

	spclust_scipy_row* rows = (spclust_scipy_row*)malloc(n * sizeof(spclust_scipy_row));
	if (rows) {
		check(spclust_hac_to_scipy(h, SPCLUST_DISTANCE_COMPLEMENT, 0.0, 0, NAN,
								   rows, n, &n, &comps) == SPCLUST_OK, "scipy 导出：取行");
		/* complement 映射：dist = C - sim，C 自动取 91 */
		check(fabs(rows[0].distance - 0.0) < 1e-9, "最大相似度 -> 距离 0");
		check(rows[0].left < rows[0].right, "子节点升序");
		int mono = 1;
		for (size_t i = 1; i < n; i++)
			if (rows[i].distance < rows[i - 1].distance)
				mono = 0;
		check(mono, "距离单调");
		free(rows);
	}
	check(spclust_hac_to_scipy(h, 99, 0.0, 0, NAN, NULL, 0, &n, &comps) ==
		  SPCLUST_ERR_INVALID_ARGUMENT, "非法距离映射 -> 参数错误");
	spclust_hac_free(h);
}

int main(void) {
	printf("C API 测试\n");
	printf("--------------------------------------------------------------------\n");
	test_version_and_names();
	test_error_paths();
	test_tree();
	test_cuts();
	test_reset_and_instances();
	test_scipy_export();
	printf("--------------------------------------------------------------------\n");
	printf("PASS=%d FAIL=%d\n", g_pass, g_fail);
	return g_fail == 0 ? 0 : 1;
}

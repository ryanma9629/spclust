/* SPDX-License-Identifier: MIT
 * Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
 *
 * C API test. This file is **real C** (compiled by a C compiler), so it verifies two
 * things at once:
 *   1. include/spclust/spclust_c.h can be included directly by a C compiler;
 *   2. the extern "C" symbol names and the ABI work from C.
 *
 * Build: see `make test-c` in the Makefile (test_c_api.o with $(CC), linked with $(CXX)).
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

/* 6 objects, 3 clusters {0,1} {2,3} {4,5}, all values distinct (no ties introduced) */
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
	check(spclust_version() != NULL && strlen(spclust_version()) > 0, "the version string is non-empty");
	check(strcmp(spclust_version(), "0.3.3") == 0, "the version matches the header macros");
	check(SPCLUST_VERSION_MAJOR == 0 && SPCLUST_VERSION_MINOR == 3,
		  "the version macros exist in the header");
	check(spclust_linkage_from_name("average") == SPCLUST_LINKAGE_AVERAGE, "the name average");
	check(spclust_linkage_from_name("max-similarity") == SPCLUST_LINKAGE_MAX_SIMILARITY,
		  "the name max-similarity");
	check(spclust_linkage_from_name("single") == SPCLUST_LINKAGE_MAX_SIMILARITY, "the alias single");
	check(spclust_linkage_from_name("min") == SPCLUST_LINKAGE_MIN_NONZERO, "the alias min");
	check(spclust_linkage_from_name("nope") == -1, "an unrecognized name returns -1");
	check(spclust_linkage_from_name(NULL) == -1, "a NULL name returns -1");
}

static void test_error_paths(void) {
	check(spclust_hac_new(7) == NULL, "bad linkage -> NULL");
	check(spclust_last_error() != NULL, "there is an error description after a bad linkage");

	/* A NULL handle must be safe (an error code, never a crash) */
	check(spclust_hac_add_edge(NULL, 0, 1, 1.0) == SPCLUST_ERR_INVALID_ARGUMENT,
		  "add_edge(NULL) returns an argument error");
	check(spclust_hac_build(NULL, 0.0) == SPCLUST_ERR_INVALID_ARGUMENT, "build(NULL)");
	check(spclust_hac_cut(NULL, 1.0) == NULL, "cut(NULL) returns NULL");
	check(spclust_hac_num_objects(NULL) == 0, "num_objects(NULL) returns 0");
	spclust_hac_free(NULL);        /* must not crash */
	spclust_cut_free(NULL);

	spclust_hac* h = spclust_hac_new(SPCLUST_LINKAGE_AVERAGE);
	check(h != NULL, "create a handle");
	check(spclust_hac_add_edge(h, 3, 3, 5.0) == SPCLUST_ERR_INVALID_ARGUMENT,
		  "self-loop -> argument error");
	check(spclust_hac_add_edge(h, 0, 1, -1.0) == SPCLUST_ERR_INVALID_ARGUMENT,
		  "negative similarity -> argument error");
	check(spclust_hac_add_edge(h, 0, 1, NAN) == SPCLUST_ERR_INVALID_ARGUMENT,
		  "NaN -> argument error");
	check(strlen(spclust_last_error()) > 0, "the error description is non-empty");

	/* After an error, valid input must still work */
	check(spclust_hac_add_edge(h, 0, 1, 5.0) == SPCLUST_OK, "still usable after an error");
	check(spclust_hac_linkage(h) == SPCLUST_LINKAGE_AVERAGE, "the linkage id");
	check(strcmp(spclust_hac_linkage_name(h), "average") == 0, "the linkage name");
	spclust_hac_free(h);
}

static void test_add_edges_batch(void) {
	/* The batch entry must be exactly equivalent to the per-edge loop. */
	uint32_t a[N_EDGES], b[N_EDGES];
	double s[N_EDGES];
	for (size_t i = 0; i < N_EDGES; i++) {
		a[i] = EDGES[i].a;
		b[i] = EDGES[i].b;
		s[i] = EDGES[i].s;
	}

	spclust_hac* one = spclust_hac_new(SPCLUST_LINKAGE_AVERAGE);
	spclust_hac* batch = spclust_hac_new(SPCLUST_LINKAGE_AVERAGE);
	check(one && batch, "two handles for the batch comparison");
	if (!one || !batch)
		return;
	for (size_t i = 0; i < N_EDGES; i++)
		spclust_hac_add_edge(one, a[i], b[i], s[i]);
	check(spclust_hac_add_edges(batch, a, b, s, N_EDGES) == SPCLUST_OK, "add_edges returns OK");
	check(spclust_hac_build(one, 0.0) == SPCLUST_OK &&
		  spclust_hac_build(batch, 0.0) == SPCLUST_OK, "both build");

	size_t n1 = spclust_hac_tree_size(one), n2 = spclust_hac_tree_size(batch);
	check(n1 == n2 && n1 == N_OBJECTS - 1, "same number of rows");
	if (n1 == n2) {
		spclust_merge* r1 = malloc(n1 * sizeof *r1);
		spclust_merge* r2 = malloc(n2 * sizeof *r2);
		if (r1 && r2 && spclust_hac_tree(one, r1, n1, &n1) == SPCLUST_OK &&
			spclust_hac_tree(batch, r2, n2, &n2) == SPCLUST_OK) {
			int same = 1;
			for (size_t i = 0; i < n1; i++)
				if (r1[i].left != r2[i].left || r1[i].right != r2[i].right ||
					r1[i].similarity != r2[i].similarity || r1[i].size != r2[i].size ||
					r1[i].cohesion != r2[i].cohesion)
					same = 0;
			check(same, "batch and per-edge trees are identical");
		}
		free(r1);
		free(r2);
	}
	/* count == 0 (arrays may be NULL) and the argument checks */
	check(spclust_hac_add_edges(batch, NULL, NULL, NULL, 0) == SPCLUST_OK, "count 0 is a no-op");
	check(spclust_hac_add_edges(batch, NULL, b, s, 3) == SPCLUST_ERR_INVALID_ARGUMENT,
		  "NULL array with count > 0 -> argument error");
	check(spclust_hac_add_edges(NULL, a, b, s, 3) == SPCLUST_ERR_INVALID_ARGUMENT,
		  "NULL handle -> argument error");
	/* an invalid edge inside a batch is reported */
	spclust_hac* bad = spclust_hac_new(SPCLUST_LINKAGE_AVERAGE);
	double wrong[N_EDGES];
	for (size_t i = 0; i < N_EDGES; i++) wrong[i] = EDGES[i].s;
	wrong[3] = -1.0;
	check(spclust_hac_add_edges(bad, a, b, wrong, N_EDGES) == SPCLUST_ERR_INVALID_ARGUMENT,
		  "a negative similarity inside a batch -> argument error");
	check(spclust_last_error() != NULL, "and it has a message");
	spclust_hac_free(bad);
	spclust_hac_free(one);
	spclust_hac_free(batch);
}

static void test_tree(void) {
	spclust_hac* h = build_small(SPCLUST_LINKAGE_AVERAGE);
	check(h != NULL, "build the small case");
	if (!h)
		return;

	check(spclust_hac_num_objects(h) == N_OBJECTS, "object count");
	check(spclust_hac_tree_size(h) == N_OBJECTS - 1, "a complete tree has n-1 rows");

	/* Ask for the length first (out = NULL) */
	size_t n = 0;
	check(spclust_hac_tree(h, NULL, 0, &n) == SPCLUST_OK, "ask for the row count only");
	check(n == N_OBJECTS - 1, "the reported row count is right");

	/* Buffer too small: BUFFER_TOO_SMALL, and *n still holds the required count */
	spclust_merge* rows = (spclust_merge*)malloc(n * sizeof(spclust_merge));
	check(rows != NULL, "allocate the output buffer");
	if (!rows) { spclust_hac_free(h); return; }
	check(spclust_hac_tree(h, rows, n - 1, &n) == SPCLUST_ERR_BUFFER_TOO_SMALL,
		  "one row short -> BUFFER_TOO_SMALL");
	check(n == N_OBJECTS - 1, "*n still holds the required count when the buffer is short");

	check(spclust_hac_tree(h, rows, n, &n) == SPCLUST_OK, "fetch the whole tree");
	/* Structural invariant: child ids must be below the parent id (parent = n + row index) */
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
	check(ok_ids, "child id < parent id");
	check(rows[0].size == 2, "the first merge has size = 2");
	check(last_size == N_OBJECTS, "the last merge has size = n");
	check(ok_height, "similarity is non-increasing along the merge sequence");
	check(rows[0].similarity < 92.0 && rows[0].similarity > 90.0, "the similarity of the first merge");

	/* The three intra-cluster merges come first: 91 / 88 / 85 */
	check(fabs(rows[0].similarity - 91.0) < 1e-9 &&
		  fabs(rows[1].similarity - 88.0) < 1e-9 &&
		  fabs(rows[2].similarity - 85.0) < 1e-9, "the first three merges are the intra-cluster pairs");

	free(rows);
	spclust_hac_free(h);
}

static void test_cuts(void) {
	spclust_hac* h = build_small(SPCLUST_LINKAGE_AVERAGE);
	if (!h)
		return;

	/* cut(50): the three intra-cluster pairs are > 50 and every cross pair is < 50 -> 3 clusters */
	spclust_cut* c = spclust_hac_cut(h, 50.0);
	check(c != NULL, "cut(50)");
	if (c) {
		check(spclust_cut_n_groups(c) == 3, "cut(50) gives 3 clusters");
		check(spclust_cut_n_members(c) == N_OBJECTS, "the member count is n");
		const size_t* off = spclust_cut_offsets(c);
		const uint32_t* mem = spclust_cut_members(c);
		const double* coh = spclust_cut_cohesion(c);
		check(off != NULL && mem != NULL && coh != NULL, "all three views are non-NULL");
		char seen[N_OBJECTS];
		memset(seen, 0, sizeof seen);
		int ok = (off[0] == 0);
		for (size_t g = 0; g < spclust_cut_n_groups(c); g++) {
			if (off[g + 1] - off[g] != 2)
				ok = 0;                          /* every cluster has 2 members */
			if (fabs(coh[g] - (g == 0 ? 91.0 : (g == 1 ? 88.0 : 85.0))) > 1e-9)
				ok = 0;
			for (size_t k = off[g]; k < off[g + 1]; k++) {
				if (mem[k] >= N_OBJECTS || seen[mem[k]])
					ok = 0;
				seen[mem[k]] = 1;
			}
		}
		check(ok, "2 members per cluster, correct cohesion, members are a permutation of 0..n-1");
		check(off[spclust_cut_n_groups(c)] == N_OBJECTS, "the last offset equals the member count");
		spclust_cut_free(c);
	}

	/* cut_n_clusters */
	spclust_cut* k2 = spclust_hac_cut_n_clusters(h, 2);
	check(k2 && spclust_cut_n_groups(k2) == 2, "cut_n_clusters(2) gives 2 clusters");
	spclust_cut_free(k2);
	spclust_cut* k1 = spclust_hac_cut_n_clusters(h, 1);
	check(k1 && spclust_cut_n_groups(k1) == 1 &&
		  spclust_cut_n_members(k1) == N_OBJECTS, "cut_n_clusters(1) gives a single cluster");
	spclust_cut_free(k1);

	/* cut_min_size */
	spclust_cut* m2 = spclust_hac_cut_min_size(h, 2);
	check(m2 && spclust_cut_n_groups(m2) == 3, "cut_min_size(2) gives 3 clusters of size >= 2");
	spclust_cut_free(m2);
	spclust_cut* m1 = spclust_hac_cut_min_size(h, 1);
	check(m1 && spclust_cut_n_groups(m1) == N_OBJECTS, "cut_min_size(1) = everything separate");
	spclust_cut_free(m1);

	/* A very high threshold -> everything separate */
	spclust_cut* hi = spclust_hac_cut(h, 1000.0);
	check(hi && spclust_cut_n_groups(hi) == N_OBJECTS, "cut(1000) = everything separate");
	spclust_cut_free(hi);

	spclust_hac_free(h);
}

static void test_reset_and_instances(void) {
	spclust_hac* h = build_small(SPCLUST_LINKAGE_AVERAGE);
	if (!h)
		return;
	check(spclust_hac_reset(h) == SPCLUST_OK, "reset");
	check(spclust_hac_num_objects(h) == 0, "the object count is zero after reset");
	check(spclust_hac_tree_size(h) == 0, "the tree is empty after reset");

	/* Two handles at once must not affect each other */
	spclust_hac* a = build_small(SPCLUST_LINKAGE_AVERAGE);
	spclust_hac* b = build_small(SPCLUST_LINKAGE_MIN_NONZERO);
	check(a && b, "both handles are created");
	if (a && b) {
		check(spclust_hac_tree_size(a) == N_OBJECTS - 1, "a's tree is unaffected by b");
		check(spclust_hac_linkage(a) == SPCLUST_LINKAGE_AVERAGE &&
			  spclust_hac_linkage(b) == SPCLUST_LINKAGE_MIN_NONZERO, "the two handles have independent linkages");
		spclust_hac_free(b);
		check(spclust_hac_tree_size(a) == N_OBJECTS - 1, "a is still usable after b is freed");
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
							   NULL, 0, &n, &comps) == SPCLUST_OK, "scipy export: ask for the row count");
	check(n == N_OBJECTS - 1, "scipy row count = n-1");
	check(comps == 1, "a complete tree has 1 connected component");

	spclust_scipy_row* rows = (spclust_scipy_row*)malloc(n * sizeof(spclust_scipy_row));
	if (rows) {
		check(spclust_hac_to_scipy(h, SPCLUST_DISTANCE_COMPLEMENT, 0.0, 0, NAN,
								   rows, n, &n, &comps) == SPCLUST_OK, "scipy export: fetch the rows");
		/* complement mapping: dist = C - sim, with C taken from the tree (91) */
		check(fabs(rows[0].distance - 0.0) < 1e-9, "the largest similarity maps to distance 0");
		check(rows[0].left < rows[0].right, "children are ascending");
		int mono = 1;
		for (size_t i = 1; i < n; i++)
			if (rows[i].distance < rows[i - 1].distance)
				mono = 0;
		check(mono, "distances are monotone");
		free(rows);
	}
	check(spclust_hac_to_scipy(h, 99, 0.0, 0, NAN, NULL, 0, &n, &comps) ==
		  SPCLUST_ERR_INVALID_ARGUMENT, "bad distance mapping -> argument error");
	spclust_hac_free(h);
}

int main(void) {
	printf("C API test\n");
	printf("--------------------------------------------------------------------\n");
	test_version_and_names();
	test_error_paths();
	test_add_edges_batch();
	test_tree();
	test_cuts();
	test_reset_and_instances();
	test_scipy_export();
	printf("--------------------------------------------------------------------\n");
	printf("PASS=%d FAIL=%d\n", g_pass, g_fail);
	return g_fail == 0 ? 0 : 1;
}

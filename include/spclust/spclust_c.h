/* SPDX-License-Identifier: MIT
 * Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
 *
 * The C interface of spclust - for FFI (Rust / Go / Python ctypes / C projects).
 *
 * This header is **plain C**: it only uses <stddef.h> and <stdint.h> and can be
 * included directly by a C compiler. The API and its semantics are identical to
 * the C++ spclust::HAC; the contract lives in include/spclust/hac.h, and the
 * essentials are:
 *   * the similarity matrix is sparse, an absent pair counts as similarity 0;
 *   * object ids must be dense 0..n-1;
 *   * if a pair appears more than once the largest value wins;
 *   * self-loops, negative values and NaN are reported as errors (a non-zero
 *     status, never a crash);
 *   * ties are broken by id, not by insertion order: among pairs tied on similarity the
 *     **smaller (left, right) id pair** wins, so the same edges with the same ids give a
 *     byte-identical tree in any order;
 *   * one handle holds one independent state, separate handles never affect each
 *     other; a single handle is not thread-safe.
 *
 * Error handling: C++ exceptions are always caught inside the library and turned
 * into the status codes below, so no exception ever crosses the C boundary.
 * Details are available from spclust_last_error() (thread-local; valid until the
 * next call from the same thread).
 *
 * Typical use:
 *
 *     spclust_hac* h = spclust_hac_new(SPCLUST_LINKAGE_AVERAGE);
 *     if (!h) { fprintf(stderr, "%s\n", spclust_last_error()); return 1; }
 *     if (spclust_hac_add_edge(h, 0, 1, 90.0) != SPCLUST_OK) { ... }
 *     if (spclust_hac_build(h, 0.0) != SPCLUST_OK) { ... }
 *
 *     size_t n = spclust_hac_tree_size(h);              // ask for the length first
 *     spclust_merge* rows = malloc(n * sizeof *rows);   // allocate your own buffer
 *     if (spclust_hac_tree(h, rows, n, &n) != SPCLUST_OK) { ... }
 *     for (size_t i = 0; i < n; i++)
 *         printf("%u %u %g %zu\n", rows[i].left, rows[i].right,
 *                rows[i].similarity, rows[i].size);
 *     free(rows);
 *
 *     spclust_cut* c = spclust_hac_cut(h, 50.0);  // owned by the caller; free with spclust_cut_free
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

/* ------------------------------------------------------------------ version */

#define SPCLUST_VERSION_MAJOR 0
#define SPCLUST_VERSION_MINOR 3
#define SPCLUST_VERSION_PATCH 3

/* Runtime version string, e.g. "0.3.2" */
const char* spclust_version(void);

/* ------------------------------------------------------------------ status */

typedef enum {
	SPCLUST_OK = 0,
	SPCLUST_ERR_INVALID_ARGUMENT = 1,  /* self-loop / negative / NaN / bad linkage id ... */
	SPCLUST_ERR_OUT_OF_MEMORY = 2,     /* allocation failed (C++ bad_alloc) */
	SPCLUST_ERR_STATE = 3,             /* calls out of order (e.g. tree() before build()) */
	SPCLUST_ERR_BUFFER_TOO_SMALL = 4,  /* the output buffer is too small */
	SPCLUST_ERR_UNKNOWN = 5,           /* some other C++ exception */
} spclust_status;

/* Description of the last error. **Thread-local**: each thread has its own, and
 * it stays valid until the next call into this library from that thread.
 * Returns NULL when there is no error. */
const char* spclust_last_error(void);

/* ------------------------------------------------------------------ linkage */

typedef enum {
	SPCLUST_LINKAGE_AVERAGE = 0,         /* mean over all |A|x|B| pairs (missing edges count as 0) */
	SPCLUST_LINKAGE_MAX_SIMILARITY = 1,  /* maximum (single linkage in distance terms; Kruskal path) */
	SPCLUST_LINKAGE_MIN_NONZERO = 2,     /* minimum, but only over pairs where both sides have an edge */
} spclust_linkage;

/* Name ("average" / "max-similarity" / "min-nonzero") -> id; -1 if unrecognized */
int spclust_linkage_from_name(const char* name);

/* ------------------------------------------------------------------ handle */

typedef struct spclust_hac spclust_hac;

/* Create a clusterer. linkage is one of the enum values above.
 * Returns NULL on failure (see spclust_last_error). */
spclust_hac* spclust_hac_new(int linkage);

/* Release. Passing NULL is safe. */
void spclust_hac_free(spclust_hac* h);

/* Reserve capacity (optional). */
spclust_status spclust_hac_reserve(spclust_hac* h, size_t n_objects, size_t n_edges);

/* Feed one edge. */
spclust_status spclust_hac_add_edge(spclust_hac* h, uint32_t a, uint32_t b,
									double similarity);

/* Feed `count` edges in one call: a[i], b[i], sim[i] for i in 0..count-1.
 * Equivalent to calling spclust_hac_add_edge for each i, but one foreign call instead of
 * count of them - which is what makes this usable from another language at all, since a
 * cross-language call costs microseconds while the work per edge is nanoseconds.
 * On the first invalid edge it returns that edge's status; **edges before it have already
 * been added** (this is not transactional). To locate the offender, scan the arrays with
 * the same rules (self-loop, negative, NaN).
 * `a`, `b` and `sim` may be NULL when count == 0. */
spclust_status spclust_hac_add_edges(spclust_hac* h, const uint32_t* a, const uint32_t* b,
									 const double* sim, size_t count);

/* Compute the merge tree. floor <= 0 computes the whole tree. */
spclust_status spclust_hac_build(spclust_hac* h, double floor);

/* Clear all state so the same handle can be reused. */
spclust_status spclust_hac_reset(spclust_hac* h);

/* Number of objects (= number of leaves). */
size_t spclust_hac_num_objects(const spclust_hac* h);

/* Number of rows in the merge tree. */
size_t spclust_hac_tree_size(const spclust_hac* h);

/* Linkage id / name of this handle. */
int spclust_hac_linkage(const spclust_hac* h);
const char* spclust_hac_linkage_name(const spclust_hac* h);

/* ------------------------------------------------------------------ linkage matrix */

typedef struct {
	uint32_t left;
	uint32_t right;
	double similarity;
	double cohesion;
	size_t size;
} spclust_merge;

/* Fetch the merge tree. out may be NULL (then only the row count is written to
 * *count); otherwise out must hold at least cap rows and *count receives the
 * actual row count. If the buffer is too small the function returns
 * SPCLUST_ERR_BUFFER_TOO_SMALL (*count still holds the required row count). */
spclust_status spclust_hac_tree(const spclust_hac* h, spclust_merge* out,
								size_t cap, size_t* count);

/* ------------------------------------------------------------------ cutting */

/* Opaque result. Its memory is owned by the library; free the whole result with
 * spclust_cut_free. */
typedef struct spclust_cut spclust_cut;

/* Cut at similarity min_similarity */
spclust_cut* spclust_hac_cut(const spclust_hac* h, double min_similarity);
/* Exactly k clusters */
spclust_cut* spclust_hac_cut_n_clusters(const spclust_hac* h, size_t k);
/* Merge until no cluster is smaller than min_size */
spclust_cut* spclust_hac_cut_min_size(const spclust_hac* h, size_t min_size);

void spclust_cut_free(spclust_cut* c);

size_t spclust_cut_n_groups(const spclust_cut* c);
size_t spclust_cut_n_members(const spclust_cut* c);
/* Object ids, n_members of them, stored contiguously group by group */
const uint32_t* spclust_cut_members(const spclust_cut* c);
/* Length n_groups+1; group g is members[offsets[g] .. offsets[g+1]) */
const size_t* spclust_cut_offsets(const spclust_cut* c);
/* Length n_groups; the mean intra-cluster similarity of each group */
const double* spclust_cut_cohesion(const spclust_cut* c);

/* ------------------------------------------------------------------ scipy export */

typedef enum {
	SPCLUST_DISTANCE_COMPLEMENT = 0,  /* C - sim (default; the most faithful mapping) */
	SPCLUST_DISTANCE_NEGATE = 1,      /* -sim */
	SPCLUST_DISTANCE_RECIPROCAL = 2,  /* 1/sim */
} spclust_distance_map;

typedef struct {
	uint32_t left;
	uint32_t right;
	double distance;
	uint32_t size;
} spclust_scipy_row;

/* Export as a scipy.cluster.hierarchy linkage matrix (height is a **distance**).
 * Semantics match include/spclust/scipy.h; see there for details.
 * When bridge != 0 a forest is turned into a single tree (scipy requires exactly
 * n-1 rows); bridge_distance == NaN picks the value automatically. */
spclust_status spclust_hac_to_scipy(const spclust_hac* h, int distance_map,
									double reference, int bridge, double bridge_distance,
									spclust_scipy_row* out, size_t cap,
									size_t* count, size_t* n_components);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SPCLUST_C_H */

# Design and algorithm

This document covers the internal structures, the algorithm and a few semantic conventions.
For the interface see the [README](../README.md), for history and open issues
[known-issues](known-issues.md), and for measurements [performance](performance.md).

## Data structures

The input is an **edge list** rather than a matrix: `(obj1, obj2, sim)`, where an absent pair
counts as similarity 0 and is not stored.

```
vec_dim1 (the row vector)         vec_dim2 (the column vector)
 row i: { head -> [cell chain] , max -> this row's largest cell }
 col j: { head -> [cell chain] }
```

- A `sp_unit`'s four links (`left/right/up/down`) are **indices**, not pointers, and the row
  and column chains share the same cells. All cells live in one contiguous block, with
  `sizeof(sp_unit)` at 32 bytes.
- **Only the upper triangle is stored**: `add_grp()` enforces `group1 < group2` (the smaller
  becomes the row, the larger the column).
- Every row keeps a `max` pointer. This is the key optimisation: `get_max()` only scans the
  column of per-row maxima, so it is O(number of groups) rather than O(nnz).

## Algorithm

```
build(floor)                        src/hac.cpp
  add_edge(a, b, sim)               register a singleton group per new id; insert the edge
  while (groups > 1):
      unit = get_max()              the globally largest cell (heap top, amortised O(1))
      if (no cell || sim < floor) break
      merge_grp(group1, group2, sim)
```

`merge_grp()` has four steps:

1. Walk group1's and group2's row and column chains, collecting the **row/column indices that
   need recomputing into one reused small array** (rather than scanning the whole index space)
   and marking the rows whose `max` pointer is invalidated. This is where sparsity pays off.
2. Compute the new group's similarity to each neighbour according to the linkage, appending the
   new group id at the end (`grp.size()`). The two similarities recorded along the way are also
   stored in `marks`, so the next step needs no chain lookups.
3. Remove and release group1's and group2's entire row and column with four-way pointer
   splicing, releasing the centre cell last.
4. Rescan the invalidated rows to repair their `max` pointers, and record this merge in
   `merges` (which is the linkage matrix).

Single linkage (`Linkage::MaxSimilarity`) takes a different route: its inter-cluster
similarity is the maximum over all **original edges** between two clusters, so the whole tree
is equivalent to "sort the edges by similarity descending, then union-find" (Kruskal), never
building the matrix at all. See [performance](performance.md).

## Linkage methods

There are three, given when constructing a `HAC`. The names follow **similarity** semantics,
and note that the comments in the original code were the wrong way round (it called max
"complete" and min "single"); behaviour is what counts here:

| `Linkage` | Inter-cluster similarity | In distance terms | Notes |
|---|---|---|---|
| `Average` | the mean over all `|A|x|B|` pairs (missing edges count as 0) | average | larger clusters get diluted |
| `MaxSimilarity` | the maximum | **single** | it chains; its values can only be original edges, so it uses Kruskal |
| `MinNonzero` | the minimum, but only where both sides have an edge | a conservative complete | missing edges do not participate; otherwise min would be identically 0 |

`MinNonzero` is **not** standard complete linkage: on a sparse matrix treating missing edges
as 0 makes the minimum identically 0 and meaningless, so only pairs where **both sides have an
edge** participate - hence the explicit name. One consequence is that it can leave
**components that never merge** (two clusters with no edge between them are not even a
candidate), so the tree may be incomplete; `cut_n_clusters` / `cut_min_size` then return the
finest cut the tree supports. Its intra-cluster tightness is also initialised to the +inf of
the empty set (`empty_cohesion`), because 0 is the absorbing element of min and would make the
report identically 0.

## State and instances

This used to be **34 file-static variables plus a pile of static functions**: one clustering
per process, and changing the linkage meant recompiling. All of that state now lives in
`spclust::detail::HACImpl`, and the public header exposes a single `class HAC` (holding a
`unique_ptr`):

- **One object, one state**: a process can hold any number of clusterings, interleaving edges
  and `build()` calls, and the results are identical to running each alone (tested explicitly
  in `tests/test_cut.cpp`, including reset reuse).
- **No mutable globals in the library**: `.data`/`.bss` of `libspclust.a` according to `nm`
  contains only C++ runtime constants (string literals, jump tables), no algorithm state. So
  multi-threading is safe as long as no single object is shared.
- **The linkage is chosen at run time**: `HAC(Linkage)`. Dispatch happens only at the public
  entry points (`add_edge` / `build` / `cut` / ...), while the inner loops stay in one function
  where the linkage is known, leaving no branch in the hot path - a measured cost of 0.2-0.6%
  (alternating A/B, min/quartile/median agreeing).
- **The merge workspace is released when `build()` returns**: the node pool and the O(n) buffers
  (row/column vectors, marks, cached similarities, heap, group table) are handed back rather than
  kept at capacity, so a process holding the result owns roughly the tree - about 32 bytes per
  object - instead of several times that. At n = 1M that is 30.5 MB instead of 163.5 MB measured in
  allocated bytes. Re-building on the same object re-allocates the workspace, measured to cost
  nothing outside noise.
- **Cutting is a pure function**: the union-find parent array is a local of `cut()`. It used to
  be a file static, which made `cut()` nominally `const` while mutating global state and let
  concurrent cuts corrupt each other.

## Tie breaking

When several candidate pairs share the highest similarity, the winner is the one with the
**smaller (left, right) id pair** - first the smaller left id, then the smaller right id. The
same rule holds for all three linkages, and it is implemented in three places:

- `entry_less`, the heap comparator: on equal similarity the smaller **row** (= smaller left id)
  is the larger heap key;
- `add_grp`, when a cell is inserted and has to displace the row's current max;
- the max-repair loop, which rescans a row whose `max` was invalidated.

All three go through one helper, `beats(sim_a, col_a, sim_b, col_b)`. Two of the three used to
take the *first* candidate they met (`>` instead of `>=` plus an id comparison), which meant the
chain order - i.e. the order edges were added in - decided the outcome. Measured at the time:
four to five of six input permutations produced a different tree. With ids deciding, the tree is
a function of the edge set alone, which is what `tests/test_cut.cpp` now asserts.

## Cutting

The tree is kept in full (the `floor` argument of `build()` only decides how far down to
compute it), so a threshold is only "where to cut" and changing the cut never re-clusters.
All three entry points reduce to "take the first `use` merges":

| Entry point | Semantics | Tree needed |
|---|---|---|
| `cut(t)` | every merge with `similarity >= t` (ties come along) | only down to `t` |
| `cut_n_clusters(k)` | exactly k clusters = the first `n-k` merges | deep enough; otherwise the finest available |
| `cut_min_size(m)` | the **shallowest** cut where every cluster is >= m | deep enough |

`cut_min_size` is not "discard small clusters" but "keep merging until there are none". Note
that under all three semantics **two clusters can have no edge at all** (under `MinNonzero`
such a pair is not even a candidate), so disconnected components may never merge - in which
case it returns the finest cut the tree supports.

## Tunables

**The library itself has no compile-time options** (the linkage is a constructor argument).
What remains belongs to the command-line tool, in `tools/cli_config.h`, wrapped in `#ifndef`
so `-D` on the command line wins:

```bash
make MIN_SIM=60 DEBUG=2          # the default cut threshold and the debug level
```

- `MIN_SIM_CLUST`: the default of `--cut` (50)
- `MAX_SIM_READ`: the input bound on similarity (65535). **This is historical** (it used to be
  an `unsigned short`); the library accepts any finite non-negative `double`
- `DEBUG`: 0 = silent, 1 = phases and progress, 2 = per-merge detail
- `DEBUG_FREQ_READDB`: how often read progress is printed

The linkage, the cut mode and `--floor` are all run-time arguments; see `spclust_cli --help`.

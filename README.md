# spclust

Agglomerative hierarchical clustering over a **sparse similarity matrix**.

It repeatedly merges the pair of clusters with the highest global similarity until the
largest similarity drops below a given lower bound. A **cross-linked list** stores the
sparse upper triangle, so memory is O(nnz) rather than O(n²) - the main advantage over a
dense implementation. The merge tree is kept in full, so a threshold is only "where to cut
the tree" and changing it never requires re-clustering.

```
library:  include/spclust/hac.h + src/hac.cpp     no third-party dependencies, integer ids only
user:     tools/                                  all IO and presentation (CSV, id <-> name, printing)
```

## Quick start

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure      # differential + CLI contract tests
./build/spclust_cli sim.csv                     # run the sample data
```

Standard library only, no third-party dependencies: a C++17 compiler is enough. A
Makefile is also available as a development shortcut:

```bash
make            # build ./spclust
make lib        # build just the static library, libspclust.a
make run        # run sim.csv
make test       # differential + cut property + CLI + scipy interop + C API
make test-cut   # cut property test only
make test-scipy # scipy interop only (skips itself when scipy is absent)
make test-c     # C API: compile with a C compiler, then call it from Python via ctypes
make lib-shared # libspclust.so, for FFI
make bench      # performance regression guard (fixed inputs + a recorded baseline)
make check      # run sim.csv once per linkage
make asan       # AddressSanitizer + UBSan
make asan-cases # run every case in tests/cases under ASan/UBSan
make leaks      # check for leaks
make clean
```

## Using it as a library

```cmake
find_package(spclust REQUIRED)
target_link_libraries(your_target PRIVATE spclust::spclust)
```

```cpp
#include <spclust/hac.h>

spclust::HAC hac(spclust::Linkage::Average);   // the linkage is a **run-time** argument

hac.reserve(n_objects, n_edges);               // optional: reserve capacity
for (auto& e : edges)
    hac.add_edge(e.a, e.b, e.similarity);      // ids must be dense 0..n-1

hac.build(0.0);                                // 0 = compute the whole tree

for (const spclust::Merge& m : hac.tree())
    ...                                        // this is the linkage matrix: left right similarity size

// Three ways to cut, switchable at will, none of them re-clusters
spclust::Cut c = hac.cut(50.0);                // cut at similarity 50
spclust::Cut k = hac.cut_n_clusters(12);       // exactly 12 clusters
spclust::Cut s = hac.cut_min_size(5);          // merge until no cluster is smaller than 5

for (const auto& members : c.members)
    ...                                        // object ids per cluster (Cut also carries cohesion)

hac.reset();                                   // reuse the same object for the next clustering

// Any number of independent clusterings can coexist in one process, with no locking.
// (Just never use the **same** HAC object from two threads.)
spclust::HAC other(spclust::Linkage::MaxSimilarity);
```

**Read the contract at the top of `include/spclust/hac.h` before using it** (an absent edge
counts as 0, ids must be dense, repeated pairs take the maximum, self-loops throw, the merge
order of tied similarities is not specified, one object holds one independent state that is
not safe to read concurrently with a write). How the linkage shapes the clusters is in
[docs/design.md](docs/design.md#linkage-methods).

## Input format

One `obj1,obj2,sim` per line, or read from standard input (`spclust_cli < data.csv`, or
`spclust_cli -`).

- `sim`: the CLI accepts a real number in 0..65535, fractions included (e.g. `87.5`). **The
  upper bound is the CLI's input contract** (a historical artefact); the library itself
  accepts any finite non-negative `double`. Edges below the cut threshold are only
  candidates and never trigger a merge.
- Separators: comma when the line contains one, otherwise whitespace (spaces / tabs)
- A `#` starts a comment; blank lines are ignored; an optional `obj1,obj2,sim` header is
  skipped automatically
- An out-of-range or non-numeric `sim`, too few fields, or a self-loop (`obj1 == obj2`)
  prints a warning and skips that line
- **Repeated pairs take the maximum** (the semantics are defined and pinned by tests)

## Output

By default it prints the flat partition obtained by cutting the merge tree at `--cut`
(50 by default):

```
Objects: 16, groups: 3
Group 0 (size=5, sim=90.5): a1 a2 a3 a4 a5
Group 1 (size=4, sim=87.6667): b1 b2 b3 b4
Group 2 (size=3, sim=82): c1 c2 c3
Unclustered (4): d1 d2 e1 e2
```

`sim` is the mean intra-cluster similarity (under average linkage it accumulates weighted by
cluster size) and can be read as a tightness measure. Only clusters of size >= 2 are listed;
the remaining singletons appear under `Unclustered`, members in id order (first-appearance
order in the input).

Because the whole merge tree is computed, a threshold is only "where to cut":

```
spclust_cli data.csv --cut 86         # cut at 86 (no re-clustering)
spclust_cli data.csv --n-clusters 12  # exactly 12 clusters (needs the full tree)
spclust_cli data.csv --min-size 5     # merge until no cluster is smaller than 5 (needs the full tree)
spclust_cli data.csv --linkage 1      # linkage: 0/average 1/max|single 2/min|min-nonzero
spclust_cli data.csv --tree           # dendrogram (--tree-depth N limits depth, 0 = unlimited)
spclust_cli data.csv --matrix         # linkage matrix: left right height size
spclust_cli data.csv --floor 0        # compute the tree all the way down (full tree; default = cut, cheaper)
```

Node ids follow scipy's convention: `0..n-1` are the objects and `n+i` is merge i - here the
group number is naturally the node id. **Note that `height` is a similarity (higher merges
first), while scipy uses a distance.** The semantics and cost of `--cut` / `--floor` are in
[docs/design.md](docs/design.md).

## C interface (for FFI)

`include/spclust/spclust_c.h` is a **plain C** header (only `<stddef.h>` / `<stdint.h>`), so
C, Rust, Go or Python ctypes can include or bind it directly:

```c
#include "spclust/spclust_c.h"

spclust_hac* h = spclust_hac_new(SPCLUST_LINKAGE_AVERAGE);
spclust_hac_add_edge(h, 0, 1, 90.0);
spclust_hac_build(h, 0.0);

size_t n = spclust_hac_tree_size(h);
spclust_merge* rows = malloc(n * sizeof *rows);
spclust_hac_tree(h, rows, n, &n);          /* two calls: ask for the length, then fill the buffer */

spclust_cut* c = spclust_hac_cut(h, 50.0);
for (size_t g = 0; g < spclust_cut_n_groups(c); g++) { ... }
spclust_cut_free(c);
spclust_hac_free(h);
```

Design points:

- **C++ exceptions never cross the C boundary.** Every entry point `try/catch`es and returns a
  `spclust_status` (`SPCLUST_ERR_INVALID_ARGUMENT` / `OUT_OF_MEMORY` / `BUFFER_TOO_SMALL` /
  ...). Details come from `spclust_last_error()`, which is **thread-local**.
  (`tests/test_c_api.py` tests exactly this: a self-loop passed from Python must yield an
  error code rather than aborting the process.)
- **Opaque handles** (`spclust_hac*` / `spclust_cut*`); one handle is one state, and separate
  handles never affect each other.
- **Ownership is explicit**: `spclust_hac_tree` writes into a buffer the caller allocated,
  while the flat views of a cut (`members` / `offsets` / `cohesion`) belong to the library and
  are released in one go by `spclust_cut_free`.
- The linkage, all three cut modes and the scipy export have entry points;
  `spclust_linkage_from_name` accepts `average` / `max-similarity` / `min-nonzero` (plus the
  `max`/`single`/`min` aliases).

```bash
make lib        # libspclust.a
make lib-shared # libspclust.so (for FFI; on macOS -shared also produces a Mach-O, which ctypes can dlopen)
make test-c     # compile the test with a C compiler, then call it from Python via ctypes
```

## scipy interop

```bash
spclust_cli data.csv --scipy --scipy-bridge > Z.txt
```
```python
from scipy.cluster.hierarchy import fcluster, dendrogram, to_tree
import numpy as np
Z = np.loadtxt("Z.txt")          # a scipy linkage matrix, ready to use
fcluster(Z, 3, criterion='maxclust')
```

In scipy `height` is a **distance** (smaller merges first) while this library works in
**similarities** (larger merges first), so the export has to choose a mapping
(`--scipy-distance`):

- `complement` (the default) `dist = C - sim`, with `C` the largest similarity. **This is the
  most faithful one**: this library's "an absent edge is similarity 0" corresponds exactly to
  "distance = C", so `average` linkage equals scipy's `method='average'` applied to the same
  dense distance matrix - on cases **without ties** it matches scipy row for row
  (`tests/test_scipy.py`, 23 combinations).
- `negate` `dist = -sim`: order-preserving and parameter-free, but distances are negative.
- `reciprocal` `dist = 1/sim`: only meaningful when similarities behave like affinities.

**Three differences worth knowing** (semantics rather than bugs; all pinned by
`tests/test_scipy.py`):

1. **Ties.** Equal similarities never guaranteed a merge order. That includes **derived
   values**: with 28 distinct input values, intermediate similarities still collide (this is
   what `15_many_merges.csv` does). The two trees then differ, and both are legal.
2. **`MinNonzero` is not scipy's `complete`.** scipy reaches the maximum distance
   (= similarity 0) on missing pairs, so on sparse input almost everything looks far apart,
   while this library takes the minimum only over pairs where both sides have an edge. **On
   dense input the two agree** (measured 6/6).
3. **Forests.** Under `MinNonzero`, two clusters with no edge between them are not even a
   candidate, so components may never merge and the tree is incomplete (fewer than n-1 rows).
   scipy requires exactly n-1 rows, so either add `--scipy-bridge` to append bridging rows
   (under `complement` the bridging distance is `C`, exactly what scipy itself would compute)
   or use another linkage. Without `--scipy-bridge` it warns on stderr.

## Documentation

| File | Contents |
|---|---|
| [docs/design.md](docs/design.md) | data structures, algorithm, linkage methods, tunables |
| [docs/known-issues.md](docs/known-issues.md) | bugs that were fixed (and how they were found) + open issues |
| [docs/performance.md](docs/performance.md) | measurements, optimisations that measurement rejected, methodology |
| [docs/testing.md](docs/testing.md) | cut properties / differential / CLI / scipy interop / C API (with FFI) / performance guard |
| [CHANGELOG.md](CHANGELOG.md) | changes that matter to users |

## Scale of the testing and performance work

- Differential test: an independent reference implementation (with no Lance-Williams
  recurrence) compared merge by merge, plus 22 hand-written boundary cases and 900
  multi-threshold cut cross-checks. Currently `PASS=90 WARN=12 KNOWN=0 FAIL=0` (a WARN is an
  unspecified merge order caused by a tie; both sides are legal).
- Cut property test (`tests/test_cut.cpp`, library level): the expectations are not a second
  copy of the union-find but come from **collecting leaf sets bottom-up from the merge tree** -
  a different code path. It covers three linkages x 21 random graphs x tied thresholds,
  exactly k clusters, minimum cluster size, cohesion provenance and monotonicity, plus
  multiple instances and reset reuse. Currently `PASS=18281 FAIL=0`.
- scipy interop (`tests/test_scipy.py`): the `--scipy` export is handed to the **real scipy**.
  23 cases x 2 linkages, all format-valid, and on the tie-free cases **row-for-row identical**
  to scipy's own `linkage`; tie cases are skipped and labelled as the contract allows.
- ASan/UBSan clean, 0 leaks; the node pool is poisoned by hand, so "reading a released node"
  is still caught by ASan.
- Cumulative performance (against the original baseline; see
  [docs/performance.md](docs/performance.md)): 800k edges **21.9x**, 2M edges **51.1x**, peak
  memory **-55 to -61%**, single linkage **300x**.

## License

MIT, see [LICENSE](LICENSE).

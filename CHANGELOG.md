# Changelog

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and version
numbers follow [Semantic Versioning](https://semver.org/).

This file records changes that matter to users. Per-change implementation detail (including
the measurements behind every performance change and the methodology corrections) lives in
the git history, where each commit is independently revertible.

## [0.3.3] - 2026-10

### Changed

- **Ties are now specified by id instead of being "unspecified".** Among pairs tied on the highest
  similarity the one with the **smaller (left, right) id pair** wins, for all three linkages. This
  makes the tree a function of the edge set alone: adding the same edges with the same ids in any
  order now produces a byte-identical tree. If you relied on the old wording - or if your
  reproducibility came from a fixed input order - note that the *rule* is new even though the
  *tree* is usually the same.
- `build()` now releases its workspace, so an object kept for cutting holds roughly the tree
  (~32 bytes per object) rather than several times that: 30.5 MB instead of 163.5 MB at n = 1M,
  measured in allocated bytes. Re-building on the same object re-allocates it, measured to cost
  nothing outside noise.

### Fixed

- Ties were actually resolved by **insertion order**, not by id: `add_grp`'s max insertion, the
  max-repair loop and `record_sim_g1/g2` all took the first candidate they met. Feeding the same
  edges in six shuffled orders produced a different tree in four to five of them (Kruskal was
  already order-independent). All three sites now share one `beats()` helper, and
  `tests/test_cut.cpp` asserts order independence as a property - verified to fail against the
  pre-fix code.
- `record_sim_g1/g2` recorded a derived similarity with "first write wins" while the documented
  semantic is that duplicates take the maximum. For duplicate pairs with different values the
  derived similarity therefore depended on insertion order. They now take the maximum, which is
  what the documentation always said.
- `free_mem()` used `clear()`, which keeps capacity, so the whole O(n) workspace stayed resident
  after `build()`.

### Notes

- Re-verified after this change: differential `PASS=102 WARN=0` (was `PASS=90 WARN=12`), cut
  property `PASS=18288 FAIL=0`, CLI 11/11, C API 64+31, scipy interop 46 combinations, ASan/UBSan
  clean, ctest 6/6, and the benchmark still matches the recorded baseline.

## [0.3.2] - 2026-10

### Added

- **C interface** (`include/spclust/spclust_c.h`): a plain C header for FFI (Rust / Go /
  Python ctypes / C projects). Opaque handles, two-call output (ask for the length, then
  fill the buffer), `spclust_status` error codes, and a thread-local
  `spclust_last_error()`. It covers the linkage, all three cut modes and the scipy export.
  **C++ exceptions never cross the C boundary**: every entry point try/catches into an error
  code, so bad input from another language yields a status, never an aborted host process.
- `tests/test_c_api.c` (compiled by a **C compiler**, which also proves the header is plain
  C) plus `tests/test_c_api.py` (ctypes against the shared library, checking dlopen, symbol
  names, struct layout and the exception boundary).
- `make test-c` and `make lib-shared` (`libspclust.so`).

## [0.3.1] - 2026-10

### Added

- **scipy interop** (`include/spclust/scipy.h` + `src/scipy.cpp`): `to_scipy()` exports the
  merge tree as a `scipy.cluster.hierarchy` linkage matrix; from the command line it is
  `spclust_cli data.csv --scipy [--scipy-distance complement|negate|reciprocal]
  [--scipy-bridge]`. The distance mapping has to be chosen by the caller (a similarity is not
  a distance); the default `complement` (`C - sim`) is the most faithful one, because
  "an absent edge is similarity 0" corresponds exactly to "distance = C", which makes average
  linkage equivalent to scipy's `method='average'` on the same dense distance matrix.
- `tests/test_scipy.py`: uses the **real scipy** as the ruler (every earlier test used an
  oracle we wrote ourselves). 23 cases x 2 linkages: all format-valid, and on the
  tie-free, complete-tree cases **row-for-row identical** to scipy's own `linkage`, heights
  and sizes included. Tie cases are skipped as the contract allows and labelled. Without
  numpy/scipy the test skips itself and returns 0.
- `make test-scipy`.

### One finding about ties (now documented)

Distinct input values can still tie: derived intermediate similarities collide.
`tests/cases/15_many_merges.csv` (28 distinct input values) has two candidates both equal to
82.5 at step 3; the oracle and we picked one branch and scipy picked the other, after which
all three trees differ - and all three are legal. That is what "the merge order of tied
similarities is unspecified" actually means.

### Known differences (not bugs; both pinned by tests)

- `Linkage::MinNonzero` is not scipy's `complete`: scipy reaches the maximum distance
  (= similarity 0) on missing pairs, so on sparse input almost everything looks far apart,
  while this library takes the minimum only over pairs where both sides have an edge. **On
  dense input the two agree** (measured 6/6).
- A forest (an incomplete tree, fewer than n-1 rows) is invisible to scipy's
  `is_valid_linkage`, which simply treats it as a smaller n. The CLI therefore warns on
  stderr, and `--scipy-bridge` appends the missing merges.

## [0.3.0] - 2026-10

This release **changed the public API** (free functions -> a `spclust::HAC` object), because
the previous interface had a fundamental limit: the state was 34 file-static variables, so
**one process could run exactly one clustering**, and changing the linkage meant recompiling
the whole library.

### Changed

- **`class spclust::HAC` replaces the free functions** (`reserve` / `add_edge` / `build` /
  `num_objects` / `tree` / `cut` / `reset` / `linkage_name` are now members). One object is
  one state: a process can hold any number of clusterings and interleaving them changes
  nothing. The library has no mutable globals left (no algorithm state in `.data`/`.bss`
  according to `nm libspclust.a`), so it is thread-safe as long as no single object is shared.
- **The linkage went from the compile-time `#define LINKAGE_MTD` to a constructor argument
  `spclust::Linkage`**: `Average` / `MaxSimilarity` / `MinNonzero`. Three binaries used to be
  needed to exercise three linkages; one now suffices (`spclust_cli --linkage 0|1|2`).
  Dispatch happens only at the public entry points and leaves no branch in the hot path -
  measured cost 0.2-0.6% (alternating A/B, min/quartile/median agreeing).
- **The linkage names were settled by behaviour.** The comments in the original code were
  wrong (it called max "complete" and min "single"); behaviour wins: max is single linkage
  (through Kruskal) and min takes the minimum only over pairs where both sides have an edge,
  which is why it is named `MinNonzero` rather than complete.
- `cut()` became a pure function: the union-find parent array is now a local. It used to be a
  file static, so `cut()` was nominally `const` while mutating global state and concurrent
  cuts corrupted each other.
- Command line: added `--linkage`, `--n-clusters <K>` and `--min-size <M>`. The latter two
  need the full tree, so without an explicit `--floor` they use `--floor 0` automatically.
- The CMake `SPCLUST_LINKAGE` option was removed (no longer needed), and
  `SPCLUST_BUILD_TOOLS` / `SPCLUST_BUILD_TESTS` now use `option()`, so `-D...=OFF` really
  takes effect (it used to be `set()`, which overwrote a user's OFF).

### Added

- `Cut cut_n_clusters(std::size_t k)`: exactly k clusters (= the first n-k merges).
- `Cut cut_min_size(std::size_t m)`: merge until no cluster is smaller than m (the
  **shallowest** cut satisfying that, not "discard small clusters").
- `tests/test_cut.cpp`: a library-level cut property test. Its expectations are derived
  independently by **collecting leaf sets bottom-up from the merge tree** (not by writing the
  union-find a second time), covering three linkages x 21 random graphs x tied thresholds,
  k clusters, minimum cluster size, cohesion provenance and monotonicity, plus **multiple
  instances and reset reuse**. Currently `PASS=18281 FAIL=0`.
- `make test-cut`.

### Fixed

- The linkage names and the flow diagram in `docs/design.md` (they still said `LINKAGE_MTD`
  and the old function names).

## [0.2.0] - 2026-10

Separating the library from its user (still with the free-function API).

### Added

- **Library and user split**: `include/spclust/hac.h` + `src/hac.cpp` are the library itself
  (no dependencies, integer ids only, prints nothing, reads no files); `tools/` is its first
  user and owns all IO and presentation.
- **CMake**: `find_package(spclust)` + `spclust::spclust` + `install` + `ctest`, with
  `SPCLUST_BUILD_TOOLS` / `SPCLUST_BUILD_TESTS` options. The `Makefile` stays as a
  development shortcut.
- Public API: `reserve` / `add_edge` / `build` / `num_objects` / `tree` / `cut` / `reset` /
  `linkage_name`. Cutting (`Cut`) moved into the library as well; it had been sitting in the
  CLI's anonymous namespace while being core functionality.
- `docs/`: four documents covering design, fixed/open issues, measurements and testing; the
  README only explains how to use the thing.
- GitHub Actions CI: two platforms x (build + ctest + ASan/UBSan over all cases + a
  library-only build).
- Full merge-tree output: `--cut <X>` (cut at similarity X), `--floor <X>` (how deep to
  compute), `--tree` / `--tree-depth <N>` (dendrogram), `--matrix` (the linkage matrix
  `left right height size`). The threshold is no longer a stopping condition but a cut
  position, so changing it never re-clusters.
- Performance guard: `make bench` (fixed inputs + recorded baseline + a machine-state canary)
  and `make bench-update`.
- CLI interface tests (`tests/test_cli.py`): exit-code contract, cut agreement with the
  oracle, `--matrix` structural invariants, `--tree` completeness, stdin equivalence.
- Differential test suite: `tests/oracle.py` (an independent implementation written from the
  definition, with no Lance-Williams recurrence) + `tests/run_tests.py` (merge-by-merge
  comparison) + 22 hand-written boundary cases.
- `LICENSE` (MIT) and an SPDX header on every source file.

### Changed

- The library no longer stores object names: `add_edge(a, b, sim)` takes integer ids and the
  id-to-name mapping belongs to the CLI. It also has no standard output any more (the
  per-merge debug log used to go to stdout).
- The library no longer calls `exit()`: allocation failure throws `std::bad_alloc` and invalid
  arguments throw `std::invalid_argument`.
- Similarities went from `unsigned short` to `double`. The original truncated with integer
  division at every step, so error accumulated downwards, collapsing distinct similarities
  into ties and **changing the final cluster count**.
- Single linkage (`LINKAGE_MTD == 1`) now goes through Kruskal: its inter-cluster similarity
  is the maximum over all **original edges**, so the whole tree equals "sort edges by
  similarity descending, then union-find" with no similarity matrix at all. Measured 63-300x
  depending on scale, with peak memory dropping from 15.8 GB to 50 MB (n=80000, 800k edges).
- Duplicate-pair semantics settled as **the maximum wins** (previously undefined behaviour).
- The reader became **streaming in small batches**: the whole edge list is no longer buffered
  in memory, cutting peak memory by 40-49%.
- Nodes shrank from 48 to 32 bytes (links and group numbers stored as `uint32` indices),
  cutting peak memory by a further ~18% while running 1.4-1.5x faster.
- `readdb` no longer calls `exit(1)` internally: it returns 0 on failure (with the reason on
  stderr) and lets the caller decide the exit code.
- Headers no longer contain `using namespace std;`, and `grp` is no longer a cross-file global.
- The cohesion report for `LINKAGE_MTD == 2` used to be identically 0 (the absorbing element
  of `min(...)`); it now reports the smallest intra-cluster edge.

### Fixed

- A heap-use-after-free in `get_max()`: when a whole row had similarity 0 the max pointer was
  never reassigned and kept pointing at a freed node. The default average linkage plus one
  input edge with `sim=0` is enough to trigger it. (The differential test reported PASS for
  that input - only ASan caught it.)
- The unused `MIN_ACCT_NUM` / `MIN_SIM_WG` macros were removed.

### Performance

Cumulative against the original baseline: 800k edges **21.9x**, 2M edges **51.1x**, peak
memory **-55 to -61%**, single linkage **300x**. The main contributors: removing four O(n)
sweeps per merge, turning `get_max` into a max-heap, hashing in `init`, eliminating
`get_unit`'s linear search, streaming input, 32-byte nodes, and Kruskal for single linkage.

## [0.1.0] - 2026-10

First usable state: the 2022 Windows/ODBC version made portable (CSV / standard input), plus
the differential test suite.

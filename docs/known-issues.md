# Fixed issues and known issues

The first half covers **bugs that were fixed** (with detail, because for several of them the way
they were found is the reason the test suite exists). The second half covers problems that are
**still there**. Per-change implementation detail lives in the git history.

## Fixed

### Average linkage stored similarities in an `unsigned short`, truncating at every step

The code before the fix:

```cpp
unsigned short sim = (sim1 * group1_size + sim2 * group2_size) / (group1_size + group2_size);
```

Integer division truncated and the result went back into an `unsigned short`, so **every merge
lost a little** - and what was lost went on to be used by the next step. `make test` reported 12
`TRUNC` results at the time:

```
rand_n44d80s1.csv (44 objects) step 27: oracle=41869.6  spclust=41868   (12 of 15 merges low)
rand_n10d30s3.csv: truncation collapsed 262.5 and 262.0 into one value, 262 -> a tie
                   -> the merge order changed -> oracle got 3 clusters, spclust got 4
```

The point is that this was not merely a display problem: **the criterion itself** (take the
globally largest similarity) was being applied to truncated values.

**Fix**: `sp_sim::sim` / `sp_unit::sim` / `grp_s::sim`, the related locals and the parameters of
`add_grp`/`merge_grp` all became `double`; `parse_sim` in `readdb.cpp` switched to `strtod`,
which incidentally started accepting fractions such as `87.5`. The `max`/`min` linkages only
compare and never average, so they were unaffected.

**A positive side effect**: `if (sim > 0)` used to turn a real but small average into 0 and drop
the candidate; now it does not, so the candidate set matches the definition more closely.

**Verification**: after the fix `TRUNC=0` and all 12 cases became `PASS`; across 146 per-merge
similarity comparisons the largest relative error is `4.3e-13` (pure print precision). `TRUNC` was
also reclassified as a **failure**, so if this bug ever came back `make test` would go red
immediately.

### A heap-use-after-free in `get_max()` (when a whole row had similarity 0)

The max-repair loop used `if (u->sim > max)` with `max` initialised to 0, so **a row whose
similarities were all 0 was never reassigned**. If that row's old max pointed into the column being
removed, `vec_dim1[i].max` kept a freed pointer and the next `get_max()` read it.

Trigger (`tests/cases/22_zero_row_max.csv`): a row's max is a similarity-0 cell in the removed
column while another 0-valued cell on that row **survives** (so the row does not become empty -
the empty case was already handled). The default average linkage plus one input edge with `sim=0`
is enough; `LINKAGE_MTD=1` is not needed.

Fix: scan the row for its true maximum and assign unconditionally, so the pointer is always either
a live cell or NULL. Tie behaviour is unchanged (still "the first cell on the chain attaining the
maximum"), and a 0-valued max can never be chosen by `get_max()` (which uses a strict `>`).

**How it was found is the interesting part**: the differential test reported **PASS** for that
input - the just-freed memory happened to read back as 0 and the merge order came out the same.
**Only ASan caught it.** Comparing results and checking memory safety are two independent jobs, and
both are needed.

### Three O(n) sweeps per merge (performance)

The old `merge_grp` allocated three O(n) mark vectors per merge and then did

```cpp
for (size_t i = 0; i < dim1_size; i++)   // dim1_size ~ 2n
    if (update_dim1[i]) { ... }          // only O(degree) are marked
```

i.e. it used an O(n) full sweep to filter out O(degree) indices. Together with `get_max()` and the
max-repair loop that was four full sweeps per merge. Measured (n=40000, 18293 merges):

```
update_dim1 loop   926,710,622 iterations
update_dim2 loop   927,123,378 iterations
max repair loop    926,7xx,xxx iterations
get_max()          926,769,374 iterations
------------------------------------------------
total              ~3.7e9 iterations, about 0.8 ns each
```

Fix: iterate only the marked indices, with the mark buffers reused across merges (resetting only
the positions touched). Every step became O(degree) instead of O(n). Two things came along for
free:

- **The two update loops had byte-identical bodies**, both doing `get_unit(x, group1)` and
  `get_unit(x, group2)` for the same index x. So when x appeared in both passes the old code
  inserted two **duplicate cells** with equal values for the same `(x, grp_size)`. Merging them
  into one pass removes those (at n=40000 mallocs dropped from 1,121,925 to 1,097,456, 2.2% fewer).
- No more O(n) allocation per merge (the three mark vectors).

**Measured 2.2x faster with byte-identical output** (28 inputs x 3 linkages, including the full
DEBUG=2 merge log). But **the complexity class did not change**, because `get_max()` was still
O(n) per merge, i.e. O(n²) - see the next item.

### `get_max()` scanning O(groups) per merge (performance)

After the previous change `get_max()` was **57-78%** of the clustering phase: it swept every row on
every merge, O(n²) in total. A max-heap makes it O(1).

A trap was hit here and is worth recording: the first version was an "index heap with in-place
updates" (one entry per row, keys read live from `vec_dim1[row].max`, sifting in place when a row's
max changed). **That is wrong.** One merge changes the keys of k rows, so the heap array's order no
longer matches the current keys - and `sift_up`/`sift_down` both assume "the heap is valid
everywhere except the one place being repaired", so sifting each position in place **cannot**
restore it. The differential test caught it: for `LINKAGE=2`, `rand_n44d25s1.csv` step 21, the
oracle took 167 while spclust took 63 (the heap top was not the global maximum).

The correct approach is to **freeze the key inside the heap entry** (`struct SimEntry { sim; row;
version; }`): the heap order is then always consistent with the entries and the heap is always
valid. When a row's max changes its version is bumped and a new entry pushed, with old entries
invalidated by version. The extra benefits:

- no bookkeeping of "where an entry sits in the heap", and no need for idempotent updates
- no dependence on call order: it can be maintained **in place** in `add_grp`, when removing
  group1/group2, and in the max rescan
- heap operations never touch `vec_dim1`, so the max pointers that dangle temporarily mid-merge
  are irrelevant to the heap
- when invalidated entries outnumber live rows by 4x the heap is rebuilt once (keeping only current
  entries plus Floyd heapify); the cost scales with the heap size, so it is O(1) amortised while
  memory stays O(live rows)

Tie behaviour matches the old sweep: compare similarity first, and on a tie the smaller row index
wins.

### `map<string,size_t>` string lookups in `init_grp` (performance)

The sampled profile showed **the init phase had grown to the same order as the clustering phase**,
with its cost concentrated in name lookups:

```
init_grp()                 728     <- the function itself
_platform_memcmp           497     <- string comparison
__tree...__emplace_unique  287     <- map lookup
(for comparison: malloc/free-related frames ~50 in total, get_max 25)
```

The original code **walked the tree four times per edge**: `find(obj1)` + `find(obj2)` +
`op[](obj1)` + `op[](obj2)`, each O(log n) string comparisons. Switching to
`unordered_map` and reusing the iterator `insert` returns leaves two hash lookups per edge (with no
string comparison on a hash hit):

```cpp
pair<obj_id_map::iterator, bool> o1 = map_obj_grp.insert({sims[i].obj1, map_obj_grp.size()});
...
add_grp(o1.first->second, o2.first->second, sims[i].sim);   // no op[] re-lookup
```

Ids are still assigned in first-appearance order, so clustering and output are completely unchanged
(confirmed by the differential test and a byte-for-byte comparison).

Measured init: 400k edges 0.127 -> 0.024 s and 800k edges 0.283 -> 0.053 s, both **5.3x**.

### `get_unit()`'s linear search (performance) - the single biggest win

Attribution with `-fno-inline` showed it was **about 51% of the whole program on its own** (965 of
1890 samples):

```
get_unit(...)          965
<deduplicated_symbol>  281
merge_grp(...)         101
all malloc/free frames 177
add_grp(...)            62
heap_sift_down(...)     61
```

Every index needing an update called `get_unit()` twice, which walks a chain linearly (mean chain
length 8.4; 36.7 million steps at 800k edges). But those two edges **are already on the four chains
the collection phase walks anyway**:

For index x and group g (g is either group1 or group2):

- when `x > g` the cell `(g, x)` is on the chain of **row g**;
- when `x < g` the cell `(x, g)` is on the chain of **column g**.

The collection phase walks exactly row group1, row group2, column group1 and column group2. So the
values are recorded into `sim_g1[]` / `sim_g2[]` while collecting, and the update phase reads them
directly, leaving `get_unit()` with one call per merge (to fetch the merged cell itself).

Two implementation details:

- `-1` is the sentinel for "no such edge" - safe because similarities are non-negative (the reader
  limits them to `[0,65535]`);
- **the first write wins**: `get_unit()` walks from the head of the chain (`add_grp` inserts at the
  head, i.e. the most recently inserted cell) and takes the first match, so the first write
  corresponds to the same cell. Behaviour stays identical on duplicate-pair input.

The cost is two O(groups) `double` arrays, reset only where touched; measured peak memory
**+0.5 to +1.5%**.

Clustering speedups: n=20k/40k/80k/200k -> 1.76x / 1.62x / 1.78x / **1.88x**.

### The reader buffered a `vector<sp_sim> sims` (memory)

The reader used to store **every edge** in a global `vector<sp_sim>` for the initialiser to walk -
and nothing used it after "read everything, build everything" was done.

Measured cost (1M edges, n=100k, identical topology, only the object-name length changed):

| Object names | Peak RSS |
|---|---|
| 8 characters (SSO, no extra allocation) | 202 MB |
| 64 characters (heap, about 65 bytes per string) | 371 MB |

**The 169 MB difference is that buffer.** Switching to streaming (`readdb(path, handler)`, calling
back per parsed edge) removed it entirely and halved peak RSS:

| Edges | n | Peak RSS before | After | Change |
|---|---|---|---|---|
| 200k | 20k | 44.4 MB | 26.6 MB | **-40%** |
| 400k | 40k | 86.3 MB | 46.4 MB | **-46%** |
| 800k | 80k | 174.0 MB | 88.8 MB | **-49%** |
| 2M | 200k | 401.3 MB | 222.9 MB | **-45%** |

**But purely streaming has a trap**: interleaving structure building into the parsing loop lets the
growing structure pollute the parse phase's cache working set. Measured at 2M edges, "read + build"
went from 0.313 s to 0.416 s (+33%), 7% slower overall, while clustering was unaffected
(0.884 -> 0.887) - so the problem really was that phase's locality.

Hence **small-batch buffering**: fill 65536 edges, then call back in one go, keeping both phases
batched.

| Approach | Time (2M edges, median) | Peak RSS |
|---|---|---|
| Read everything first (original) | 1.303 s | 401.3 MB |
| Purely streaming | 1.393 s | 219.4 MB |
| **Streaming + 64k batches** | **1.301 s** | **223.0 MB** |

The best of both: 44% less memory at the same speed (costing one batch buffer, 3.6 MB). The now
unused `struct sp_sim` was deleted, and names are moved from the reader into the map, saving one
copy per edge.

### Group numbers from `size_t` to `uint32_t` (memory + speed)

Group numbers live in `[0, 2*objects)`, so `uint32` suffices, and `sizeof(sp_unit)` dropped from 56
to 48 bytes.

| Scale | Peak RSS 56B | 48B | Change | Time (alternating median) | Speedup |
|---|---|---|---|---|---|
| 800k edges | 88.8 MB | 76.6 MB | **-13.7%** | 0.464 -> 0.385 s | **1.21x** |
| 2M edges | 222.9 MB | 192.4 MB | **-13.7%** | 1.474 -> 1.370 s | **1.10x** |

(The cost: an object limit of about 2.1 billion, where the address space would run out first.)

### A probe: node size is a first-order factor

To decide whether a big pointer-to-index refactor was worth it, an **amplification probe** was run:
pad the struct from 48 to 128 bytes (2.67x, same field offsets, no other code changed).

| Version | 800k edges | 2M edges |
|---|---|---|
| 48B (current) | 0.427 s | 1.485 s |
| 128B | 0.498 s (**1.167x slower**) | 1.511 s (1.018x slower) |

Conclusion: **node size really does matter**, while malloc/free does not (the node pool measured
1.00x). The two are not contradictory: the pool does not change `sizeof(sp_unit)`, so it does not
change how many cache lines a node spans.

### Links from pointers to indices: 32-byte nodes (speed + memory)

The four links and both group numbers became `uint32` indices, with the whole pool in one
contiguous block:

```
sizeof(sp_unit):  4*8(pointer) + 2*4(groups) + 8(sim) = 48   ->   4*4(links) + 2*4 + 8 = 32 bytes
```

| Scale | Time (alternating median) | Speedup | Peak RSS | Change |
|---|---|---|---|---|
| 800k edges | 0.393 -> 0.266 s | **1.49x** | 76.6 -> 63.3 MB | **-17.3%** |
| 2M edges | 1.337 -> 0.954 s | **1.40x** | 189.0 -> 155.2 MB | **-17.9%** |

(Min, lower quartile and median all agreed, so this is not noise.)

This change had **already been rejected once**. The "pad to 128 bytes" probe had suggested size was
a second-order factor, but those two measurements used only 7 and 5 repetitions, right in the
noise - see methodology lesson 2 below. Re-measured with 11 runs, size turned out to be first-order,
which turned the refactor from "not worth it" into "worth it".

**Two implementation traps, both recorded:**

1. **An empty link must be 0, and node 0 must be reserved as a sentinel.**
   `vec_dim1`/`vec_dim2` are **value-initialised (zeroed)** when they grow. When `head`/`max` were
   raw pointers, zero was exactly NULL; with indices, 0 became a **valid node index**, so every
   empty row pointed at node 0 and the chains fell apart (as an infinite loop). Now `SP_NONE = 0`
   and node 0 is never handed out, which is consistent with value-initialisation and needs fewer
   branches.

2. **Manage the pool with `malloc`/`realloc`, not `std::vector`.**
   Growing a vector means "new block + copy + release old", and the allocator does not return large
   blocks to the OS promptly, so every move left a resident intermediate block: growth by 1.25x cost
   **150 MB** (for 26 MB of nodes in total!) and even 2x cost 85 MB. `realloc` can usually extend
   large blocks in place with neither a copy nor an intermediate block, taking the peak down to
   63 MB.

**ASan's use-after-free detection was not lost.** The pool's memory stays accessible, so ASan could
not see "a node was freed and then read" - the very bug class that had appeared twice in this code
base. So a node is poisoned with `__asan_poison_memory_region` on release and unpoisoned on
allocation; when the pool moves the whole block is unpoisoned first and re-poisoned from the free
list afterwards. **Free indices live in their own array rather than borrowing fields from the
nodes** - otherwise walking the free list would itself hit poisoned memory.

This defence was verified: deliberately reading a released node inside `get_max()` makes ASan report
`use-after-poison ... in get_max()` immediately (see also case 22 in `tests/`, a use-after-free that
was really hit).

### Single linkage moved to Kruskal (speed + memory) - 300x

Single linkage's inter-cluster similarity is the maximum over all **original edges** between two
clusters (the max is not diluted), so the whole merge tree equals **sorting the edges by similarity
descending and merging with a union-find** (Kruskal / MST). `LINKAGE_MTD == 1` therefore needs no
similarity matrix at all: m edges plus one sort.

| n | Edges | Greedy (old) | Kruskal (new) | Speedup | Greedy peak RSS | Kruskal |
|---|---|---|---|---|---|---|
| 20,000 | 200k | 1.767 s | 0.0277 s | **63.7x** | 2261 MB | **16.5 MB** |
| 40,000 | 400k | 7.428 s | 0.0581 s | **127.9x** | 6251 MB | **28.1 MB** |
| 80,000 | 800k | 35.416 s | 0.1180 s | **300.2x** | **15765 MB** | **50.3 MB** |

The speedup grows with scale (the old path was O(n²), the new one O(m log m)) and memory differs by
two orders of magnitude - that 15.8 GB was the similarity graph being inflated into a near-complete
graph.

**Merge order under ties may differ from the old implementation**, but under single linkage the
clusters at any threshold equal "the connected components after keeping only edges with
`sim >= X`", which is independent of processing order. So the partition at every threshold, and
every cluster's cohesion (the largest intra-cluster edge), are completely unchanged. That was tested
as an acceptance criterion: **23 inputs x 10 thresholds = 230 combinations, each identical to the
old greedy implementation**. The differential suite still passes as usual (`max` 31 PASS results
plus the corresponding part of the 900 cut checks).

### `grp[].objs` was only used for `size()`: replaced by a member count (memory -28%)

A live group used to carry a `vector<string>` member list. But member names now come from the
**merge tree** (`merges`) plus `object_names` (the CLI looks names up by id after cutting), so that
list was **only used to get `size()`**; `tools/spclust_cli.cpp` does not even read `grp`. It was
replaced by an `n_members` count.

What that removed was not only the copies but also a piece of overlooked resident memory:
`grp[g].objs.clear()` clears elements while **keeping capacity**, and `grp` grows to 2n slots - those
retired slots held the member lists' capacity forever.

| | Median wall clock (n=200k / 1M edges) | Peak RSS |
|---|---|---|
| Before | 3.511 s | 183.6 MB |
| After | 3.285 s | **132.0 MB** |
| Change | 1.069x | **-28.1%** |

Two small cleanups came along:

- After reading, the name-to-id map (`unordered_map`) is released - it is never needed again, and in
  a long-running library scenario this returns about 7 MB of resident memory. (**It does not lower
  the peak**: the peak occurs during reading, when the map and the structure are both present.)
- A **sampling attribution that did not hold** was confirmed: the profile showed
  `vector<unsigned int>::__base_destruct_at_end` at 7%, which looked like free-list churn, but
  replacing it with a raw array plus a top index measured **zero gain**. Sampling points a
  direction; only A/B measurement establishes a gain.

### A batch of small optimisations: only two of five survived

The intent was to squeeze the remaining odds and ends (expecting 5-10% in total). **After A/B-ing
each one, only two held**:

| Change | Expected | Measured | Verdict |
|---|---|---|---|
| pre-size `mark_flags`/`sim_g1`/`sim_g2` to the whole id space | 1-3% | +1.1-2.1% (bench_mid) | **kept** |
| `add_grp`: seven field writes -> one aggregate assignment | 1-3% | (measured with the above) | **kept** |
| `sort(updated)`/`sort(dirty_max)` -> insertion sort | 1-2% | **-21%** | **reverted** |
| heap compaction threshold 4x live -> 2x live | 1-3% | 0% (within noise) | **reverted** |
| hand-rolled scanning in the reader instead of `split_fields` | 3-5% | ceiling only ~2%, not done | **not done** |

**The insertion-sort mistake is worth recording**: the mean `get_unit` chain length of 8.4 was
mistaken for the number of indices affected per merge. The latter measures **38** (1,402,219
iterations over 36,564 merges), and insertion sort is O(k²) - at k ~ 38 that is 7x slower than
introsort, i.e. 21% slower overall. That number is now in a code comment so nobody "optimises" it
again.

**Why the reader was left alone**: read plus structure building is 15-22% of the total (bench_large:
0.114 s of 0.749 s), which looks significant, but `add_grp` is about 7% of it and name hashing
about 4%, leaving roughly 4% for parsing itself (the string construction in `split_fields`/`trim`).
Halving parsing would therefore gain about 2% - right at the edge of what repeated measurement can
resolve on this machine, and not worth rewriting the parse path and risking behaviour drift.

### A final round of library hygiene and tidiness

| Item | Change |
|---|---|
| `readdb` error handling | from an internal `exit(1)` to **returning 0 on failure** (reason on stderr), leaving the exit code to the caller - otherwise it cannot be used as a library |
| cohesion under min linkage | a singleton's initial `sim` became `EMPTY_COHESION`: the min recurrence is `min(...)` and 0 is its absorbing element, so the reported cohesion used to be identically 0. Now `LINKAGE=2` reports the smallest intra-cluster edge (sim.csv: `85 / 82 / 76`, previously `0 / 0 / 0`), with the **partition completely unchanged** |
| unused macros removed | `MIN_ACCT_NUM` / `MIN_SIM_WG` |
| `--tree-depth <N>` | the dendrogram depth limit became a parameter (default 6, 0 = unlimited, implies `--tree`) |
| headers no longer pollute namespaces | `src/hac.cpp` dropped `using namespace std;` and qualifies types as `std::`; `include/spclust/hac.h` / `readdb.h` were qualified too; each .cpp does its own `using` |
| no more wiring globals across files | `grp` moved from "defined in main.cpp, `extern` in spclust.cpp" into `src/hac.cpp` as a file static |

The exit code contract is now pinned by the CLI test as well:

| Situation | Exit code |
|---|---|
| `--help` | 0 |
| unknown option / extra argument / `--cut`, `--tree-depth` without a value | 2 |
| input file cannot be opened / IO error / no records at all | 1 |

**Two infrastructure problems were fixed as a side effect, both self-inflicted in that round:**

1. `make test` / `make bench` depended on `linkages` but **not on `./spclust`**, so the CLI test
   could run against a stale binary. It happened: a patch adding a `readdb` return-value check to
   `tools/spclust_cli.cpp` did not apply, and the CLI test **still passed** (it was testing the
   previous build). Now `test`/`bench` both depend on `$(BIN)`.
2. That patch failed because its `oldText` no longer matched (the call site had changed from
   `readdb(argc>1?...)` to `readdb(path, ...)`) and **the assertion was missing**, so it silently
   did nothing. Lesson: patch scripts always get an assertion.

## Known issues

Everything below comes from measurement or static analysis. The first two **have been reproduced by
`make test`** (with the corresponding marks); the rest are static analysis only.

### 1. Duplicate pairs: semantics defined (the maximum wins), so the two sides no longer disagree

When the same `(obj1, obj2)` appears more than once with different `sim`, the semantics are **the
larger value wins**. This is not a guess: taking the maximum and swapping the two lines both merge
the same pair with value 90; internally the merge value comes from `get_max()` (the row maximum) and
derived similarities are recorded by maximum too (`record_sim_g1/g2`). The oracle models the same
semantics, so **this divergence no longer exists** (the six known divergences that were registered
became XPASS and were deleted).

Both cells still go into the structure (`get_unit` may return either, but only for the "release
later" mark, so the result is unaffected); merging neither crashes nor leaks. Cases:
`tests/cases/07_duplicate_edge.csv`, `19_duplicate_order.csv`.

Changing this also unified the implementation: derived similarities used to take "the first on the
chain" (i.e. the most recently inserted), so "the merge value of this pair" took the maximum while
"the new similarities it derives" took the last inserted - two semantics for one duplicate input.
Everything takes the maximum now.

### 2. Self-loops (`obj1 == obj2`) break the structural premise, and the reader drops them

All of `add_grp`'s logic rests on `group1 < group2`, and a self-loop makes `group1 == group2`, so the
row and column chains stomp on each other. `readdb.cpp` prints a warning and skips such lines, which
is **necessary**, not optional.

### 3. `max`/`min` (now `Linkage::MaxSimilarity` / `MinNonzero`) are not standard linkages - both handled

`1` (commented "complete") takes `max`, which in standard distance terms is **single linkage**; `2`
(commented "single") takes `min` and requires both sides to exist, a variant of **complete linkage**.
Note that this naming matches the sparseAHC documentation exactly - in the "similarity matrix"
branch it is conventional and not necessarily a typo.

Neither is a problem any more:

- `1` (single linkage) used to densify on sparse graphs: the max is not diluted, so `sim(i,G)` is
  positive as long as cluster i has **any** edge to G, turning the cluster-level similarity graph
  into the **transitive closure** of the input (a near-complete graph) - measured at n=40000 /
  400k edges it built 93,418,345 cells (234x the edge count) in 7.35 s, and at n=80000 / 800k edges
  peak memory was **15.8 GB**. It now goes through **Kruskal**: single linkage's inter-cluster
  similarity is the maximum over all **original edges**, so the whole tree equals "edges sorted by
  similarity descending plus union-find" and no matrix is built. See "Fixed" above.

- `2`'s `grp[].sim` was identically 0 (see item 4), affecting only the cohesion report.

Separately, `add_grp(i, grp_size, max(sim1, sim2))` originally had no `if (sim > 0)` guard. Adding
it is **correct** (a value of 0 means there is no edge at all between the clusters, and `get_max()`
never picks a cell with similarity <= 0, so no result changes) but **it was not the cause of the
performance problem** - measured, that guard removes only 52,147 / 93,418,345 = 0.06% of the cells.
Attributing the 8x slowdown to it was wrong.

### 4. Under `min` (now `Linkage::MinNonzero`), `grp[].sim` is identically 0

A singleton's `sim` starts at 0, and `min(min(...), max_sim)` is 0 forever once it meets a 0. This
affects only the cohesion report, not the merge order.

### 5. "An absent edge is 0" inherently dilutes inter-cluster similarity

Not a bug, but worth documenting: average linkage's denominator is all `|A|·|B|` pairs, so the larger
the clusters the more diluted they get. Case `05_chain_bridge.csv` is an instance: `a3-b1 = 58` is
the only above-threshold edge between the two groups, but `a1-a3` is absent, so after one merge the
similarity of `{a1,a2}` to `a3` is only `(0+94)/2 = 47` and the bridge sinks below the threshold
before it can trigger a merge. sparseAHC uses the same denominator, so this is not an implementation
difference but inherent to the model; if your cluster shape should be decided by connectivity rather
than by a mean, you have to choose a different normalisation (for example dividing by the number of
edges that actually exist).

### 6. Static analysis only; the tests do not distinguish these (results are unaffected)

- `grp` / `vec_dim1` / `vec_dim2` grow monotonically to about 2n-1; merged slots are cleared but
  never reclaimed.
- `MIN_ACCT_NUM` and `MIN_SIM_WG` were defined but never used. -> deleted.
- `using namespace std;` lived in headers and `grp` was wired up with globals plus `extern`. ->
  fixed: the headers no longer say `using namespace std`, and all 34 file statics moved into
  `spclust::detail::HACImpl` (one object, one state). See
  [design](design.md#state-and-instances).

### Verification record

- `make test` now runs two parts: the differential test (algorithm) and the CLI test (interface
  contract, 11 checks all passing, including the exit-code contract).
- `make bench`: the performance guard, matching the baseline on two consecutive runs (0.90-1.09x);
  a synthetic `-O0` regression is caught.
- `make asan` / `make asan-cases`: no memory errors and no UB on `sim.csv` and all 22 cases.
- `make leaks`: `0 leaks for 0 total leaked bytes` (`free_mem()` releases every list node).
- `make test`: `cut` = 20 / 50 / 70 / 200 all give `RESULT: PASS`, with `TRUNC=0` and `KNOWN=0`.
- Heap correctness: `LINKAGE=2` on `rand_n44d25s1.csv` once caught "the heap top is not the global
  maximum"; after the fix that case is byte-identical along with every other input.
- The differential suite now has **no KNOWN divergences at all**: the one that was registered
  (duplicate-pair semantics) disappeared once the semantics were fixed to "the maximum wins" and the
  oracle was aligned, and the six entries were deleted.
- Equivalence: 22 cases + 7 random graphs + `sim.csv`, 31 inputs x 3 linkages, full DEBUG output
  (including the similarity of every merge) is **byte-identical**; only the printed whole seconds
  differ.
- Performance: alternating runs with the minimum wall clock taken (>= 11 runs, min/lower
  quartile/median agreeing). 800k edges 5.906 -> 0.270 s (**21.9x**); 2M edges 47.204 -> 0.925 s
  (**51.1x**).
- Memory: replacing `grp[].objs` with a member count cut the peak by a further **28%**
  (183.6 -> 132.0 MB at 1M edges); streaming the reader cut peak RSS by **40-49%** (2M edges:
  401 -> 223 MB) at no time cost; shrinking nodes from 56 to 48 bytes cut a further **13.7%**
  (223 -> 192 MB) while running 1.10-1.21x faster.

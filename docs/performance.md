# Measurements and methodology

Every figure here was measured, and the occasions where the measurement was **wrong** are
recorded too - machine noise is enough to fool a single run. That is why the protocol
(alternating runs, >= 11 repetitions, min/lower-quartile/median agreeing) and the baseline are
pinned in the repository (`make bench`, see [testing](testing.md)).

## Measurements

Compiler options: `-O2`.

**A caveat for this whole section**: absolute wall clock depends heavily on machine load - the
same code can differ by 2x between a busy/hot run and an idle one (the canary mechanism is
described in [testing](testing.md)). **The ratios come from alternating runs and are
trustworthy**; absolute numbers are only comparable on the same machine under the same load.
The canonical baseline is `tests/bench_baseline.txt`, recorded with the machine idle.

### Memory: linear, consistent with sparsity, and far below where it started

| Edges | n | Peak RSS | Bytes/edge | vs baseline |
|---|---|---|---|---|
| 200k | 20k | 20.7 MB | 104 | -50% |
| 400k | 40k | 34.3 MB | 86 | -60% |
| 800k | 80k | 64.2 MB | 80 | -61% |
| 2M | 200k | 155.2 MB | 78 | -61% |

Bytes per edge converges to a constant, which is what a cross-linked list storing only the
non-zeros should do. Three memory changes stack up: the reader no longer buffers every edge
(-45 to -49%), nodes went 56 -> 48 bytes (-13.7%), and 48 -> 32 bytes with the pool switched to
`realloc` (-17.9%).

### Time (by phase; alternating repetitions, minimum taken, >= 11 runs)

Input: a random sparse graph, n = 80000, 800,000 edges, mean degree 20, similarities 0..100,
`MIN_SIM_CLUST=50`.

| Version | Read | Init | Cluster | Total | This step | Cumulative |
|---|---|---|---|---|---|---|
| baseline (after the UAF fix) | 0.050 | 0.306 | 5.574 | 5.930 | - | 1.00x |
| P0-1 incremental index walk | 0.049 | 0.316 | 2.475 | 2.840 | 2.07x | 2.07x |
| P0-2 max-heap | 0.049 | 0.311 | 0.531 | 0.891 | 3.09x | 6.41x |
| init via `unordered_map` | 0.050 | 0.060 | 0.512 | 0.622 | 1.41x | 9.01x |
| no `get_unit` linear search | 0.049 | 0.058 | 0.295 | 0.402 | 1.49x | 13.46x |
| streaming reader + batching | ~0.05 (merged) | - | ~0.29 | ~0.40 | ~1.00x | ~14x |
| group ids as `uint32` (48B) | ~0.06 | | ~0.26 | ~0.36 | ~1.15x | ~16x |
| links as indices (32B pool) | ~0.05 | | ~0.27 | ~0.33 | ~1.10x | ~18x |

(The first three rows are per-version measurements; in the last three the streaming reader
merged "read" into "init", so only totals are estimated.)

Cumulative, alternating runs with the minimum taken:

| n | Edges | Time, baseline | Time, now | Speedup | Memory, baseline | Memory, now | Memory |
|---|---|---|---|---|---|---|---|
| 20,000 | 200k | 0.435 s | 0.050 s | 8.7x | 41.3 MB | 20.7 MB | **-50%** |
| 40,000 | 400k | 1.565 s | 0.111 s | 14.2x | 85.5 MB | 34.3 MB | **-60%** |
| 80,000 | 800k | 5.906 s | 0.270 s | **21.9x** | 163.9 MB | 64.2 MB | **-61%** |
| 200,000 | 2M | 47.204 s | 0.925 s | **51.1x** | 396.4 MB | 155.2 MB | **-61%** |

### Where the time goes now

Sampled profile (`sample`, 1 ms, attributed with `-fno-inline`, n = 200k / 2M edges):

| Frame | Samples | |
|---|---|---|
| `<deduplicated_symbol>` | 387 | unsymbolised runtime and merged code blocks |
| `merge_grp()` itself | 189 | row/column removal, the update loop skeleton |
| `add_grp()` | 131 | cell initialisation, two chain insertions, the heap push |
| `__hash_table::emplace_unique` | 101 | name-to-id lookup in init (already optimised once) |
| `set_flag()` | 90 | marking plus recording the cached values |
| `std::basic_string::__is_long` | 89 | string copies (building the name table, member lists) |
| malloc/free-related frames, total | 223 | about 13% (but the node pool measured no gain - see below) |
| `heap_sift_down()` | 57 | |
| `__write_nocancel` | 51 | writing the final few statistics even with `DEBUG=0` |

It is flat now: nothing exceeds 25% and there is no single hot spot.

### The operation counts are fully linear; what remains super-linear is cache behaviour

Update-loop iterations scale strictly linearly with n (the 36.7 million `get_unit` steps were
equally linear and have since been eliminated):

| Edges | Update iterations | Factor | (the old `get_unit` chain walk) |
|---|---|---|---|
| 100k | 174,906 | - | 4,576,072 |
| 200k | 351,062 | 2.0x | 9,175,748 |
| 400k | 701,373 | 2.0x | 18,305,037 |
| 800k | 1,402,219 | 2.0x | 36,673,256 |

Yet throughput per unit time falls from 176 Mops/s (n = 10k) to 58 Mops/s (n = 80k), a 3x
decay. So the ~n^1.5 seen in wall clock is **pointer chasing cache-missing on a growing working
set**, not some remaining O(n) algorithmic term. The node-pool experiment shows it is not
allocator overhead either (2.19M malloc/free calls eliminated, no change in time whatsoever).

### Other hypotheses, all ruled out by measurement rather than by argument

| Suspect | Measured result |
|---|---|
| `malloc`/`free` overhead | a full node pool gave **0.998x** (800k edges, 11 alternating runs) -> about 2%, **not the problem** |
| string copies in `grp[].objs` | 8-char vs 64-char object names: 2.915 vs 2.915 s -> **not the problem** |
| memset of the mark vector | `size_t` -> `char`: 2.91 -> 2.80 s (4%) -> not the main cause |
| `get_max` | minor before P0-1, then 57-78% of clustering, and only 6.5% after P0-2 |
| object-name to id lookup | **a main cause** (init 5.3x faster) -> optimised |
| `get_unit` linear search | **a main cause** (51% of the whole program) -> optimised (clustering 1.88x) |
| `LINKAGE_MTD=1` | 2.9 s -> **24.0 s** (8.3x) -> unconditional `add_grp` densification is real |

## Optimisations that were tried and gave nothing

### Replacing malloc/free with a node arena - measured 1.00x, reverted

The argument looked strong: with 800k edges `add_grp` is called **2,194,304 times**, i.e. 2.19
million malloc/free pairs, each 56-byte cell carrying allocator metadata. So a "4096 cells per
block plus a free list" pool was implemented, with a `-DSPCLUST_USE_SYSTEM_MALLOC` switch so ASan
could still detect use-after-free.

**It bought nothing.** Alternating runs, minimum taken:

| Scale | Without the pool | With the pool | Speedup |
|---|---|---|---|
| 800k edges (11 runs each) | 0.5626 s | 0.5638 s | 0.998x |
| 1.6M edges (7 runs each) | 1.3924 s | 1.3895 s | 1.002x |

Peak memory did not move either (0.0%) - the peak is set by the same set of live cells
regardless of how they are allocated.

There are two layers to why, and both are now clear:

- malloc/free-related frames are only about **2%** of the sampled profile; the small-object
  allocator is simply fast;
- more fundamentally, **the pool does not change `sizeof(sp_unit)`**, so it does not change how
  many cache lines a node spans. And node size is a first-order factor (padding a node to 128
  bytes costs 1.17x; squeezing it to 32 bytes gains 1.4-1.5x). So "change the allocator" does
  nothing while "shrink the node" does - and those two results do not contradict each other.

The cost, on the other hand, was real: about 60 extra lines, a second build configuration, and
**in a default build ASan loses the ability to detect "a cell was freed and then read"** - the
very bug class that had already appeared twice in this code base. Zero benefit, non-zero cost:
reverted.

**Methodology lessons (two mistakes were made on this one, both recorded):**

1. **Alternate first, then repeat.** The first measurement of 1.31-1.60x was fake - it ran all
   of A and then all of B, and load drift fooled it; two builds of the same code differed by 22%
   in single measurements.
2. **Too few repetitions produce confidently wrong conclusions.** On the 56 -> 48 byte change,
   7 runs said +7.7% and 5 runs said -2%, while 15 runs settled at **+20.9%** (min, lower
   quartile and median all agreeing). With background activity on the machine, too few runs never
   reach the true minimum, so the conclusion flips around in the noise.
   **The rule now: >= 11 alternating runs, and min/lower-quartile/median must agree before any
   conclusion is drawn.**

   Which means the node-pool conclusion above, although it reproduced 1.00x three times and has a
   mechanism behind it (the pool does not change `sizeof`), strictly speaking should be re-run
   with >= 11 repetitions before being treated as settled.

3. **Profiles come from `sample`, not from instrumentation.** Instrumenting `get_unit`/`add_grp`
   with `clock()` - functions called millions of times - added thirteen million `clock()` calls;
   total time went from 0.72 s to 3.00 s and the share data was meaningless (use the system
   sampler instead: `sample`, `perf record`, or similar).

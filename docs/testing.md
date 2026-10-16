# Testing

Six layers, each watching something different:

| Layer | Where | Watches | Size |
|---|---|---|---|
| cut property test | `tests/test_cut.cpp` | the library's three cut entry points + multiple instances | `PASS=18281` |
| differential test | `tests/run_tests.py` + `tests/oracle.py` | whether the merge sequence is right (the algorithm) | `PASS=90 WARN=12 FAIL=0` |
| CLI contract test | `tests/test_cli.py` | exit codes, options, output invariants | `11/11` |
| scipy interop | `tests/test_scipy.py` | format and semantics really line up (real scipy as the ruler) | 46 combinations |
| C API | `tests/test_c_api.c` + `tests/test_c_api.py` | the header is plain C, the ABI and error contract, FFI can dlopen it | 64 + 31 checks |
| performance guard | `tests/bench.py` | that it has not got slower (anti-regression) | 5 cases |

The first five are inside `make test` / `ctest`; the sixth is `make bench`.

## Cut property test (library level)

`tests/test_cut.cpp` looks only at the library, never at the CLI. The point is that the
**expectations must not be a second copy of the union-find inside the library** - if both sides
made the same mistake, the test would see nothing. A different route is used instead: every node
of a merge tree is naturally a leaf set, so "cut after the first i merges" is "the leaf sets of the
nodes with no parent in that prefix". All three cut modes can be computed that way directly and
then compared, one by one, with what the library returns.

Coverage: three linkages x 21 random graphs (2-40 objects, dense/sparse/with duplicate edges),
plus tied thresholds, `k > n`, `cut(0)`, monotonicity (every cluster of `k-1` must be an exact
union of clusters from `k`), and the rule that a cluster's cohesion equals the value recorded by
the merge that created it. Finally there is the **multiple instances** section: two datasets
across three objects, interleaved edges and interleaved `build()` calls, must reproduce the solo
runs byte for byte - the property that mattered most after moving 34 file statics into an object.

One structural piece of evidence comes along for free: `.data`/`.bss` of `libspclust.a`
(according to `nm`) contains no algorithm state, only constants such as string literals and jump
tables, so `HAC` objects cannot secretly share anything.

## C API (real C compiler + real FFI)

`tests/test_c_api.c` is **real C** (compiled by `$(CC)`), so it incidentally proves the header
really is plain C and that the `extern "C"` symbols work from C. `tests/test_c_api.py` then calls
a shared build of the same library from Python via ctypes, verifying dynamic loading, symbol
names, struct layout and - most importantly - that **C++ exceptions do not cross the C boundary**:
a self-loop passed from Python must yield `SPCLUST_ERR_INVALID_ARGUMENT` rather than aborting the
process.

## scipy interop (a foreign ruler)

The layers above all use an oracle we wrote ourselves. `tests/test_scipy.py` swaps in scipy - a
completely differently sourced, widely used implementation - and checks two things: **format
validity** (`is_valid_linkage` / `is_monotonic` / `to_tree` / `fcluster` all usable directly) and
**semantic agreement** (on cases with no ties and a complete tree, row for row identical to
`scipy.cluster.hierarchy.linkage`, heights and sizes included).

The key prerequisite is choosing the distance mapping correctly: this library's "an absent edge
is similarity 0" corresponds to "distance = C - similarity" (C = the largest similarity), which
makes `average` linkage equivalent to scipy's `method='average'` on the same dense distance
matrix.

Two **known differences** are tested here as well, so they cannot quietly turn into bugs:
`MinNonzero` must equal scipy's `complete` on dense input (6/6) and must differ on sparse input;
and a forest (an incomplete tree) must warn on stderr, with `--scipy-bridge` producing exactly
n-1 rows and a last size of exactly n.

Without numpy/scipy the whole test skips and returns 0, so `make test` works without them.

## Differential test

The core (`src/hac.cpp`) is hand-written `malloc`/`free` with four-way pointers and a cross-linked
list, and correctness of that kind of code cannot be established by reading it. So `make test`
performs a **three-way cross-check**:

```
1. oracle(brute)  tests/oracle.py  recomputes from the raw dense matrix by definition at every
                                   step. Most trustworthy, O(n^4).
2. oracle(fast)   tests/oracle.py  maintains the aggregates incrementally. O(n^3).
3. spclust                         the subject under test.
```

`oracle.py` contains **no** Lance-Williams recurrence and no linked lists - only "compute by
definition and take the global argmax" - which makes it independent. All three implementations
use exact integer arithmetic (`Fraction`), so no floating-point error is introduced.

The rule: spclust and the oracle must agree on the **merge sequence** (member sets plus
similarity) and on the final partition. Disagreements fall into six classes:

| Mark | Meaning |
|---|---|
| `PASS` | exactly consistent with the definition |
| `WARN` | a tie for the maximum (or a tie created by truncation) makes the merge order unspecified; both sides may be legal |
| `TRUNC` | the deviation is fully explained by `unsigned short` integer truncation (see known issue 1) |
| `KNOWN` | a divergence listed in `tests/known_divergences.txt` and accepted as explicable |
| `DRIFT` / `FAIL` / `ORACLE` | unexplained deviation / different results / the harness itself is inconsistent - all three fail `make test` |

The oracle's own credibility rests on brute and fast cross-checking each other (on small inputs
they must agree merge by merge; the number of such checks is printed in the report).

Cases are `tests/cases/*.csv` (22 hand-written boundary cases, each with a comment saying what it
probes) plus 12 random sparse graphs (with distinct similarities, to avoid ties). The whole set
runs in about a second.

```bash
make test                      # default min_sim=50
make test MIN_SIM=70           # rerun at another threshold (the oracle follows, and the binary is rebuilt)
python3 tests/run_tests.py --only average --verbose
python3 tests/run_tests.py --quick --case 15_   # only matching cases
python3 tests/oracle.py sim.csv --linkage average --cut 50    # run the oracle alone (--cut cuts the tree)
```

Current results (`cut=50`, 22 hand-written + 12 random x 3 linkages = 102 comparisons):

```
PASS=90  WARN=12  KNOWN=0  TRUNC=0  DRIFT=0  XPASS=0  FAIL=0  ORACLE=0
oracle self-checks (brute/fast agree + cut == early stop): 198
multi-threshold cut consistency checks: 900
RESULT: PASS
```

What is compared is the **full merge tree** (the runner passes `--floor 0` to spclust and
`min_sim=0` to the oracle), which is much stronger than comparing only the part above the cut
threshold. On top of that, ten thresholds are cut separately for cross-validation (900 = 90 PASS
results x 10 thresholds). **KNOWN is now empty** - the one divergence that was ever registered
(duplicate-pair semantics) disappeared once the semantics were settled and both sides were
aligned. WARN rose from 8 to 12 because a tree reaching down to low similarities has more ties,
which is normal.

The conclusion: **the linked-list maintenance in `src/hac.cpp` is reliable** - across 44 objects,
43 consecutive merges and seven input shapes (dense complete graph, star, long chain, several
connected components, zero-similarity edges, reversed input, fractional similarities) the 64
comparisons under max/min linkage show only WARNs (ties) and the registered duplicate-pair
divergence.

After the fix, average linkage's residual against the exact definition is purely floating-point:

```
excluding the two duplicate-pair semantic cases, 146 per-merge similarity comparisons:
  largest absolute error = 3.3e-11
  largest relative error = 4.3e-13   (the DEBUG log keeps 12 significant digits, so this is
                                      essentially the print precision)
```

`TRUNC` now serves as a **regression guard**: if anyone put similarities back to `unsigned
short`, those 12 `TRUNC` results would reappear, and `TRUNC` counts as **failure**.

### CLI interface test (`tests/test_cli.py`)

The differential test above verifies the **algorithm** (parsing the merge log and comparing with
the oracle), while the command-line interface - `--cut/--floor/--tree/--matrix`, exit codes, the
stdin mode - had no coverage at all. The two also tend to mask each other: the program and the
runner were implicitly coupled through the merge log's format, so one changed output line could
make the runner silently under-count (it happened here: the Kruskal path forgot to print that line
and the runner reported "0 merges").

The second part of `make test` tests only the interface contract, not the algorithm again:

| Check | Contents |
|---|---|
| exit codes and error paths | `--help` -> 0; unknown option / extra argument / `--cut` without a value -> 2; missing file / empty input -> 1 |
| `--cut` matches the oracle | 22 cases x 5 thresholds; the CLI partition must equal the oracle's full tree cut at the same place |
| cut == early stop | `--floor 0 --cut X` must equal `--floor X --cut X` (monotonicity at the CLI level) |
| `--floor` / `--cut` interaction | with only `--cut 10` given, the floor must drop to 10 automatically, or the cut is unreachable |
| `--matrix` structural invariants | each row's child ids must be smaller and used as a child only once; size >= 2; similarity non-increasing; root sizes summing to the object count |
| `--tree` completeness | every object appears exactly once; the root count matches `--matrix` |
| stdin equals file | file / `-` / no argument must all give the same result |
| duplicate pairs are order-independent | swapping two input lines changes nothing (maximum-wins semantics) |

Writing this test **immediately found a real interface problem**: `--nope` was treated as a
filename, so it was reported as "cannot open file" (exit code 1), conflating a usage error with an
input error. Unknown options now return 2.

### Performance guard (`tests/bench.py`, `make bench`)

This project has gained about 39x in time and -55% in memory in total, but **nothing was
preventing it from regressing**, and machine noise is enough to fool a single measurement (noise
fooled it twice during development and sampling attribution once). So the protocol and the baseline
are pinned:

- **Fixed inputs**: the parameters in `tests/gen_bench.py` are frozen. This is necessary - runtime
  depends strongly on the similarity distribution (with the same n and edge count, a range of
  0..100 and 1..65535 measured over 3x apart), so changing the seed or the range makes before and
  after incomparable.
- **The canary only catches "obvious" state changes.** It was tested: a baseline recorded on a
  busy/hot machine, rerun idle, showed everything "2x faster" (the canary read 0.57x and correctly
  refused to conclude). But another drift was only about 20%, which the canary (0.89x) let
  through, so the long cases showed "18-34% faster" - while an alternating A/B against the
  pre-refactor binary showed the real difference was 1.00-1.01x. **So `make bench` is a smoke
  alarm, not a measuring instrument**: do not treat 1.2x as a regression, re-record the baseline
  with `make bench-update` whenever the machine state changes, and to compare two versions
  precisely, alternate them as described above.
- **Protocol**: each case is repeated N times and the median taken; **if min / lower quartile /
  median differ by more than 15% it is marked noisy** and no conclusion is drawn. Tolerance is
  25% (time) / 15% (memory); only exceeding **3x time / 2x memory** fails (`--strict` fails on
  tolerance instead).
- **Baseline**: `tests/bench_baseline.txt`, re-recorded with `make bench-update` (after a machine
  change or an intentional performance change).

**The machine-state canary.** During development a baseline was recorded while the machine was
busy/hot (that round had just run a 5 GB scipy comparison); rerunning it idle showed everything
"2x faster", which looked like a huge optimisation but was the environment changing. The fastest
case, `bench_small`, is now used as a probe: if it deviates from the baseline by more than +-30%,
**no comparison with the baseline is concluded at all**, and the user is told to rerun or re-record.

The baseline is recorded **with the machine idle** (the file records which machine produced it, and
`make bench` warns when the numbers are not directly comparable). The current baseline:

| Case | Configuration | Median wall clock | Peak RSS |
|---|---|---|---|
| bench_small | n=2000, 20k edges | 0.009 s | 4.7 MB |
| bench_mid | n=20000, 200k edges, sim 1..65535 | 0.106 s | 20.7 MB |
| bench_large | n=50000, 500k edges, sim 1..65535 | 0.402 s | 47.7 MB |
| bench_lowval | n=20000, 200k edges, sim 0..100 | 0.049 s | 20.5 MB |
| bench_single | n=50000, 500k edges (single linkage / Kruskal) | 0.129 s | 30.2 MB |

Two consecutive normal runs sit within 0.95-1.03x of the baseline, with the canary at
0.98-1.00x.

The guard itself was verified too: rebuilding at `-O0` (a synthetic regression) had three cases
reported as `REGRESSION` (5.1-6.1x) and exiting with code 1, while consecutive idle runs stayed
within 0.90-1.09x in time and 0.98-1.03x in memory.

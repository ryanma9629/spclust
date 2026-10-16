# spclust - build
#
#   make            build ./spclust
#   make run        build and run sim.csv
#   make test       differential + cut properties + CLI + scipy interop + C API
#   make lib-shared build libspclust.so (for FFI/ctypes)
#   make check      run sim.csv once per linkage for eyeball comparison
#   make asan       build and run under AddressSanitizer + UBSan (worth doing once)
#   make asan-cases run every case under tests/cases under ASan/UBSan
#   make leaks      check for leaks with the leaks tool on macOS
#   make clean
#   make help
#
# Overridable variables:
#   make CXX=g++-15                    use another compiler
#   make run LINKAGE=1 MIN_SIM=60      LINKAGE is passed to --linkage now (run time, no recompile)
#   make test MIN_SIM=70               rerun the differential test at another threshold (the oracle follows)
#   make CXXFLAGS="-std=c++11 -O2"     take over the compile options entirely
#
# LINKAGE: 0=average  1=max-similarity (= single linkage in distance terms, Kruskal)  2=min-nonzero
#   The linkage is a **run-time** argument (spclust::Linkage), so this variable is only
#   passed to the CLI - one binary per linkage is no longer needed.
#
# Note: MIN_SIM_CLUST / DEBUG are the CLI's own compile-time defaults (wrapped in #ifndef
# in tools/cli_config.h) and -D wins; the library itself has no compile-time options.

CC       ?= cc
CXX      = c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
LDFLAGS  ?=
LDLIBS   ?=
PYTHON   ?= python3
BENCH_ARGS ?=

LINKAGE ?= 0
MIN_SIM ?= 50
DEBUG   ?= 1

DEFS = -DMIN_SIM_CLUST=$(MIN_SIM) -DDEBUG=$(DEBUG)

# Layout: the library (src/ + include/) is separate from its user (tools/).
# The build system is CMakeLists.txt; this Makefile is a development shortcut.
LIB_SRCS = src/hac.cpp src/scipy.cpp src/spclust_c.cpp
CLI_SRCS = tools/spclust_cli.cpp tools/readdb.cpp
ALL_SRCS = $(CLI_SRCS) $(LIB_SRCS)
# Every header, not just the one the library used to be: leaving a public header out
# means editing it does not trigger a rebuild, and a stale library then reports a stale
# version or misses a new entry point.
HDRS = $(wildcard include/spclust/*.h) tools/readdb.h tools/cli_config.h
INCS = -Iinclude -Itools
DATA = sim.csv
BIN  = spclust

LIB_A    = libspclust.a
# Shared library for FFI (the ctypes test needs it). On macOS -shared also produces a
# Mach-O dylib, so the name is uniformly .so and both can dlopen it.
SHARED_A = libspclust.so
CUT_TEST = cut_property
C_API_TEST = c_api_test

.PHONY: all lib lib-shared run test test-cut test-scipy test-c test-py check asan asan-cases leaks bench bench-update clean help

all: $(BIN)

$(BIN): $(ALL_SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) $(INCS) $(DEFS) $(LDFLAGS) -o $@ $(ALL_SRCS) $(LDLIBS)

# Build only the static library (shows the library layer is independent of tools/)
$(LIB_A): $(LIB_SRCS) $(HDRS)
	@rm -rf .libobj && mkdir -p .libobj
	$(CXX) $(CXXFLAGS) $(INCS) $(DEFS) -c $(LIB_SRCS)
	@mv *.o .libobj/
	ar rcs $@ .libobj/*.o
	rm -rf .libobj

lib: $(LIB_A)

$(SHARED_A): $(LIB_SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) -fPIC -shared $(INCS) $(DEFS) -o $@ $(LIB_SRCS) $(LDLIBS)

lib-shared: $(SHARED_A)

run: $(BIN)
	./$(BIN) $(DATA) --linkage $(LINKAGE)

# Differential test (CLI side): an independent oracle vs spclust merge by merge; one binary, three linkages
# Property test (library side): the three cut entry points, expectations derived from the tree bottom-up
test: $(BIN) test-cut test-c test-py
	$(PYTHON) tests/run_tests.py --bin ./$(BIN) --min-sim $(MIN_SIM)
	$(PYTHON) tests/test_cli.py --bin ./$(BIN)
	$(SCIPY_PYTHON) tests/test_scipy.py --bin ./$(BIN)

$(CUT_TEST): tests/test_cut.cpp $(LIB_SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) $(INCS) -o $@ tests/test_cut.cpp $(LIB_SRCS) $(LDLIBS)

test-cut: $(CUT_TEST)
	./$(CUT_TEST)

# scipy interop: hand the --scipy export to the real scipy (skips and returns 0 without scipy)
SCIPY_PYTHON ?= $(PYTHON)
test-scipy: $(BIN)
	$(SCIPY_PYTHON) tests/test_scipy.py --bin ./$(BIN)

# Python binding (ctypes over the C API). Needs numpy for the array path; the test skips
# itself and returns 0 when numpy or the shared library is missing.
test-py: $(SHARED_A)
	PYTHONPATH=python $(PYTHON) tests/test_python.py

# C API: compile the test with a **C compiler** (proving the header is plain C), then call
# it from Python via ctypes (dlopen, symbol names, struct layout, the exception boundary).
test-c: $(LIB_A) $(SHARED_A)
	$(CC) -std=c11 -Wall -Wextra -Iinclude -c tests/test_c_api.c -o .c_api.o
	$(CXX) $(CXXFLAGS) -o $(C_API_TEST) .c_api.o $(LIB_A) $(LDLIBS)
	rm -f .c_api.o
	./$(C_API_TEST)
	$(PYTHON) tests/test_c_api.py --lib ./$(SHARED_A)

# Performance guard: inputs from the fixed generator, alternating repetitions, median vs baseline.
# Over tolerance only warns; only 3x time / 2x memory fails (machine noise makes hard limits unreliable).
bench: $(BIN)
	$(PYTHON) tests/bench.py --bin ./$(BIN) $(BENCH_ARGS)

# Re-record the baseline after an intentional performance change (and on a new machine)
bench-update: $(BIN)
	$(PYTHON) tests/bench.py --bin ./$(BIN) --update

# Run sim.csv once per linkage for eyeball comparison
check: $(BIN)
	@for m in 0 1 2; do \
		echo "=============== --linkage $$m ==============="; \
		./$(BIN) $(DATA) --linkage $$m | grep -vE "^(Reading|Done\.|Initializing|Clustering)"; \
	done

# The core is hand-written malloc/free and four-way pointers, so ASan is worth running
# Note: ASan on macOS does not support detect_leaks (LeakSanitizer is Linux-only),
#       use `make leaks` below for leak checking (the leaks tool shipped with macOS)
asan: $(ALL_SRCS) $(HDRS)
	$(CXX) -std=c++17 -g -O1 -fno-omit-frame-pointer -Wall -Wextra $(INCS) \
		-fsanitize=address,undefined $(DEFS) -o $(BIN)-asan $(ALL_SRCS) $(LDLIBS)
	./$(BIN)-asan $(DATA)

# Run every test case under ASan/UBSan (duplicate pairs, self-loops, threshold edges)
asan-cases: $(ALL_SRCS) $(HDRS)
	$(CXX) -std=c++17 -g -O1 -fno-omit-frame-pointer -Wall -Wextra $(INCS) \
		-fsanitize=address,undefined $(DEFS) -o $(BIN)-asan $(ALL_SRCS) $(LDLIBS)
	@fail=0; \
	for f in tests/cases/*.csv; do \
		out=$$(./$(BIN)-asan $$f 2>&1) || { echo "non-zero exit code: $$f"; echo "$$out" | head -20; fail=1; }; \
		case "$$out" in \
			*AddressSanitizer*|*"runtime error"*) \
				echo "ASan/UBSan reported a problem: $$f"; \
				echo "$$out" | grep -m3 -E "AddressSanitizer|runtime error"; \
				fail=1;; \
		esac; \
	done; \
	if [ $$fail -eq 0 ]; then echo "ASan/UBSan: all $$(ls tests/cases/*.csv | wc -l | tr -d ' ') cases passed"; else exit 1; fi

# Leak detection shipped with macOS (not Linux's LSan).
# Note: leaks occasionally reports "Process ... is not debuggable" once the process has exited; rerun.
leaks: $(BIN)
	MallocStackLogging=1 leaks --atExit -- ./$(BIN) $(DATA)

clean:
	rm -f $(BIN) $(BIN)-asan $(LIB_A) $(SHARED_A) $(CUT_TEST) $(C_API_TEST) *.o
	rm -rf .libobj
	rm -rf $(BIN)*.dSYM
	find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true

help:
	@echo "make [all|run|test|test-cut|test-scipy|test-c|lib|lib-shared|bench|check|asan|asan-cases|leaks|clean|help]"
	@echo "  LINKAGE=$(LINKAGE)  MIN_SIM=$(MIN_SIM)  DEBUG=$(DEBUG)"
	@echo "  e.g. make run LINKAGE=1 MIN_SIM=60"
	@echo "       make test MIN_SIM=70"
	@echo "       python3 tests/run_tests.py --only average --verbose"
	@echo "       make bench BENCH_ARGS=--quick     # performance guard (skip the largest cases)"

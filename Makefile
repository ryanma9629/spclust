# spclust —— 构建
#
#   make            构建 ./spclust
#   make run        构建并跑 sim.csv
#   make test       差分测试 + 切树属性测试 + CLI 测试 + scipy 互操作
#   make check      用三种 linkage 各跑一遍 sim.csv 做对比
#   make asan       AddressSanitizer + UBSan 构建并运行（强烈建议跑一次）
#   make asan-cases 用 ASan/UBSan 跑一遍 tests/cases 下全部用例
#   make leaks      macOS 自带 leaks 工具查内存泄漏
#   make clean
#   make help
#
# 可覆盖的变量：
#   make CXX=g++-15                    换编译器
#   make run LINKAGE=1 MIN_SIM=60      LINKAGE 现在传给 --linkage（运行时选，不用重编译）
#   make test MIN_SIM=70               用另一个阈值重跑差分测试（oracle 会同步）
#   make CXXFLAGS="-std=c++11 -O2"     完全接管编译选项
#
# LINKAGE: 0=average  1=max-similarity（= 距离语义的单链接，走 Kruskal）  2=min-nonzero
#   链接方式已经是**运行时**的（spclust::Linkage），所以这个变量只是传给命令行工具，
#   不再需要为每种链接各编一个二进制。
#
# 注：MIN_SIM_CLUST / DEBUG 是 CLI 自己的编译期默认值（tools/cli_config.h 里用 #ifndef
#     包住），命令行 -D 优先；库本体没有任何编译期选项。

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

# 目录结构：库（src/ + include/）与使用者（tools/）分开。
# 构建系统见 CMakeLists.txt；这个 Makefile 是开发时的快捷方式。
LIB_SRCS = src/hac.cpp src/scipy.cpp
CLI_SRCS = tools/spclust_cli.cpp tools/readdb.cpp
ALL_SRCS = $(CLI_SRCS) $(LIB_SRCS)
HDRS = include/spclust/hac.h tools/readdb.h tools/cli_config.h
INCS = -Iinclude -Itools
DATA = sim.csv
BIN  = spclust

LIB_A    = libspclust.a
CUT_TEST = cut_property

.PHONY: all lib run test test-cut test-scipy check asan asan-cases leaks bench bench-update clean help

all: $(BIN)

$(BIN): $(ALL_SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) $(INCS) $(DEFS) $(LDFLAGS) -o $@ $(ALL_SRCS) $(LDLIBS)

# 只把库编成静态库（演示"库"这一层是独立的、不依赖 tools/）
$(LIB_A): $(LIB_SRCS) $(HDRS)
	@rm -rf .libobj && mkdir -p .libobj
	$(CXX) $(CXXFLAGS) $(INCS) $(DEFS) -c $(LIB_SRCS)
	@mv *.o .libobj/
	ar rcs $@ .libobj/*.o
	rm -rf .libobj

lib: $(LIB_A)

run: $(BIN)
	./$(BIN) $(DATA) --linkage $(LINKAGE)

# 差分测试（CLI 侧）：独立 oracle 与 spclust 逐次合并比对，一个二进制覆盖三种 linkage
# 属性测试（库侧）：三个切树接口，期望值从合并树自底向上独立算出
test: $(BIN) test-cut
	$(PYTHON) tests/run_tests.py --bin ./$(BIN) --min-sim $(MIN_SIM)
	$(PYTHON) tests/test_cli.py --bin ./$(BIN)
	$(SCIPY_PYTHON) tests/test_scipy.py --bin ./$(BIN)

$(CUT_TEST): tests/test_cut.cpp $(LIB_SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) $(INCS) -o $@ tests/test_cut.cpp $(LIB_SRCS) $(LDLIBS)

test-cut: $(CUT_TEST)
	./$(CUT_TEST)

# scipy 互操作：把 --scipy 的导出交给真的 scipy 读（没有 scipy 时自动跳过、返回 0）
SCIPY_PYTHON ?= $(PYTHON)
test-scipy: $(BIN)
	$(SCIPY_PYTHON) tests/test_scipy.py --bin ./$(BIN)

# 性能回归护栏：固定生成器的输入 + 交替重复取中位，与 tests/bench_baseline.txt 比。
# 超出容差只警告，超过 3x 时间 / 2x 内存才失败（机器噪声大，硬判定会误报）。
bench: $(BIN)
	$(PYTHON) tests/bench.py --bin ./$(BIN) $(BENCH_ARGS)

# 有意改变性能之后重新记录基线（换机器也要重记）
bench-update: $(BIN)
	$(PYTHON) tests/bench.py --bin ./$(BIN) --update

# 三种 linkage 各跑一遍 sim.csv，用于肉眼对比
check: $(BIN)
	@for m in 0 1 2; do \
		echo "=============== --linkage $$m ==============="; \
		./$(BIN) $(DATA) --linkage $$m | grep -vE "^(Reading|Done\.|Initializing|Clustering)"; \
	done

# spclust.cpp 里全是手写 malloc/free 和四向指针，asan 很值得跑
# 注意：macOS 的 ASan 不支持 detect_leaks（LeakSanitizer 只在 Linux 上可用），
#       查内存泄漏请用下面的 make leaks（macOS 自带的 leaks 工具）
asan: $(ALL_SRCS) $(HDRS)
	$(CXX) -std=c++17 -g -O1 -fno-omit-frame-pointer -Wall -Wextra $(INCS) \
		-fsanitize=address,undefined $(DEFS) -o $(BIN)-asan $(ALL_SRCS) $(LDLIBS)
	./$(BIN)-asan $(DATA)

# 用 ASan/UBSan 跑一遍全部测试用例（含重复 pair、自环、阈值边界）
asan-cases: $(ALL_SRCS) $(HDRS)
	$(CXX) -std=c++17 -g -O1 -fno-omit-frame-pointer -Wall -Wextra $(INCS) \
		-fsanitize=address,undefined $(DEFS) -o $(BIN)-asan $(ALL_SRCS) $(LDLIBS)
	@fail=0; \
	for f in tests/cases/*.csv; do \
		out=$$(./$(BIN)-asan $$f 2>&1) || { echo "退出码非 0: $$f"; echo "$$out" | head -20; fail=1; }; \
		case "$$out" in \
			*AddressSanitizer*|*"runtime error"*) \
				echo "ASan/UBSan 报告问题: $$f"; \
				echo "$$out" | grep -m3 -E "AddressSanitizer|runtime error"; \
				fail=1;; \
		esac; \
	done; \
	if [ $$fail -eq 0 ]; then echo "ASan/UBSan: 全部 $$(ls tests/cases/*.csv | wc -l | tr -d ' ') 个用例通过"; else exit 1; fi

# macOS 自带的泄漏检测（不是 Linux 的 LSan）。
# 注：leaks 偶尔会因为进程已经退出而报 "Process ... is not debuggable"，多跑一次即可。
leaks: $(BIN)
	MallocStackLogging=1 leaks --atExit -- ./$(BIN) $(DATA)

clean:
	rm -f $(BIN) $(BIN)-asan $(LIB_A) $(CUT_TEST) *.o
	rm -rf .libobj
	rm -rf $(BIN)*.dSYM
	find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true

help:
	@echo "make [all|run|test|test-cut|test-scipy|bench|check|asan|asan-cases|leaks|clean|help]"
	@echo "  LINKAGE=$(LINKAGE)  MIN_SIM=$(MIN_SIM)  DEBUG=$(DEBUG)"
	@echo "  e.g. make run LINKAGE=1 MIN_SIM=60"
	@echo "       make test MIN_SIM=70"
	@echo "       python3 tests/run_tests.py --only average --verbose"
	@echo "       make bench BENCH_ARGS=--quick     # 性能护栏（跳过最大用例）"

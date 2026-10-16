# spclust —— macOS / POSIX 构建
#
#   make            构建 ./spclust
#   make run        构建并跑 sim.csv
#   make test       差分测试：用独立 oracle 交叉验证 spclust 的正确性
#   make check      用三种 linkage 各跑一遍 sim.csv 做对比
#   make asan       AddressSanitizer + UBSan 构建并运行（强烈建议跑一次）
#   make asan-cases 用 ASan/UBSan 跑一遍 tests/cases 下全部用例
#   make leaks      macOS 自带 leaks 工具查内存泄漏
#   make clean
#   make help
#
# 可覆盖的变量：
#   make CXX=g++-15                    换编译器
#   make run LINKAGE=1 MIN_SIM=60      覆盖 main.h 里的同名宏（无需改源码）
#   make test MIN_SIM=70               用另一个阈值重跑差分测试（oracle 会同步）
#   make CXXFLAGS="-std=c++11 -O2"     完全接管编译选项
#
# LINKAGE: 0=average  1=max（原注释 complete）  2=min（原注释 single，仅当两条边都在时才更新）
#
# 注：LINKAGE_MTD / MIN_SIM_CLUST / DEBUG 通过命令行 -D 传入，
#     main.h 里对它们用了 #ifndef，所以命令行优先。

CXX      = c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
LDFLAGS  ?=
LDLIBS   ?=
PYTHON   ?= python3
BENCH_ARGS ?=

LINKAGE ?= 0
MIN_SIM ?= 50
DEBUG   ?= 1

DEFS = -DLINKAGE_MTD=$(LINKAGE) -DMIN_SIM_CLUST=$(MIN_SIM) -DDEBUG=$(DEBUG)

SRCS = main.cpp readdb.cpp spclust.cpp
HDRS = main.h readdb.h spclust.h
DATA = sim.csv
BIN  = spclust

# 测试用的三种 linkage 二进制：LINKAGE_MTD 是编译期常量，只能各编一个。
# 它们固定用 DEBUG=2，因为差分测试需要解析 "Group i and group j are merged..." 日志。
# 二进制名里带上 MIN_SIM，这样 `make test MIN_SIM=70` 不会误用 50 编出来的旧二进制。
LBIN        = $(BIN)-linkage-$(MIN_SIM)
TRACE_DEBUG = 2

.PHONY: all run test check linkages asan asan-cases leaks bench bench-update clean help

all: $(BIN)

$(BIN): $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) $(DEFS) $(LDFLAGS) -o $@ $(SRCS) $(LDLIBS)

$(LBIN)-%: $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) -DLINKAGE_MTD=$* -DMIN_SIM_CLUST=$(MIN_SIM) -DDEBUG=$(TRACE_DEBUG) \
		$(LDFLAGS) -o $@ $(SRCS) $(LDLIBS)

linkages: $(LBIN)-0 $(LBIN)-1 $(LBIN)-2

run: $(BIN)
	./$(BIN) $(DATA)

# 差分测试：独立 oracle 与 spclust 逐次合并比对，覆盖三种 linkage
test: $(BIN) linkages
	$(PYTHON) tests/run_tests.py --bin-prefix ./$(LBIN)- --min-sim $(MIN_SIM)
	$(PYTHON) tests/test_cli.py --bin ./$(BIN)

# 性能回归护栏：固定生成器的输入 + 交替重复取中位，与 tests/bench_baseline.txt 比。
# 超出容差只警告，超过 3x 时间 / 2x 内存才失败（机器噪声大，硬判定会误报）。
bench: $(BIN) linkages
	$(PYTHON) tests/bench.py --bin ./$(BIN) --bin-linkage ./$(LBIN)-1 $(BENCH_ARGS)

# 有意改变性能之后重新记录基线（换机器也要重记）
bench-update: $(BIN) linkages
	$(PYTHON) tests/bench.py --bin ./$(BIN) --bin-linkage ./$(LBIN)-1 --update

# 三种 linkage 各跑一遍 sim.csv，用于肉眼对比
check: linkages
	@for m in 0 1 2; do \
		echo "=============== LINKAGE_MTD=$$m ==============="; \
		./$(LBIN)-$$m $(DATA) | grep -vE "^(Reading|Done\.|Initializing|Clustering)"; \
	done

# spclust.cpp 里全是手写 malloc/free 和四向指针，asan 很值得跑
# 注意：macOS 的 ASan 不支持 detect_leaks（LeakSanitizer 只在 Linux 上可用），
#       查内存泄漏请用下面的 make leaks（macOS 自带的 leaks 工具）
asan: $(SRCS) $(HDRS)
	$(CXX) -std=c++17 -g -O1 -fno-omit-frame-pointer -Wall -Wextra \
		-fsanitize=address,undefined $(DEFS) -o $(BIN)-asan $(SRCS) $(LDLIBS)
	./$(BIN)-asan $(DATA)

# 用 ASan/UBSan 跑一遍全部测试用例（含重复 pair、自环、阈值边界）
asan-cases: $(SRCS) $(HDRS)
	$(CXX) -std=c++17 -g -O1 -fno-omit-frame-pointer -Wall -Wextra \
		-fsanitize=address,undefined $(DEFS) -o $(BIN)-asan $(SRCS) $(LDLIBS)
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
	rm -f $(BIN) $(BIN)-asan $(BIN)-linkage-*
	rm -rf $(BIN)*.dSYM
	find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true

help:
	@echo "make [all|run|test|bench|check|linkages|asan|asan-cases|leaks|clean|help]"
	@echo "  LINKAGE=$(LINKAGE)  MIN_SIM=$(MIN_SIM)  DEBUG=$(DEBUG)"
	@echo "  e.g. make run LINKAGE=1 MIN_SIM=60"
	@echo "       make test MIN_SIM=70"
	@echo "       python3 tests/run_tests.py --only average --verbose"
	@echo "       make bench BENCH_ARGS=--quick     # 性能护栏（跳过最大用例）"

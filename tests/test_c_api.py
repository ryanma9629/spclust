#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
#
# C API 的 FFI 测试：用 ctypes 从 Python 调它。
#
# 为什么除了 test_c_api.c（真 C 编译器）还要这一层：C 接口的主要用途就是给**别的运行时**
# 调（Rust / Go / Python / Julia...）。真 C 那一层验证的是头文件与 ABI；这一层验证的是
# 共享库能被动态加载、extern "C" 的符号名在外部运行时里找得到、指针/结构体布局对得上、
# 并且**跨语言调用不会因为异常而炸掉进程**（比如自环必须返回错误码，而不是 terminate）。
#
# 需要一个共享库：make test-c 会顺手编一个（见 Makefile）。
#   python3 tests/test_c_api.py --lib /tmp/libspclust.so

import argparse
import ctypes
import sys


class Merge(ctypes.Structure):
    _fields_ = [("left", ctypes.c_uint32),
                ("right", ctypes.c_uint32),
                ("similarity", ctypes.c_double),
                ("cohesion", ctypes.c_double),
                ("size", ctypes.c_size_t)]


class ScipyRow(ctypes.Structure):
    _fields_ = [("left", ctypes.c_uint32),
                ("right", ctypes.c_uint32),
                ("distance", ctypes.c_double),
                ("size", ctypes.c_uint32)]


EDGES = [(0, 1, 91.0), (2, 3, 88.0), (4, 5, 85.0),
         (0, 2, 31.0), (1, 3, 22.0), (0, 4, 13.0), (2, 5, 7.0), (1, 5, 3.0)]
N = 6

OK = 0
ERR_INVALID_ARGUMENT = 1


def main():
    ap = argparse.ArgumentParser(description="C API 的 ctypes/FFI 测试")
    ap.add_argument("--lib", required=True, help="共享库路径（libspclust.so/.dylib）")
    args = ap.parse_args()

    lib = ctypes.CDLL(args.lib)
    pass_ = 0
    fail = 0

    def check(ok, what):
        nonlocal pass_, fail
        if ok:
            pass_ += 1
        else:
            fail += 1
            print(f"  ✗ {what}")

    # ---- 签名 ----
    lib.spclust_version.restype = ctypes.c_char_p
    lib.spclust_last_error.restype = ctypes.c_char_p
    lib.spclust_linkage_from_name.argtypes = [ctypes.c_char_p]
    lib.spclust_linkage_from_name.restype = ctypes.c_int
    lib.spclust_hac_new.argtypes = [ctypes.c_int]
    lib.spclust_hac_new.restype = ctypes.c_void_p
    lib.spclust_hac_free.argtypes = [ctypes.c_void_p]
    lib.spclust_hac_add_edge.argtypes = [ctypes.c_void_p, ctypes.c_uint32,
                                        ctypes.c_uint32, ctypes.c_double]
    lib.spclust_hac_add_edge.restype = ctypes.c_int
    lib.spclust_hac_build.argtypes = [ctypes.c_void_p, ctypes.c_double]
    lib.spclust_hac_build.restype = ctypes.c_int
    lib.spclust_hac_num_objects.argtypes = [ctypes.c_void_p]
    lib.spclust_hac_num_objects.restype = ctypes.c_size_t
    lib.spclust_hac_tree_size.argtypes = [ctypes.c_void_p]
    lib.spclust_hac_tree_size.restype = ctypes.c_size_t
    lib.spclust_hac_tree.argtypes = [ctypes.c_void_p, ctypes.POINTER(Merge),
                                    ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]
    lib.spclust_hac_tree.restype = ctypes.c_int
    lib.spclust_hac_cut.argtypes = [ctypes.c_void_p, ctypes.c_double]
    lib.spclust_hac_cut.restype = ctypes.c_void_p
    lib.spclust_hac_cut_n_clusters.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    lib.spclust_hac_cut_n_clusters.restype = ctypes.c_void_p
    lib.spclust_cut_free.argtypes = [ctypes.c_void_p]
    for fn in ("spclust_cut_n_groups", "spclust_cut_n_members"):
        getattr(lib, fn).argtypes = [ctypes.c_void_p]
        getattr(lib, fn).restype = ctypes.c_size_t
    lib.spclust_cut_members.argtypes = [ctypes.c_void_p]
    lib.spclust_cut_members.restype = ctypes.POINTER(ctypes.c_uint32)
    lib.spclust_cut_offsets.argtypes = [ctypes.c_void_p]
    lib.spclust_cut_offsets.restype = ctypes.POINTER(ctypes.c_size_t)
    lib.spclust_hac_to_scipy.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_double,
        ctypes.POINTER(ScipyRow), ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t)]
    lib.spclust_hac_to_scipy.restype = ctypes.c_int

    check(lib.spclust_version().decode() == "0.3.2", "版本字符串")
    check(lib.spclust_linkage_from_name(b"single") == 1, "名字解析")

    # ---- 非法 linkage：必须返回 NULL 而不是抛异常穿出 C 边界 ----
    check(not lib.spclust_hac_new(7), "非法 linkage -> NULL")
    check(lib.spclust_last_error() is not None, "last_error 有内容")

    h = lib.spclust_hac_new(0)
    check(bool(h), "建 handle")
    for a, b, s in EDGES:
        check(lib.spclust_hac_add_edge(h, a, b, s) == OK, f"add_edge({a},{b})")
    check(lib.spclust_hac_build(h, 0.0) == OK, "build")
    check(lib.spclust_hac_num_objects(h) == N, "num_objects")

    # **关键**：自环在 C++ 里是抛异常，从 Python 调必须拿到错误码。
    # 如果异常真的穿过去了，这里会直接 abort 掉整个进程（测试就崩了）。
    check(lib.spclust_hac_add_edge(h, 2, 2, 1.0) == ERR_INVALID_ARGUMENT,
          "自环 -> 错误码（异常没有穿过 C 边界）")
    check(b"NUL" not in (lib.spclust_last_error() or b""), "last_error 是 C 字符串")

    # ---- 树：两段式（先问长度再填缓冲）----
    n = lib.spclust_hac_tree_size(h)
    check(n == N - 1, "树有 n-1 行")
    got = ctypes.c_size_t(0)
    check(lib.spclust_hac_tree(h, None, 0, ctypes.byref(got)) == OK, "只问行数")
    check(got.value == n, "问出来的长度对")
    rows = (Merge * n)()
    check(lib.spclust_hac_tree(h, rows, n, ctypes.byref(got)) == OK, "取树")
    check(rows[0].size == 2 and rows[-1].size == N, "第一/最后一行的 size")
    check(abs(rows[0].similarity - 91.0) < 1e-9, "第一次合并的相似度（结构体布局对）")

    # ---- 切树：offsets/members 的视图 ----
    c = lib.spclust_hac_cut(h, 50.0)
    check(bool(c), "cut(50)")
    ng = lib.spclust_cut_n_groups(c)
    nm = lib.spclust_cut_n_members(c)
    check(ng == 3 and nm == N, "3 个簇、成员总数 n")
    off = lib.spclust_cut_offsets(c)
    mem = lib.spclust_cut_members(c)
    groups = [sorted(mem[off[g]:off[g + 1]]) for g in range(ng)]
    check(groups == [[0, 1], [2, 3], [4, 5]], f"分组内容 {groups}")
    lib.spclust_cut_free(c)

    c = lib.spclust_hac_cut_n_clusters(h, 2)
    check(lib.spclust_cut_n_groups(c) == 2, "cut_n_clusters(2)")
    lib.spclust_cut_free(c)

    # ---- scipy 导出 ----
    cnt = ctypes.c_size_t(0)
    comps = ctypes.c_size_t(0)
    check(lib.spclust_hac_to_scipy(h, 0, 0.0, 0, float("nan"), None, 0,
                                   ctypes.byref(cnt), ctypes.byref(comps)) == OK,
          "scipy 导出：问行数")
    check(cnt.value == N - 1 and comps.value == 1, "行数 n-1、1 个分量")
    srows = (ScipyRow * cnt.value)()
    check(lib.spclust_hac_to_scipy(h, 0, 0.0, 0, float("nan"), srows, cnt.value,
                                   ctypes.byref(cnt), ctypes.byref(comps)) == OK,
          "scipy 导出：取行")
    check(srows[0].left < srows[0].right, "子节点升序")

    lib.spclust_hac_free(h)
    lib.spclust_hac_free(None)          # 不该崩

    print("--------------------------------------------------------------------")
    print(f"PASS={pass_} FAIL={fail}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

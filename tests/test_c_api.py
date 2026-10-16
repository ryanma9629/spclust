#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
#
# FFI test for the C API: call it from Python via ctypes.

# Why this layer exists on top of test_c_api.c (a real C compiler): the main use of a C
# interface is being called from **other runtimes** (Rust / Go / Python / Julia ...). The C
# layer verifies the header and the ABI; this layer verifies that the shared library can be
# loaded dynamically, that the extern "C" symbol names are findable from a foreign runtime, that
# pointer/struct layout matches, and - most importantly - that **a cross-language call cannot
# kill the process through an exception** (a self-loop must return an error code, not terminate).

# Needs a shared library: make test-c builds one (see the Makefile).
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
    ap = argparse.ArgumentParser(description="ctypes/FFI test for the C API")
    ap.add_argument("--lib", required=True, help="path to the shared library (libspclust.so)")
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

    # ---- signatures ----
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

    check(lib.spclust_version().decode() == "0.3.2", "version string")
    check(lib.spclust_linkage_from_name(b"single") == 1, "name lookup")

    # ---- bad linkage: must return NULL rather than throw through the C boundary ----
    check(not lib.spclust_hac_new(7), "bad linkage -> NULL")
    check(lib.spclust_last_error() is not None, "last_error has content")

    h = lib.spclust_hac_new(0)
    check(bool(h), "create handle")
    for a, b, s in EDGES:
        check(lib.spclust_hac_add_edge(h, a, b, s) == OK, f"add_edge({a},{b})")
    check(lib.spclust_hac_build(h, 0.0) == OK, "build")
    check(lib.spclust_hac_num_objects(h) == N, "num_objects")

    # **The key one**: a self-loop throws in C++, so a call from Python must yield an error code.
    # If the exception really escaped, this would abort the whole process (and the test would die).
    check(lib.spclust_hac_add_edge(h, 2, 2, 1.0) == ERR_INVALID_ARGUMENT,
          "self-loop -> error code (the exception did not cross the C boundary)")
    check(b"NUL" not in (lib.spclust_last_error() or b""), "last_error is a C string")

    # ---- tree: the two-call pattern (ask for the length, then fill the buffer) ----
    n = lib.spclust_hac_tree_size(h)
    check(n == N - 1, "the tree has n-1 rows")
    got = ctypes.c_size_t(0)
    check(lib.spclust_hac_tree(h, None, 0, ctypes.byref(got)) == OK, "ask for the length only")
    check(got.value == n, "the reported length is right")
    rows = (Merge * n)()
    check(lib.spclust_hac_tree(h, rows, n, ctypes.byref(got)) == OK, "fetch the tree")
    check(rows[0].size == 2 and rows[-1].size == N, "sizes of the first and last rows")
    check(abs(rows[0].similarity - 91.0) < 1e-9, "similarity of the first merge (struct layout matches)")

    # ---- cutting: the offsets/members views ----
    c = lib.spclust_hac_cut(h, 50.0)
    check(bool(c), "cut(50)")
    ng = lib.spclust_cut_n_groups(c)
    nm = lib.spclust_cut_n_members(c)
    check(ng == 3 and nm == N, "3 clusters and n members in total")
    off = lib.spclust_cut_offsets(c)
    mem = lib.spclust_cut_members(c)
    groups = [sorted(mem[off[g]:off[g + 1]]) for g in range(ng)]
    check(groups == [[0, 1], [2, 3], [4, 5]], f"cluster contents {groups}")
    lib.spclust_cut_free(c)

    c = lib.spclust_hac_cut_n_clusters(h, 2)
    check(lib.spclust_cut_n_groups(c) == 2, "cut_n_clusters(2)")
    lib.spclust_cut_free(c)

    # ---- scipy export ----
    cnt = ctypes.c_size_t(0)
    comps = ctypes.c_size_t(0)
    check(lib.spclust_hac_to_scipy(h, 0, 0.0, 0, float("nan"), None, 0,
                                   ctypes.byref(cnt), ctypes.byref(comps)) == OK,
          "scipy export: ask for the row count")
    check(cnt.value == N - 1 and comps.value == 1, "n-1 rows and 1 component")
    srows = (ScipyRow * cnt.value)()
    check(lib.spclust_hac_to_scipy(h, 0, 0.0, 0, float("nan"), srows, cnt.value,
                                   ctypes.byref(cnt), ctypes.byref(comps)) == OK,
          "scipy export: fetch the rows")
    check(srows[0].left < srows[0].right, "children are ascending")

    lib.spclust_hac_free(h)
    lib.spclust_hac_free(None)          # must not crash

    print("--------------------------------------------------------------------")
    print(f"PASS={pass_} FAIL={fail}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

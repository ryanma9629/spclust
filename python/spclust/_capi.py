# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
"""ctypes binding of the spclust C API (include/spclust/spclust_c.h).

Kept separate from the friendly API in ``__init__`` so that the raw layer stays a
mechanical translation of the C header and is easy to audit against it.

Nothing here holds the GIL during a call: ctypes releases it around every foreign
call by default, which is why ``build()`` on a large input does not block other
Python threads.
"""

import ctypes
import os
import sys
from ctypes.util import find_library

__all__ = ["load_library", "Lib", "SpclustError", "STATUS_NAMES"]


class SpclustError(Exception):
    """Raised for any non-OK spclust_status, carrying the library's message."""


# status codes, mirroring the enum in the header
STATUS_NAMES = {
    1: "invalid argument",
    2: "out of memory",
    3: "bad state",
    4: "buffer too small",
    5: "unknown error",
}

# linkage ids
AVERAGE = 0
MAX_SIMILARITY = 1
MIN_NONZERO = 2

# distance mappings for the scipy export
COMPLEMENT = 0
NEGATE = 1
RECIPROCAL = 2


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


def _candidate_paths(explicit):
    """Where to look for the shared library, in order.

    An explicit path always wins; then $SPCLUST_LIBRARY; then a copy shipped inside
    this package (which is what an installed wheel would contain); then the build
    directories of a source checkout; then the system loader.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    if explicit:
        yield explicit
    env = os.environ.get("SPCLUST_LIBRARY")
    if env:
        yield env
    for name in ("libspclust.so", "libspclust.dylib", "spclust.dll", "libspclust.dll"):
        yield os.path.join(here, name)
    root = os.path.dirname(os.path.dirname(here))          # the repository root
    for sub in ("build", "build-libonly", "cmake-build-release", "."):
        for name in ("libspclust.so", "libspclust.dylib"):
            yield os.path.join(root, sub, name)
    found = find_library("spclust")
    if found:
        yield found


def load_library(path=None):
    """Load libspclust and return a ready-to-use :class:`Lib`.

    Raises :class:`SpclustError` listing every location that was tried, because
    "cannot find the shared library" is the first thing that goes wrong and a bare
    ctypes message does not say where to look.
    """
    tried = []
    for cand in _candidate_paths(path):
        tried.append(cand)
        if not os.path.exists(cand):
            continue
        try:
            return Lib(ctypes.CDLL(cand))
        except OSError as exc:                              # pragma: no cover
            raise SpclustError(f"failed to load {cand}: {exc}") from exc
    raise SpclustError(
        "cannot find the spclust shared library. Tried:\n  "
        + "\n  ".join(tried)
        + "\n\nBuild it with `make lib-shared` (or `cmake --build build`), or set "
          "SPCLUST_LIBRARY=/path/to/libspclust.so"
    )


class Lib:
    """Thin typed wrapper over the C API. One instance is enough per process."""

    def __init__(self, cdll):
        self.dll = cdll
        d = cdll
        d.spclust_version.restype = ctypes.c_char_p
        d.spclust_last_error.restype = ctypes.c_char_p
        d.spclust_linkage_from_name.argtypes = [ctypes.c_char_p]
        d.spclust_linkage_from_name.restype = ctypes.c_int
        d.spclust_hac_new.argtypes = [ctypes.c_int]
        d.spclust_hac_new.restype = ctypes.c_void_p
        d.spclust_hac_free.argtypes = [ctypes.c_void_p]
        d.spclust_hac_reserve.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_size_t]
        d.spclust_hac_reserve.restype = ctypes.c_int
        d.spclust_hac_add_edge.argtypes = [ctypes.c_void_p, ctypes.c_uint32,
                                           ctypes.c_uint32, ctypes.c_double]
        d.spclust_hac_add_edge.restype = ctypes.c_int
        d.spclust_hac_add_edges.argtypes = [ctypes.c_void_p,
                                            ctypes.POINTER(ctypes.c_uint32),
                                            ctypes.POINTER(ctypes.c_uint32),
                                            ctypes.POINTER(ctypes.c_double),
                                            ctypes.c_size_t]
        d.spclust_hac_add_edges.restype = ctypes.c_int
        d.spclust_hac_build.argtypes = [ctypes.c_void_p, ctypes.c_double]
        d.spclust_hac_build.restype = ctypes.c_int
        d.spclust_hac_reset.argtypes = [ctypes.c_void_p]
        d.spclust_hac_reset.restype = ctypes.c_int
        d.spclust_hac_num_objects.argtypes = [ctypes.c_void_p]
        d.spclust_hac_num_objects.restype = ctypes.c_size_t
        d.spclust_hac_tree_size.argtypes = [ctypes.c_void_p]
        d.spclust_hac_tree_size.restype = ctypes.c_size_t
        d.spclust_hac_linkage.argtypes = [ctypes.c_void_p]
        d.spclust_hac_linkage.restype = ctypes.c_int
        d.spclust_hac_linkage_name.argtypes = [ctypes.c_void_p]
        d.spclust_hac_linkage_name.restype = ctypes.c_char_p
        d.spclust_hac_tree.argtypes = [ctypes.c_void_p, ctypes.POINTER(Merge),
                                       ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]
        d.spclust_hac_tree.restype = ctypes.c_int
        d.spclust_hac_cut.argtypes = [ctypes.c_void_p, ctypes.c_double]
        d.spclust_hac_cut_n_clusters.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        d.spclust_hac_cut_min_size.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        for fn in (d.spclust_hac_cut, d.spclust_hac_cut_n_clusters, d.spclust_hac_cut_min_size):
            fn.restype = ctypes.c_void_p
        d.spclust_cut_free.argtypes = [ctypes.c_void_p]
        d.spclust_cut_n_groups.argtypes = [ctypes.c_void_p]
        d.spclust_cut_n_groups.restype = ctypes.c_size_t
        d.spclust_cut_n_members.argtypes = [ctypes.c_void_p]
        d.spclust_cut_n_members.restype = ctypes.c_size_t
        d.spclust_cut_members.argtypes = [ctypes.c_void_p]
        d.spclust_cut_members.restype = ctypes.POINTER(ctypes.c_uint32)
        d.spclust_cut_offsets.argtypes = [ctypes.c_void_p]
        d.spclust_cut_offsets.restype = ctypes.POINTER(ctypes.c_size_t)
        d.spclust_cut_cohesion.argtypes = [ctypes.c_void_p]
        d.spclust_cut_cohesion.restype = ctypes.POINTER(ctypes.c_double)
        d.spclust_hac_to_scipy.argtypes = [
            ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_double,
            ctypes.POINTER(ScipyRow), ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t)]
        d.spclust_hac_to_scipy.restype = ctypes.c_int

    # -- error handling ---------------------------------------------------
    def last_error(self):
        msg = self.dll.spclust_last_error()
        return msg.decode("utf-8", "replace") if msg else None

    def check(self, status, what="call"):
        if status != 0:
            detail = self.last_error() or STATUS_NAMES.get(status, f"status {status}")
            raise SpclustError(f"{what}: {detail}")

    # -- misc -------------------------------------------------------------
    def version(self):
        return self.dll.spclust_version().decode()

    def linkage_from_name(self, name):
        return self.dll.spclust_linkage_from_name(name.encode())


_lib = None


def get_library(path=None):
    """Process-wide lazily loaded library (so importing this module is cheap)."""
    global _lib
    if _lib is None:
        _lib = load_library(path)
    return _lib

# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
"""Agglomerative hierarchical clustering over a sparse similarity matrix.

This is the Python face of the C++ library: thin, and deliberately shaped like the
rest of the ecosystem (``scipy.cluster.hierarchy``, scikit-learn) so it drops into an
existing pipeline.

    import numpy as np, spclust
    h = spclust.HAC(linkage="average")
    h.add_edges(rows, cols, sims)      # sparse: only the edges you have
    h.build()                          # full merge tree
    labels = h.cut(50.0).labels_       # cut wherever you like, no re-clustering
    Z = h.to_scipy()                   # scipy linkage matrix, for fcluster/dendrogram

Why this exists: ``scipy.cluster.hierarchy.linkage`` needs a dense distance matrix, so
its memory is O(n^2) before it does any work - about 9.7 GB at n = 36000, and ~320 GB
at n = 10^5. This library stores only the edges, so memory is O(nnz): the same
clustering at n = 12800 measured 18.7 MB against scipy's 5005 MB, and 22 ms against
1.29 s.

Semantics (see include/spclust/hac.h for the full contract):
  * object ids must be **dense 0..n-1**; an absent pair counts as similarity 0;
  * a repeated pair takes the **largest** value;
  * self-loops, negative values and NaN raise :class:`SpclustError`;
  * **ties are broken by id**, not by insertion order: among pairs tied on similarity
    the smaller ``(left, right)`` id pair wins, so the same edges with the same ids
    give the same tree in any order - reproducible across machines and versions;
  * one :class:`HAC` is one independent state; separate objects never interact, and
    while a call is running the GIL is released, but a single object is not safe to
    use from two threads at once.
"""

import ctypes
import os
import sys

from ._capi import (AVERAGE, COMPLEMENT, MAX_SIMILARITY, MIN_NONZERO, NEGATE, RECIPROCAL,
                    Lib, Merge, ScipyRow, SpclustError, get_library, load_library)

try:                                    # numpy is optional: everything also works with
    import numpy as _np                 # plain lists / array.array
except ImportError:                     # pragma: no cover
    _np = None

__all__ = [
    "HAC", "Cut", "SpclustError", "load_library", "library_path", "version",
    "AVERAGE", "MAX_SIMILARITY", "MIN_NONZERO",
    "COMPLEMENT", "NEGATE", "RECIPROCAL",
    "LINKAGES", "LINKAGE_NAMES",
]

LINKAGE_NAMES = {AVERAGE: "average", MAX_SIMILARITY: "max-similarity", MIN_NONZERO: "min-nonzero"}
LINKAGES = LINKAGE_NAMES


def version():
    """Version string of the loaded shared library, e.g. ``"0.3.3"``."""
    return get_library().version()


def library_path():
    """Path of the loaded shared library (useful when debugging a mismatch)."""
    return get_library().dll._name


def _resolve_linkage(linkage):
    if isinstance(linkage, str):
        lib = get_library()
        idx = lib.linkage_from_name(linkage)
        if idx < 0:
            raise ValueError(f"unknown linkage {linkage!r}; use one of "
                             f"{sorted(LINKAGE_NAMES.values())}")
        return idx
    if linkage in LINKAGE_NAMES:
        return int(linkage)
    raise ValueError(f"unknown linkage {linkage!r}; use one of {sorted(LINKAGE_NAMES.values())}")


def _as_u32(values):
    """Return (buffer, length) for an integer array usable as uint32."""
    if _np is not None and isinstance(values, _np.ndarray):
        arr = _np.ascontiguousarray(values)
        if arr.dtype != _np.uint32:
            if arr.dtype.kind not in "iu":
                raise TypeError("ids must be integers")
            arr = arr.astype(_np.uint32, copy=False)
        return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)), len(arr), arr
    buf = (ctypes.c_uint32 * len(values))(*values)
    return buf, len(values), buf


def _as_f64(values):
    if _np is not None and isinstance(values, _np.ndarray):
        arr = _np.ascontiguousarray(values, dtype=_np.float64)
        return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(arr), arr
    buf = (ctypes.c_double * len(values))(*values)
    return buf, len(values), buf


class Cut:
    """A flat partition: the tree cut at one place.

    Attributes (numpy arrays when numpy is available, lists otherwise):

    ``groups``    list of arrays, the object ids of each cluster
    ``cohesion``  mean intra-cluster similarity of each cluster
    ``labels_``   one cluster id per object, length ``num_objects`` (scipy/sklearn
                  convention: the label of cluster ``g`` is ``g``, unclustered objects
                  get their own group, so ``labels_`` always covers every object)
    ``members``   all member ids concatenated; ``offsets`` delimits the groups
    """

    def __init__(self, lib, handle, n_objects):
        self._lib = lib
        self._h = handle
        self.n_objects = n_objects
        self.n_groups = lib.dll.spclust_cut_n_groups(handle)
        self.n_members = lib.dll.spclust_cut_n_members(handle)

        om = lib.dll.spclust_cut_offsets(handle)
        self.offsets = [om[i] for i in range(self.n_groups + 1)]
        mm = lib.dll.spclust_cut_members(handle)
        self.members = [mm[i] for i in range(self.n_members)]
        cm = lib.dll.spclust_cut_cohesion(handle)
        self.cohesion = [cm[i] for i in range(self.n_groups)]

        if _np is not None:
            self.offsets = _np.asarray(self.offsets, dtype=_np.int64)
            self.members = _np.asarray(self.members, dtype=_np.uint32)
            self.cohesion = _np.asarray(self.cohesion, dtype=_np.float64)

        self.groups = [self.members[self.offsets[g]:self.offsets[g + 1]]
                       for g in range(self.n_groups)]

        labels = [0] * n_objects
        for g in range(self.n_groups):
            for obj in self.members[self.offsets[g]:self.offsets[g + 1]]:
                labels[int(obj)] = g
        self.labels_ = _np.asarray(labels, dtype=_np.int32) if _np is not None else labels

        lib.dll.spclust_cut_free(handle)
        self._h = None

    def __len__(self):
        return self.n_groups

    def __repr__(self):
        sizes = [len(g) for g in self.groups]
        shown = sizes[:8]
        more = "" if len(sizes) <= 8 else f", ... ({len(sizes)} groups)"
        return f"<spclust.Cut {self.n_groups} clusters, sizes {shown}{more}>"

    def sizes(self):
        return [len(g) for g in self.groups]


class HAC:
    """One hierarchical clustering.

    Parameters
    ----------
    linkage : str or int
        ``"average"`` (default), ``"max-similarity"`` (single linkage in distance
        terms) or ``"min-nonzero"``; the constants ``AVERAGE`` / ``MAX_SIMILARITY`` /
        ``MIN_NONZERO`` (or the CLI's 0/1/2) also work.
    library : Lib, optional
        An already loaded library, to avoid the lookup (see :func:`load_library`).
    """

    def __init__(self, linkage="average", library=None):
        self._lib = library or get_library()
        code = _resolve_linkage(linkage)
        h = self._lib.dll.spclust_hac_new(code)
        if not h:
            raise SpclustError(self._lib.last_error() or "cannot create HAC")
        self._h = h
        self.linkage = code
        self.linkage_name = LINKAGE_NAMES[code]
        self._built = False
        self._n = 0

    # -- lifecycle --------------------------------------------------------
    def close(self):
        """Release the handle. Idempotent, and called for you on ``del``."""
        if self._h:
            self._lib.dll.spclust_hac_free(self._h)
            self._h = None

    def __del__(self):                                  # pragma: no cover
        try:
            self.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    # -- feeding ----------------------------------------------------------
    def reserve(self, n_objects, n_edges):
        """Optional capacity hint; avoids growth while feeding."""
        self._lib.check(self._lib.dll.spclust_hac_reserve(self._h, n_objects, n_edges), "reserve")

    def add_edge(self, a, b, similarity):
        """Feed one edge."""
        self._lib.check(self._lib.dll.spclust_hac_add_edge(self._h, int(a), int(b),
                                                          float(similarity)), "add_edge")

    def add_edges(self, a, b, similarity):
        """Feed many edges in one call, from array-likes of equal length.

        This is the path to use: one foreign call instead of one per edge. Lists,
        ``array.array``, ``memoryview`` and numpy arrays all work; numpy arrays are
        used in place when they are already the right dtype and contiguous, so
        feeding a few million edges costs no extra copy.
        """
        pa, n, keep_a = _as_u32(a)
        pb, nb, keep_b = _as_u32(b)
        ps, ns, keep_s = _as_f64(similarity)
        if not (n == nb == ns):
            raise ValueError(f"length mismatch: a={n}, b={nb}, similarity={ns}")
        if n == 0:
            return
        status = self._lib.dll.spclust_hac_add_edges(self._h, pa, pb, ps, n)
        if status != 0:
            raise SpclustError("add_edges: " + self._explain_bad_edges(a, b, similarity, n))
        del keep_a, keep_b, keep_s

    def _explain_bad_edges(self, a, b, similarity, n):
        """Locate the offending edge(s) after the library rejected a batch.

        The C API reports the reason but not the index (adding an index to every
        per-edge call would cost more than it is worth), so on the error path the
        arrays are checked here instead - vectorised, and only when something is
        actually wrong.
        """
        if _np is None:
            for i in range(n):
                if a[i] == b[i]:
                    return f"self-loop at index {i} (a == b == {a[i]})"
                v = similarity[i]
                if not (v >= 0.0):
                    return f"invalid similarity at index {i}: {float(v)}"
            return self._lib.last_error() or "unknown error"
        aa = _np.asarray(a)
        bb = _np.asarray(b)
        ss = _np.asarray(similarity, dtype=_np.float64)
        bad = _np.nonzero(aa == bb)[0]
        if bad.size:
            i = int(bad[0])
            return f"self-loop at index {i} (a == b == {int(aa[i])})"
        bad = _np.nonzero(~(ss >= 0.0))[0]          # catches NaN too
        if bad.size:
            i = int(bad[0])
            return f"invalid similarity at index {i}: {float(ss[i])}"
        return self._lib.last_error() or "unknown error"

    # -- computing --------------------------------------------------------
    def build(self, floor=0.0):
        """Compute the merge tree. ``floor <= 0`` computes the whole tree.

        The GIL is released for the duration, so other Python threads keep running.
        """
        self._lib.check(self._lib.dll.spclust_hac_build(self._h, float(floor)), "build")
        self._built = True
        self._n = self._lib.dll.spclust_hac_num_objects(self._h)
        return self

    def reset(self):
        """Drop all state so this object can be reused for another clustering."""
        self._lib.check(self._lib.dll.spclust_hac_reset(self._h), "reset")
        self._built = False
        self._n = 0
        return self

    @property
    def num_objects(self):
        if self._h is None:
            return self._n
        return self._lib.dll.spclust_hac_num_objects(self._h)

    # -- results ----------------------------------------------------------
    def tree(self):
        """The merge tree as an ``(n-1, 4)`` array of ``left right similarity size``.

        With numpy this is uint32/uint32/float64/uint64 columns; without it, a list of
        ``(left, right, similarity, size)`` tuples. Note there is no cohesion column:
        use :meth:`tree_full` if you want it (the CLI's ``--matrix`` prints the first
        four, which is the scipy convention).
        """
        rows = self.tree_full()
        if _np is not None:
            # columns of tree_full are (left, right, similarity, cohesion, size); the
            # scipy/CLI convention is (left, right, height, size), so drop the cohesion
            return rows[:, [0, 1, 2, 4]] if len(rows) else rows.reshape(0, 4)
        return [(r[0], r[1], r[2], r[4]) for r in rows]

    def tree_full(self):
        """Rows of ``(left, right, similarity, cohesion, size)`` (numpy structured or tuples)."""
        n = self._lib.dll.spclust_hac_tree_size(self._h)
        buf = (Merge * n)()
        got = ctypes.c_size_t(n)
        self._lib.check(self._lib.dll.spclust_hac_tree(self._h, buf, n, ctypes.byref(got)), "tree")
        if _np is None:
            return [(b.left, b.right, b.similarity, b.cohesion, b.size) for b in buf]
        arr = _np.empty((n, 5), dtype=_np.float64)
        for i in range(n):
            arr[i] = (buf[i].left, buf[i].right, buf[i].similarity, buf[i].cohesion, buf[i].size)
        return arr

    def cut(self, min_similarity):
        """Cut the tree at ``min_similarity`` -> :class:`Cut`."""
        h = self._lib.dll.spclust_hac_cut(self._h, float(min_similarity))
        if not h:
            raise SpclustError(self._lib.last_error() or "cut failed")
        return Cut(self._lib, h, self.num_objects)

    def cut_n_clusters(self, k):
        """Exactly ``k`` clusters -> :class:`Cut` (needs a deep enough tree)."""
        h = self._lib.dll.spclust_hac_cut_n_clusters(self._h, int(k))
        if not h:
            raise SpclustError(self._lib.last_error() or "cut_n_clusters failed")
        return Cut(self._lib, h, self.num_objects)

    def cut_min_size(self, min_size):
        """Merge until no cluster is smaller than ``min_size`` -> :class:`Cut`."""
        h = self._lib.dll.spclust_hac_cut_min_size(self._h, int(min_size))
        if not h:
            raise SpclustError(self._lib.last_error() or "cut_min_size failed")
        return Cut(self._lib, h, self.num_objects)

    def to_scipy(self, distance="complement", reference=0.0, bridge=False, bridge_distance=None):
        """Export the tree as a ``scipy.cluster.hierarchy`` linkage matrix.

        scipy's third column is a **distance** (closest first) while this library works
        in similarities, so a mapping has to be chosen:

        ``"complement"`` (default) ``C - sim`` - the faithful one, because an absent edge
        (similarity 0) becomes distance ``C``, making average linkage equal to scipy's
        ``method='average'`` on the same dense distance matrix;
        ``"negate"`` ``-sim``; ``"reciprocal"`` ``1/sim``.

        ``bridge=True`` appends the merges scipy needs to see exactly ``n-1`` rows when
        the tree is a forest (MinNonzero can leave components that never merge), at
        NaN -> auto distance (``C`` for complement).
        """
        maps = {"complement": COMPLEMENT, "negate": NEGATE, "reciprocal": RECIPROCAL}
        if distance not in maps:
            raise ValueError(f"unknown distance mapping {distance!r}; use {sorted(maps)}")
        d = float("nan") if bridge_distance is None else float(bridge_distance)
        cnt = ctypes.c_size_t(0)
        comps = ctypes.c_size_t(0)
        status = self._lib.dll.spclust_hac_to_scipy(
            self._h, maps[distance], float(reference), 1 if bridge else 0, d,
            None, 0, ctypes.byref(cnt), ctypes.byref(comps))
        self._lib.check(status, "to_scipy")
        rows = (ScipyRow * cnt.value)()
        if cnt.value:
            status = self._lib.dll.spclust_hac_to_scipy(
                self._h, maps[distance], float(reference), 1 if bridge else 0, d,
                rows, cnt.value, ctypes.byref(cnt), ctypes.byref(comps))
            self._lib.check(status, "to_scipy")
        if _np is None:
            return [(r.left, r.right, r.distance, r.size) for r in rows]
        out = _np.empty((cnt.value, 4), dtype=_np.float64)
        for i in range(cnt.value):
            out[i] = (rows[i].left, rows[i].right, rows[i].distance, rows[i].size)
        return out

    def __repr__(self):
        state = f"{self.num_objects} objects, {self._lib.dll.spclust_hac_tree_size(self._h)} merges" \
            if self._built else "not built"
        return f"<spclust.HAC linkage={self.linkage_name!r} {state}>"

// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>
//
// Implementation of the C interface. Three rules:
//   1. C++ exceptions never cross the C boundary - every entry point try/catches into spclust_status;
//   2. nothing returned to C owns a library object (except the handle and spclust_cut, which have free functions);
//   3. error details live in a **thread-local** string, so reading them needs no locking.

#include "spclust/spclust_c.h"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

#include "spclust/hac.h"
#include "spclust/scipy.h"

namespace {

// Thread-local error message. C receives its c_str(), valid until this thread next calls the library.
std::string& error_slot() {
	static thread_local std::string msg;
	return msg;
}

spclust_status fail(spclust_status code, const char* what) {
	error_slot() = (what ? std::string(what) : std::string());
	return code;
}

void clear_error() { error_slot().clear(); }

spclust_status fail_from_exception() {
	try {
		throw;                       // rethrow to inspect the type
	} catch (const std::bad_alloc&) {
		return fail(SPCLUST_ERR_OUT_OF_MEMORY, "out of memory");
	} catch (const std::invalid_argument& e) {
		return fail(SPCLUST_ERR_INVALID_ARGUMENT, e.what());
	} catch (const std::exception& e) {
		return fail(SPCLUST_ERR_UNKNOWN, e.what());
	} catch (...) {
		return fail(SPCLUST_ERR_UNKNOWN, "unknown C++ exception");
	}
}

bool valid_linkage(int linkage) {
	return linkage >= SPCLUST_LINKAGE_AVERAGE && linkage <= SPCLUST_LINKAGE_MIN_NONZERO;
}

spclust::Linkage to_linkage(int linkage) {
	switch (linkage) {
	case SPCLUST_LINKAGE_MAX_SIMILARITY: return spclust::Linkage::MaxSimilarity;
	case SPCLUST_LINKAGE_MIN_NONZERO:    return spclust::Linkage::MinNonzero;
	default:                             return spclust::Linkage::Average;
	}
}

// The real contents of spclust_cut: vector<vector<uint32_t>> flattened into one dimension,
// so the C side never has to deal with nested containers.
struct CutBox {
	std::vector<uint32_t> members;
	std::vector<std::size_t> offsets;   // n_groups + 1
	std::vector<double> cohesion;
};

} // namespace

// ODR: the C header declares an incomplete type; the definition lives here.
struct spclust_cut {
	CutBox box;
};

extern "C" {

const char* spclust_version(void) {
	// Stringified from the header macros so the two versions cannot drift (the test checks it against the release literal)
#define SPCLUST_STR_(x) #x
#define SPCLUST_STR(x) SPCLUST_STR_(x)
	return SPCLUST_STR(SPCLUST_VERSION_MAJOR) "."
		   SPCLUST_STR(SPCLUST_VERSION_MINOR) "."
		   SPCLUST_STR(SPCLUST_VERSION_PATCH);
#undef SPCLUST_STR
#undef SPCLUST_STR_
}

const char* spclust_last_error(void) {
	const std::string& m = error_slot();
	return m.empty() ? NULL : m.c_str();
}

int spclust_linkage_from_name(const char* name) {
	if (!name)
		return -1;
	const char* table[] = {"average", "max-similarity", "min-nonzero"};
	// The aliases the CLI accepts are recognized too (--linkage max / single / min)
	const char* alias[] = {"max", "single", "min"};
	const int alias_id[] = {SPCLUST_LINKAGE_MAX_SIMILARITY,
							SPCLUST_LINKAGE_MAX_SIMILARITY,
							SPCLUST_LINKAGE_MIN_NONZERO};
	for (int i = 0; i < 3; i++)
		if (0 == std::strcmp(name, table[i]))
			return i;
	for (int i = 0; i < 3; i++)
		if (0 == std::strcmp(name, alias[i]))
			return alias_id[i];
	return -1;
}

spclust_hac* spclust_hac_new(int linkage) {
	if (!valid_linkage(linkage)) {
		fail(SPCLUST_ERR_INVALID_ARGUMENT, "linkage must be 0/1/2");
		return NULL;
	}
	clear_error();
	try {
		return reinterpret_cast<spclust_hac*>(new spclust::HAC(to_linkage(linkage)));
	} catch (...) {
		fail_from_exception();
		return NULL;
	}
}

void spclust_hac_free(spclust_hac* h) {
	if (!h)
		return;
	try {
		delete reinterpret_cast<spclust::HAC*>(h);
	} catch (...) {
		// a destructor should not throw; if it does, it must not escape
	}
}

spclust_status spclust_hac_reserve(spclust_hac* h, size_t n_objects, size_t n_edges) {
	if (!h)
		return fail(SPCLUST_ERR_INVALID_ARGUMENT, "handle is NULL");
	clear_error();
	try {
		reinterpret_cast<spclust::HAC*>(h)->reserve(n_objects, n_edges);
		return SPCLUST_OK;
	} catch (...) {
		return fail_from_exception();
	}
}

spclust_status spclust_hac_add_edge(spclust_hac* h, uint32_t a, uint32_t b,
									double similarity) {
	if (!h)
		return fail(SPCLUST_ERR_INVALID_ARGUMENT, "handle is NULL");
	clear_error();
	try {
		reinterpret_cast<spclust::HAC*>(h)->add_edge(a, b, similarity);
		return SPCLUST_OK;
	} catch (...) {
		return fail_from_exception();
	}
}

spclust_status spclust_hac_add_edges(spclust_hac* h, const uint32_t* a, const uint32_t* b,
									const double* sim, size_t count) {
	if (!h)
		return fail(SPCLUST_ERR_INVALID_ARGUMENT, "handle is NULL");
	if (count && (!a || !b || !sim))
		return fail(SPCLUST_ERR_INVALID_ARGUMENT, "a, b and sim must not be NULL when count > 0");
	clear_error();
	try {
		spclust::HAC* hac = reinterpret_cast<spclust::HAC*>(h);
		//Single linkage only accumulates an edge list, so it can reserve once up front and
		//skip the growth of that vector entirely. The greedy path has no such hint.
		if (spclust::Linkage::MaxSimilarity == hac->linkage())
			hac->reserve(0, count);
		for (size_t i = 0; i < count; i++)
			hac->add_edge(a[i], b[i], sim[i]);
		return SPCLUST_OK;
	} catch (...) {
		return fail_from_exception();
	}
}

spclust_status spclust_hac_build(spclust_hac* h, double floor) {
	if (!h)
		return fail(SPCLUST_ERR_INVALID_ARGUMENT, "handle is NULL");
	clear_error();
	try {
		reinterpret_cast<spclust::HAC*>(h)->build(floor);
		return SPCLUST_OK;
	} catch (...) {
		return fail_from_exception();
	}
}

spclust_status spclust_hac_reset(spclust_hac* h) {
	if (!h)
		return fail(SPCLUST_ERR_INVALID_ARGUMENT, "handle is NULL");
	clear_error();
	try {
		reinterpret_cast<spclust::HAC*>(h)->reset();
		return SPCLUST_OK;
	} catch (...) {
		return fail_from_exception();
	}
}

size_t spclust_hac_num_objects(const spclust_hac* h) {
	if (!h)
		return 0;
	return reinterpret_cast<const spclust::HAC*>(h)->num_objects();
}

size_t spclust_hac_tree_size(const spclust_hac* h) {
	if (!h)
		return 0;
	return reinterpret_cast<const spclust::HAC*>(h)->tree().size();
}

int spclust_hac_linkage(const spclust_hac* h) {
	if (!h)
		return -1;
	return (int)reinterpret_cast<const spclust::HAC*>(h)->linkage();
}

const char* spclust_hac_linkage_name(const spclust_hac* h) {
	if (!h)
		return NULL;
	return reinterpret_cast<const spclust::HAC*>(h)->linkage_name();
}

spclust_status spclust_hac_tree(const spclust_hac* h, spclust_merge* out,
								size_t cap, size_t* count) {
	if (!h || !count)
		return fail(SPCLUST_ERR_INVALID_ARGUMENT, "h or count is NULL");
	clear_error();
	try {
		const std::vector<spclust::Merge>& tr =
			reinterpret_cast<const spclust::HAC*>(h)->tree();
		*count = tr.size();
		if (!out)
			return SPCLUST_OK;
		if (cap < tr.size()) {
			for (size_t i = 0; i < cap; i++) {   // fill as much as fits, to help diagnose
				out[i].left = tr[i].left;
				out[i].right = tr[i].right;
				out[i].similarity = tr[i].similarity;
				out[i].cohesion = tr[i].cohesion;
				out[i].size = tr[i].size;
			}
			return fail(SPCLUST_ERR_BUFFER_TOO_SMALL, "output buffer too small (*count holds the required row count)");
		}
		for (size_t i = 0; i < tr.size(); i++) {
			out[i].left = tr[i].left;
			out[i].right = tr[i].right;
			out[i].similarity = tr[i].similarity;
			out[i].cohesion = tr[i].cohesion;
			out[i].size = tr[i].size;
		}
		return SPCLUST_OK;
	} catch (...) {
		return fail_from_exception();
	}
}

// ------------------------------------------------------------------ cutting

namespace {

spclust_cut* box_cut(const spclust::Cut& c) {
	std::unique_ptr<spclust_cut> out(new spclust_cut());
	out->box.offsets.reserve(c.members.size() + 1);
	out->box.offsets.push_back(0);
	std::size_t total = 0;
	for (std::size_t g = 0; g < c.members.size(); g++)
		total += c.members[g].size();
	out->box.members.reserve(total);
	for (std::size_t g = 0; g < c.members.size(); g++) {
		out->box.members.insert(out->box.members.end(),
								c.members[g].begin(), c.members[g].end());
		out->box.offsets.push_back(out->box.members.size());
	}
	out->box.cohesion = c.cohesion;
	return out.release();
}

} // namespace

spclust_cut* spclust_hac_cut(const spclust_hac* h, double min_similarity) {
	if (!h) {
		fail(SPCLUST_ERR_INVALID_ARGUMENT, "handle is NULL");
		return NULL;
	}
	clear_error();
	try {
		return box_cut(reinterpret_cast<const spclust::HAC*>(h)->cut(min_similarity));
	} catch (...) {
		fail_from_exception();
		return NULL;
	}
}

spclust_cut* spclust_hac_cut_n_clusters(const spclust_hac* h, size_t k) {
	if (!h) {
		fail(SPCLUST_ERR_INVALID_ARGUMENT, "handle is NULL");
		return NULL;
	}
	clear_error();
	try {
		return box_cut(reinterpret_cast<const spclust::HAC*>(h)->cut_n_clusters(k));
	} catch (...) {
		fail_from_exception();
		return NULL;
	}
}

spclust_cut* spclust_hac_cut_min_size(const spclust_hac* h, size_t min_size) {
	if (!h) {
		fail(SPCLUST_ERR_INVALID_ARGUMENT, "handle is NULL");
		return NULL;
	}
	clear_error();
	try {
		return box_cut(reinterpret_cast<const spclust::HAC*>(h)->cut_min_size(min_size));
	} catch (...) {
		fail_from_exception();
		return NULL;
	}
}

void spclust_cut_free(spclust_cut* c) {
	if (!c)
		return;
	try {
		delete c;
	} catch (...) {
	}
}

size_t spclust_cut_n_groups(const spclust_cut* c) {
	return c ? c->box.offsets.empty() ? 0 : c->box.offsets.size() - 1 : 0;
}

size_t spclust_cut_n_members(const spclust_cut* c) {
	return c ? c->box.members.size() : 0;
}

const uint32_t* spclust_cut_members(const spclust_cut* c) {
	return (c && !c->box.members.empty()) ? c->box.members.data() : NULL;
}

const size_t* spclust_cut_offsets(const spclust_cut* c) {
	return (c && !c->box.offsets.empty()) ? c->box.offsets.data() : NULL;
}

const double* spclust_cut_cohesion(const spclust_cut* c) {
	return (c && !c->box.cohesion.empty()) ? c->box.cohesion.data() : NULL;
}

// ------------------------------------------------------------------ scipy export

spclust_status spclust_hac_to_scipy(const spclust_hac* h, int distance_map,
									double reference, int bridge, double bridge_distance,
									spclust_scipy_row* out, size_t cap,
									size_t* count, size_t* n_components) {
	if (!h || !count)
		return fail(SPCLUST_ERR_INVALID_ARGUMENT, "h or count is NULL");
	if (distance_map < SPCLUST_DISTANCE_COMPLEMENT || distance_map > SPCLUST_DISTANCE_RECIPROCAL)
		return fail(SPCLUST_ERR_INVALID_ARGUMENT, "distance_map must be 0/1/2");
	clear_error();
	try {
		const spclust::HAC* hac = reinterpret_cast<const spclust::HAC*>(h);
		spclust::ScipyOptions o;
		o.map = (spclust::DistanceMap)distance_map;
		o.reference = reference;
		o.bridge = (bridge != 0);
		o.bridge_distance = bridge_distance;
		spclust::ScipyExport e = spclust::to_scipy(hac->tree(), hac->num_objects(), o);

		*count = e.rows.size();
		if (n_components)
			*n_components = e.n_components;
		if (!out)
			return SPCLUST_OK;
		if (cap < e.rows.size())
			return fail(SPCLUST_ERR_BUFFER_TOO_SMALL, "output buffer too small (*count holds the required row count)");
		for (size_t i = 0; i < e.rows.size(); i++) {
			out[i].left = e.rows[i].left;
			out[i].right = e.rows[i].right;
			out[i].distance = e.rows[i].distance;
			out[i].size = e.rows[i].size;
		}
		return SPCLUST_OK;
	} catch (...) {
		return fail_from_exception();
	}
}

} // extern "C"

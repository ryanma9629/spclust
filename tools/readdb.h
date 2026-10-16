#ifndef READDB_H
#define READDB_H

// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

#include "cli_config.h"

#include <cstddef>
#include <string>

// Each edge is handed to the callback as soon as it is parsed; the whole input is
// never buffered in memory.
//
// The previous design stored every edge in a std::vector<sp_sim> and let
// init_grp walk it afterwards. That buffer was substantial - measured with 1M
// edges, n = 100k, identical topology:
//     object names of 8 chars (SSO)      peak RSS 202 MB
//     object names of 64 chars (heap)    peak RSS 371 MB
// The 169 MB difference is that buffer, and nothing read it after init_grp.
// Streaming removes it entirely.
// Passing by value is deliberate: the reader moves names in, the callback moves
// them into its own container, and no string is copied anywhere.
typedef void (*edge_handler)(std::string obj1, std::string obj2, double sim);

// Read a sparse similarity edge list (obj1,obj2,sim), calling on_edge once per
// edge.
//   path == NULL or "-" reads from standard input, otherwise from that file.
// Returns the number of edges read; **0 means failure** (cannot open the file, an
// IO error, or no records at all), and the reason has already been written to
// stderr. The original implementation called exit(1) directly, which made it
// unusable as a library; now the caller decides what to do.
std::size_t readdb(const char* path, edge_handler on_edge);

#endif

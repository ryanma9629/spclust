// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

#ifndef SPCLUST_CLI_CONFIG_H
#define SPCLUST_CLI_CONFIG_H

// ---------------------------------------------------------------------------
// Compile-time options belonging to the command-line tool. **None of these are
// part of the library**: the library does not print, does not read files, and
// does not restrict the range of similarities (it accepts any finite
// non-negative double).
// They are all wrapped in #ifndef so the command line can override them:
//     make MIN_SIM=60 DEBUG=2
// ---------------------------------------------------------------------------

// How often to print read progress
#ifndef DEBUG_FREQ_READDB
#define DEBUG_FREQ_READDB 10000
#endif

// Default similarity at which to cut the tree (the default of --cut)
#ifndef MIN_SIM_CLUST
#define MIN_SIM_CLUST 50
#endif

// The CLI's input contract: similarity must be within [0, 65535]. This is
// historical (it used to be an unsigned short) and the library has no such limit;
// it is kept so the reader rejects obviously bogus input.
#ifndef MAX_SIM_READ
#define MAX_SIM_READ 65535
#endif

// Debug level: 0 = silent, 1 = phases and progress, 2 = per-merge detail
#ifndef DEBUG
#define DEBUG 1
#endif

#endif

// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

#include "readdb.h"

#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

using namespace std;

namespace {

//Trim leading/trailing whitespace (plus an optional matching pair of double quotes)
string trim(const string& s) {
	size_t b = s.find_first_not_of(" \t\r\n\"");
	if (string::npos == b)
		return string();
	size_t e = s.find_last_not_of(" \t\r\n\"");
	return s.substr(b, e - b + 1);
}

//Split a line into fields: by comma when the line contains one (CSV), otherwise by whitespace
void split_fields(const string& line, vector<string>& out) {
	out.clear();
	if (string::npos != line.find(',')) {
		size_t pos = 0;
		while (true) {
			size_t comma = line.find(',', pos);
			if (string::npos == comma) {
				out.push_back(trim(line.substr(pos)));
				break;
			}
			out.push_back(trim(line.substr(pos, comma - pos)));
			pos = comma + 1;
		}
	} else {
		istringstream ss(line);
		string tok;
		while (ss >> tok)
			out.push_back(tok);
	}
}

//Strictly parse the similarity: a decimal number within [0, MAX_SIM_READ],
//fractions allowed (e.g. 87.5). The upper bound is the CLI's input contract
//(historical; the library itself has no such limit).
bool parse_sim(const string& s, double& out) {
	if (s.empty())
		return false;
	char* end = NULL;
	errno = 0;
	double v = strtod(s.c_str(), &end);
	if (0 != errno || end == s.c_str() || '\0' != *end)
		return false;
	//Note: writing this as (v < 0 || v > MAX_SIM_READ) would let "nan" pass silently
	if (!(v >= 0.0 && v <= (double)MAX_SIM_READ))
		return false;
	out = v;
	return true;
}

//Recognize the optional header line (obj1,obj2,sim), case-insensitively
bool is_header(const vector<string>& f) {
	string a = f[0], b = f[1], c = f[2];
	for (size_t i = 0; i < a.size(); i++) a[i] = (char)tolower((unsigned char)a[i]);
	for (size_t i = 0; i < b.size(); i++) b[i] = (char)tolower((unsigned char)b[i]);
	for (size_t i = 0; i < c.size(); i++) c[i] = (char)tolower((unsigned char)c[i]);
	return "obj1" == a && "obj2" == b && "sim" == c;
}

} // namespace

//Small batches: fill one batch, then call back.
//Purely streaming (calling back per edge) interleaves structure building into the
//parsing loop, so the parse phase's cache working set is polluted by the growing
//structure - measured on 2M edges, "read + build" went from 0.313s to 0.416s,
//7% slower overall. With batching, parsing stays hot and building stays batched,
//at the cost of one batch of memory (64k edges ~ 3.6 MB).
//Comparison (2M edges, 7 alternating repetitions each, median):
//   read everything first   1.303s   401.3 MB
//   purely streaming        1.393s   219.4 MB
//   streaming + 64k batch   1.301s   223.0 MB
static const size_t BATCH_EDGES = 65536;

struct batch_edge {
	string obj1;
	string obj2;
	double sim;
};

size_t readdb(const char* path, edge_handler on_edge) {
	ifstream file;
	istream* in = &cin;
	string source = "<stdin>";

	if (NULL != path && 0 != strcmp(path, "-")) {
		file.open(path);
		if (!file.is_open()) {
			cerr << "Cannot open input file: " << path << endl;
			return 0;                 // 0 = failure, see readdb.h
		}
		in = &file;
		source = path;
	}

	clock_t s_time = clock();
	if (DEBUG >= 1)
		cout << "Reading similarity data from " << source
		     << " (one pass: parse and build together)..." << endl;

	vector<batch_edge> batch;
	batch.reserve(BATCH_EDGES);
	size_t lineno = 0;
	size_t n_read = 0;          // number of edges read so far
	size_t n_skip = 0;
	bool header_skipped = false;
	bool warned_extra = false;
	string line;
	vector<string> f;

	while (getline(*in, line)) {
		lineno++;
		if (!line.empty() && '\r' == line[line.size() - 1])
			line.erase(line.size() - 1);
		//Skip blank lines and comments
		if (string::npos == line.find_first_not_of(" \t"))
			continue;
		if ('#' == line[line.find_first_not_of(" \t")])
			continue;

		split_fields(line, f);
		if (f.size() < 3) {
			cerr << "Line " << lineno << ": expected 3 fields (obj1,obj2,sim), got "
				<< f.size() << "; skipped." << endl;
			n_skip++;
			continue;
		}
		if (f.size() > 3 && !warned_extra) {
			cerr << "Line " << lineno << ": " << (f.size() - 3)
				<< " extra field(s) ignored." << endl;
			warned_extra = true;
		}

		double sim = 0;
		if (!parse_sim(f[2], sim)) {
			//Only the very first line may be a header; anything else is an error
			if (!header_skipped && 0 == n_read && is_header(f)) {
				header_skipped = true;
				if (DEBUG >= 1)
					cout << "\tSkipped header line." << endl;
				continue;
			}
			cerr << "Line " << lineno << ": bad sim value \"" << f[2] << "\"; skipped." << endl;
			n_skip++;
			continue;
		}
		if (f[0].empty() || f[1].empty()) {
			cerr << "Line " << lineno << ": empty object name; skipped." << endl;
			n_skip++;
			continue;
		}
		//A self-loop (diagonal cell) would break the "group1 < group2" invariant,
		//so it is dropped
		if (f[0] == f[1]) {
			cerr << "Line " << lineno << ": self-pair \"" << f[0] << "\"; skipped." << endl;
			n_skip++;
			continue;
		}

		//Accumulate into a small batch first (names are moved in, not copied)
		batch_edge be;
		be.obj1 = std::move(f[0]);
		be.obj2 = std::move(f[1]);
		be.sim = sim;
		batch.push_back(std::move(be));
		if (batch.size() >= BATCH_EDGES) {
			//Batch is full: hand it to the callback to build the structure
			for (size_t k = 0; k < batch.size(); k++)
				on_edge(std::move(batch[k].obj1), std::move(batch[k].obj2), batch[k].sim);
			batch.clear();
		}
		n_read++;
		if (DEBUG >= 1 && 0 == n_read % DEBUG_FREQ_READDB)
			cout << "\t" << n_read << " records read." << endl;
	}

	//Flush the final partial batch
	for (size_t k = 0; k < batch.size(); k++)
		on_edge(std::move(batch[k].obj1), std::move(batch[k].obj2), batch[k].sim);
	batch.clear();

	if (in->bad()) {
		cerr << "I/O error while reading " << source << "." << endl;
		return 0;
	}
	if (0 == n_read) {
		cerr << "No records to read." << endl;
		return 0;
	}

	clock_t e_time = clock();
	if (DEBUG >= 1) {
		cout << "Done. (" << n_read << " edges";
		if (n_skip > 0)
			cout << ", " << n_skip << " line(s) skipped";
		cout << ", " << (double)(e_time - s_time) / CLOCKS_PER_SEC << "s)" << endl << endl;
	}

	return n_read;
}

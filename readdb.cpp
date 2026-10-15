#include "readdb.h"

#include <cctype>
#include <cerrno>
#include <fstream>
#include <sstream>

//相似性矩阵，稀疏形式
vector<sp_sim> sims;

namespace {

//去掉首尾空白（以及可有可无的成对双引号）
string trim(const string& s) {
	size_t b = s.find_first_not_of(" \t\r\n\"");
	if (string::npos == b)
		return string();
	size_t e = s.find_last_not_of(" \t\r\n\"");
	return s.substr(b, e - b + 1);
}

//把一行拆成字段：含逗号就按逗号拆（CSV），否则按空白拆
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

//严格解析相似度：必须是 [0,65535] 范围内的十进制整数
bool parse_sim(const string& s, unsigned short& out) {
	if (s.empty())
		return false;
	char* end = NULL;
	errno = 0;
	long v = strtol(s.c_str(), &end, 10);
	if (0 != errno || end == s.c_str() || '\0' != *end)
		return false;
	if (v < 0 || v > 65535)
		return false;
	out = (unsigned short)v;
	return true;
}

//识别可选表头行（obj1,obj2,sim），大小写不敏感
bool is_header(const vector<string>& f) {
	string a = f[0], b = f[1], c = f[2];
	for (size_t i = 0; i < a.size(); i++) a[i] = (char)tolower((unsigned char)a[i]);
	for (size_t i = 0; i < b.size(); i++) b[i] = (char)tolower((unsigned char)b[i]);
	for (size_t i = 0; i < c.size(); i++) c[i] = (char)tolower((unsigned char)c[i]);
	return "obj1" == a && "obj2" == b && "sim" == c;
}

} // namespace

size_t readdb(const char* path) {
	ifstream file;
	istream* in = &cin;
	string source = "<stdin>";

	if (NULL != path && 0 != strcmp(path, "-")) {
		file.open(path);
		if (!file.is_open()) {
			cerr << "Cannot open input file: " << path << endl;
			exit(1);
		}
		in = &file;
		source = path;
	}

	clock_t s_time = clock();
	if (DEBUG >= 1)
		cout << "Reading similarity data from " << source << "..." << endl;

	size_t lineno = 0;
	size_t n_skip = 0;
	bool header_skipped = false;
	bool warned_extra = false;
	string line;
	vector<string> f;

	while (getline(*in, line)) {
		lineno++;
		if (!line.empty() && '\r' == line[line.size() - 1])
			line.erase(line.size() - 1);
		//跳过空行和注释行
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

		sp_sim tmp;
		tmp.obj1 = f[0];
		tmp.obj2 = f[1];
		if (!parse_sim(f[2], tmp.sim)) {
			//只允许文件开头的第一行是表头，其余一律报错
			if (!header_skipped && 0 == sims.size() && is_header(f)) {
				header_skipped = true;
				if (DEBUG >= 1)
					cout << "\tSkipped header line." << endl;
				continue;
			}
			cerr << "Line " << lineno << ": bad sim value \"" << f[2] << "\"; skipped." << endl;
			n_skip++;
			continue;
		}
		if (tmp.obj1.empty() || tmp.obj2.empty()) {
			cerr << "Line " << lineno << ": empty object name; skipped." << endl;
			n_skip++;
			continue;
		}
		//自环（对角线单元）会破坏"始终保证 group1 < group2"的前提，直接丢弃
		if (tmp.obj1 == tmp.obj2) {
			cerr << "Line " << lineno << ": self-pair \"" << tmp.obj1 << "\"; skipped." << endl;
			n_skip++;
			continue;
		}

		sims.push_back(tmp);
		if (DEBUG >= 1 && 0 == sims.size() % DEBUG_FREQ_READDB)
			cout << "\t" << sims.size() << " records read." << endl;
	}

	if (in->bad()) {
		cerr << "I/O error while reading " << source << "." << endl;
		exit(1);
	}
	if (sims.empty()) {
		cerr << "No records to read." << endl;
		exit(1);
	}

	clock_t e_time = clock();
	if (DEBUG >= 1) {
		cout << "Done. (" << sims.size() << " edges";
		if (n_skip > 0)
			cout << ", " << n_skip << " line(s) skipped";
		cout << ", " << (e_time - s_time) / CLOCKS_PER_SEC << "s)" << endl << endl;
	}

	return sims.size();
}

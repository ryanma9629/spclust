#include "readdb.h"

#include <cctype>
#include <cerrno>
#include <fstream>
#include <sstream>
#include <utility>

using namespace std;

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

//严格解析相似度：必须是 [0,65535] 范围内的十进制数（允许小数，如 87.5）
bool parse_sim(const string& s, double& out) {
	if (s.empty())
		return false;
	char* end = NULL;
	errno = 0;
	double v = strtod(s.c_str(), &end);
	if (0 != errno || end == s.c_str() || '\0' != *end)
		return false;
	//注意不能写成 (v < 0 || v > 65535)：那样 "nan" 会静默通过
	if (!(v >= 0.0 && v <= 65535.0))
		return false;
	out = v;
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

//小批量缓冲：读满一批再回调。
//纯流式（每条边立刻回调）会把“建结构”插进解析循环，解析阶段的 cache 工作集被结构增长
//污染 —— 实测 2M 条边时“读+建”从 0.313s 涨到 0.416s，整体慢 7%。加一层小批量之后
//解析阶段是热的、建结构阶段也是批量的，内存只多一批（64k 条 ≈ 3.6 MB）。
//对照（2M 条边，各 7 次交替取中位）：
//   先全读再建      1.303s   401.3 MB
//   纯流式          1.393s   219.4 MB
//   流式 + 64k 批   1.301s   223.0 MB
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
			return 0;                 // 0 = 失败，见 readdb.h
		}
		in = &file;
		source = path;
	}

	clock_t s_time = clock();
	if (DEBUG >= 1)
		cout << "Reading similarity data from " << source
		     << " (每条边边读边建结构)..." << endl;

	vector<batch_edge> batch;
	batch.reserve(BATCH_EDGES);
	size_t lineno = 0;
	size_t n_read = 0;          // 已读入的边数
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

		double sim = 0;
		if (!parse_sim(f[2], sim)) {
			//只允许文件开头的第一行是表头，其余一律报错
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
		//自环（对角线单元）会破坏"始终保证 group1 < group2"的前提，直接丢弃
		if (f[0] == f[1]) {
			cerr << "Line " << lineno << ": self-pair \"" << f[0] << "\"; skipped." << endl;
			n_skip++;
			continue;
		}

		//先攒进一小批（字符串 move 进去，不拷贝）
		batch_edge be;
		be.obj1 = std::move(f[0]);
		be.obj2 = std::move(f[1]);
		be.sim = sim;
		batch.push_back(std::move(be));
		if (batch.size() >= BATCH_EDGES) {
			//读满一批，统一交给回调建结构
			for (size_t k = 0; k < batch.size(); k++)
				on_edge(std::move(batch[k].obj1), std::move(batch[k].obj2), batch[k].sim);
			batch.clear();
		}
		n_read++;
		if (DEBUG >= 1 && 0 == n_read % DEBUG_FREQ_READDB)
			cout << "\t" << n_read << " records read." << endl;
	}

	//冲掉最后不足一批的
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

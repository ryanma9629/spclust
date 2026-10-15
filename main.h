#ifndef MAIN_H
#define MAIN_H

// ---------------------------------------------------------------------------
// 平台相关头文件：只有 Windows 上的 ODBC 版本才需要。
// macOS / Linux 走标准库，不需要任何东西。
// ---------------------------------------------------------------------------
#ifdef _WIN32
#include <tchar.h>
#include <windows.h>
#include <sql.h>
#include <sqltypes.h>
#include <sqlext.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------
// 可调参数。全部用 #ifndef 包住，因此可以在命令行覆盖而不改源码：
//     make LINKAGE=1 MIN_SIM=60 DEBUG=2
// 或在命令行直接给编译器加 -DLINKAGE_MTD=1
// ---------------------------------------------------------------------------

//读取数据过程中调试信息显示频度
#ifndef DEBUG_FREQ_READDB
#define DEBUG_FREQ_READDB 10000
#endif

//聚类过程中链接方法
//0:average linkage
//1:complete linkage
//2:single linkage
#ifndef LINKAGE_MTD
#define LINKAGE_MTD 0
#endif

//相似度下限，低于此值将停止聚类
#ifndef MIN_SIM_CLUST
#define MIN_SIM_CLUST 50
#endif

//聚类过程中调试信息出现的频度，实际应用中建议至少设定为100
#ifndef DEBUG_FREQ_CLUST
#define DEBUG_FREQ_CLUST 100
#endif

#ifndef MIN_ACCT_NUM
#define MIN_ACCT_NUM 2
#endif

#ifndef MIN_SIM_WG
#define MIN_SIM_WG 50
#endif

//调试信息级别：0=静默 1=阶段与进度 2=每次合并的细节
#ifndef DEBUG
#define DEBUG 1
#endif

using namespace std;

//保存原始相似性的结构
struct sp_sim {
	string obj1;
	string obj2;
	double sim;
};




//聚类过程中保存相似性的节点
struct sp_unit {
	sp_unit *right;
	sp_unit *down;
	sp_unit *left;
	sp_unit *up;
	size_t group1;
	size_t group2;
	double sim;
};
//“左向量”中的元素
struct p_dim1 {
	sp_unit *head;
	sp_unit *max;
};
//"顶向量“中的元素
struct p_dim2 {
	sp_unit *head;
};
//聚类结果的元素，组内相似性+成员列表
struct grp_s {
	double sim;
	vector<string> objs;
};
#endif

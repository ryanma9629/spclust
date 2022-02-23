#ifndef MAIN_H
#define MAIN_H

#define DEBUG 1

#include <stdio.h>
#include <tchar.h>
#include <windows.h>
#include <sql.h>
#include <sqltypes.h>
#include <sqlext.h>

#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <map>

//读取数据过程中调试信息显示频度
#define DEBUG_FREQ_READDB 10000

//聚类过程中链接方法
//0:average linkage
//1:complete linkage
//2:single linkage
#define LINKAGE_MTD 0
//相似度下限，低于此值将停止聚类
#define MIN_SIM_CLUST 50
//聚类过程中调试信息出现的频度，实际应用中建议至少设定为100
#define DEBUG_FREQ_CLUST 100

#define MIN_ACCT_NUM 2
#define MIN_SIM_WG 50

using namespace std;

//保存原始相似性的结构
struct sp_sim {
	string obj1;
	string obj2;
	unsigned short sim;
};




//聚类过程中保存相似性的节点
struct sp_unit {
	sp_unit *right;
	sp_unit *down;
	sp_unit *left;
	sp_unit *up;
	size_t group1;
	size_t group2;
	unsigned short sim;
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
	unsigned short sim;
	vector<string> objs;
};
#endif
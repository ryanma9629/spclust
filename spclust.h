#ifndef SPCLUST_H
#define SPCLUST_H

#include "main.h"


uint32_t add_grp(size_t, size_t, double);
uint32_t get_max(void);
uint32_t get_unit(size_t, size_t);
void init_add_edge(string, string, double);
size_t merge_grp(size_t, size_t, double);
//聚类。floor 是"树算到多低"的相似度下限：合并到最大相似度 < floor 就停。
//floor <= 0 表示算完整棵树（一直合到没有正相似度的候选）。
//注意 MIN_SIM_CLUST 不再是停止条件，而是 main.cpp 里"切树"的默认位置。
void clust_grp(double floor);

//合并树（在 clust_grp 里填充）。left/right 是节点 id，同时也是当时的组号。
extern vector<merge_rec> merges;
//初始对象数（= 树的叶子数）。调用前是 0。
size_t tree_object_count(void);
//节点 id -> 对象名（id < tree_object_count() 时有效）
const string& object_name(size_t id);
void free_mem(void);

#endif
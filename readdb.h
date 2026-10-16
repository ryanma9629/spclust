#ifndef READDB_H
#define READDB_H

#include "main.h"

// 每条边解析出来就立刻回调，不再把整个输入攒在内存里。
//
// 原来的做法是先把全部边存进 vector<sp_sim> sims，再让 init_grp 遍历它。
// 但那份缓存占得很多 —— 实测（1M 条边、n=10 万，拓扑相同）：
//   对象名 8 字符（走 SSO）  峰值 RSS 202 MB
//   对象名 64 字符（走堆）   峰值 RSS 371 MB
// 差的 169 MB 就是它；而 init_grp 之后再也没被用过。改成流式之后这份内存直接不存在。
// 按值传参是故意的：读取端 move 进来、回调再 move 进自己的容器，全程零字符串拷贝。
typedef void (*edge_handler)(string obj1, string obj2, double sim);

// 读取稀疏相似度边表（obj1,obj2,sim），每条边调用一次 on_edge。
//   path == NULL 或 "-" 时从标准输入读取，否则从该文件读取。
// 返回读入的边数。
// 沿用原实现的约定：遇到无法打开文件、或一条记录都没有时，打印错误并 exit(1)。
size_t readdb(const char* path, edge_handler on_edge);

#endif

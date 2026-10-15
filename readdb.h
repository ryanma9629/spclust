#ifndef READDB_H
#define READDB_H

#include "main.h"

// 读取稀疏相似度边表（obj1,obj2,sim）。
//   path == NULL 或 "-" 时从标准输入读取，否则从该文件读取。
// 返回读入的边数。
// 沿用原实现的约定：遇到无法打开文件、或一条记录都没有时，打印错误并 exit(1)。
size_t readdb(const char* path);

#endif

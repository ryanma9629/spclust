// SPDX-License-Identifier: MIT
// Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

#ifndef SPCLUST_CLI_CONFIG_H
#define SPCLUST_CLI_CONFIG_H

// ---------------------------------------------------------------------------
// 命令行工具自己的编译期选项。**这些不属于库** —— 库不打印、不读文件，
// 也不规定相似度的取值范围（库接受任意有限非负 double）。
// 全部用 #ifndef 包住，可以在命令行覆盖：make MIN_SIM=60 DEBUG=2
// ---------------------------------------------------------------------------

//读取数据过程中调试信息显示频度
#ifndef DEBUG_FREQ_READDB
#define DEBUG_FREQ_READDB 10000
#endif

//切树的默认相似度位置（--cut 的默认值）
#ifndef MIN_SIM_CLUST
#define MIN_SIM_CLUST 50
#endif

//CLI 的输入契约：相似度必须在 [0, 65535] 内。这是历史遗留（原来是 unsigned short），
//库本身没有这个上限；保留它是为了让读入端能挡住明显离谱的输入。
#ifndef MAX_SIM_READ
#define MAX_SIM_READ 65535
#endif

//调试信息级别：0=静默 1=阶段与进度 2=每次合并的细节
#ifndef DEBUG
#define DEBUG 1
#endif

#endif

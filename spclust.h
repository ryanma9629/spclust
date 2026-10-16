#ifndef SPCLUST_H
#define SPCLUST_H

#include "main.h"


uint32_t add_grp(size_t, size_t, double);
uint32_t get_max(void);
uint32_t get_unit(size_t, size_t);
void init_add_edge(string, string, double);
size_t merge_grp(size_t, size_t, double);
void clust_grp(void);
void free_mem(void);

#endif
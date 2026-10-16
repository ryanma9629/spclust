#ifndef SPCLUST_H
#define SPCLUST_H

#include "main.h"


sp_unit* add_grp(size_t, size_t, double);
sp_unit* get_max(void);
sp_unit* get_unit(size_t, size_t);
void init_add_edge(string, string, double);
size_t merge_grp(size_t, size_t, double);
void clust_grp(void);
void free_mem(void);

#endif
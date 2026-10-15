#include "spclust.h"
#include "readdb.h"

//“左向量”，记录每个横向链的头、值最大节点
vector<p_dim1> vec_dim1;
//“顶向量”，记录每个纵向链的头
vector<p_dim2> vec_dim2;
//将对象名称映射到初始组
map<string, size_t> map_obj_grp;
extern vector<sp_sim> sims;
extern vector<grp_s> grp;

//将一个(group1, group2, sim)的节点加入结构，每次加入所在链表的头位置
sp_unit* add_grp(size_t group1, size_t group2, double sim) {
	//始终保证group1小于group2
	if (group1 > group2)
		swap(group1, group2);

	//”左向量“和”顶向量“扩容
	if (vec_dim1.size() <= group1)
		vec_dim1.resize(group1 + 1);
	if (vec_dim2.size() <= group2)
		vec_dim2.resize(group2 + 1);

	//待加入的节点
	sp_unit* unit = (sp_unit*)malloc(sizeof(sp_unit));
	if (NULL == unit) {
		cerr << "Malloc error." << endl;
		exit(1);
	}
	unit->group1 = group1;
	unit->group2 = group2;
	unit->sim = sim;
	unit->right = NULL;
	unit->left = NULL;
	unit->up = NULL;
	unit->down = NULL;

	//如果该行的max不存在或者比待加入的节点小，更新”左向量“的max指针
	if (NULL == vec_dim1[group1].max || vec_dim1[group1].max->sim < sim)
		vec_dim1[group1].max = unit;

	//原来的横向链表头
	sp_unit* prev_head = vec_dim1[group1].head;
	//新加入的节点成为链表头
	vec_dim1[group1].head = unit;
	//如果原链表头存在，右移一个位置，更新原链表头和新节点之间的链接关系
	if (NULL != prev_head) {
		unit->right = prev_head;
		prev_head->left = unit;
	}

	//原来的纵向链表头，余同上
	prev_head = vec_dim2[group2].head;
	vec_dim2[group2].head = unit;
	if (NULL != prev_head) {
		unit->down = prev_head;
		prev_head->up = unit;
	}
	return unit;
}

//获取全局结构中值最大的节点
sp_unit* get_max(void) {
	double max = 0;
	sp_unit* unit = NULL;
	//沿着”左向量“遍历每一行链表，每行最大值节点的指针已经记录在max中
	for (size_t i = 0; i < vec_dim1.size(); i++) {
		//这一行链表可能已经在聚合过程中被删除，需判断
		if (NULL == vec_dim1[i].max)
			continue;
		//不断更新当前找到的最大值
		else if (vec_dim1[i].max->sim > max) {
			unit = vec_dim1[i].max;
			max = unit->sim;
		}
	}
	return unit;
}

//根据(group1, group2, sim)定位该节点
sp_unit* get_unit(size_t group1, size_t group2) {
	//始终保证group1小于group2
	if (group1 > group2)
		swap(group1, group2);
	//如果越界
	if (group1 >= vec_dim1.size() || group2 >= vec_dim2.size())
		return NULL;
	//沿着vec_dim1中group1那一行链表遍历节点
	sp_unit* unit = vec_dim1[group1].head;
	while(NULL != unit) {
		if (group2 == unit->group2)
			break;
		else
			unit = unit->right;
	}
	return unit;
}

size_t init_grp(void) {
	//把每个对象映射为初始组号
	for (size_t i = 0; i < sims.size(); i++) {
		//obj1是否尚未被映射？
		if ( map_obj_grp.end() == map_obj_grp.find(sims[i].obj1)) {
			//映射为当前最大的组号+1
			map_obj_grp.insert(pair<string, size_t>(sims[i].obj1, map_obj_grp.size()));
			//”分组向量“扩容
			grp.resize(map_obj_grp.size());
			//”分组向量“记录obj1自成的新组
			grp[map_obj_grp.size() - 1].sim = 0;
			grp[map_obj_grp.size() - 1].objs.push_back(sims[i].obj1);
		}

		//obj2是否尚未被映射？
		if (map_obj_grp.end() == map_obj_grp.find(sims[i].obj2)) {
			map_obj_grp.insert(pair<string, size_t>(sims[i].obj2, map_obj_grp.size()));
			grp.resize(map_obj_grp.size());
			grp[map_obj_grp.size() - 1].sim = 0;
			grp[map_obj_grp.size() - 1].objs.push_back(sims[i].obj2);
		}
		//把组信息而不是对象信息加入全局结构
		add_grp(map_obj_grp[sims[i].obj1], map_obj_grp[sims[i].obj2], sims[i].sim);
	}
	return grp.size();
}

//将group1和group2两组聚合，形成新组grp_k，更新其它组到grp_k的相似度
size_t merge_grp(size_t group1, size_t group2, double max_sim) {
	//“左向量”的长度，由于把两个组合并过程中该向量长度会改变，用dim1_size记录合并前的长度，不能用vec_dim1.size()
	size_t dim1_size = vec_dim1.size();
	//"顶向量“的长度，由于把两个组合并过程中该向量长度会改变，用dim2_size记录合并前的长度，不能用vec_dim2.size()
	size_t dim2_size = vec_dim2.size();

	//聚合之前group1和group2两组的对象个数
	size_t group1_size = grp[group1].objs.size();
	size_t group2_size = grp[group2].objs.size();
	//”分组向量“的长度
	size_t grp_size = grp.size();

	//标记”左向量“中需要更新的行
	vector<size_t> update_dim1(dim1_size);
	//标记”顶向量“中需要更新的列
	vector<size_t> update_dim2(dim2_size);
	//标记”左向量“中需要更新最大值的行
	vector<size_t> update_max(dim1_size);

	//所有形如(group1, ..., ...)链上节点所在的列都需要更新
	sp_unit* unit = vec_dim1[group1].head;
	while (NULL != unit) {
		if (group2 != unit->group2)
			update_dim2[unit->group2] = 1;
		unit = unit->right;
	}

	//所有形如(group2, ..., ...)链上节点所在的列都需要更新，如果该链存在
	if (vec_dim1.size() > group2) {
		unit = vec_dim1[group2].head;
		while (NULL != unit) {
			update_dim2[unit->group2] = 1;
			unit = unit->right;
		}
	}

	//所有形如(..., group1, ...)链上节点所在的行都需要更新
	unit = vec_dim2[group1].head;
	while (NULL != unit) {
		update_dim1[unit->group1] = 1;
		if (unit == vec_dim1[unit->group1].max)
			update_max[unit->group1] = 1;
		unit = unit->down;
	}

	//所有形如(..., group2, ...)链上节点所在的行都需要更新
	unit = vec_dim2[group2].head;
	while (NULL != unit) {
		if (group1 != unit->group1) {
			update_dim1[unit->group1] = 1;
			if (unit == vec_dim1[unit->group1].max)
				update_max[unit->group1] = 1;
		}
		unit = unit->down;
	}

	//沿着”左向量“遍历每个需要更新的行
	for (size_t i = 0; i < dim1_size; i++) {
		if (update_dim1[i]) {
			//获得(i, group1, sim)这个节点的相似度，节点有可能不存在
			unit = get_unit(i, group1);
			double sim1 = (NULL != unit) ? unit->sim : 0;
			//获得(i, group2, sim)这个节点的相似度，节点有可能不存在
			unit = get_unit(i, group2);
			double sim2 = (NULL != unit) ? unit->sim : 0;

			//Average Linkage
			if (0 == LINKAGE_MTD) {
				double sim = (sim1 * group1_size + sim2 * group2_size) / 
					(group1_size + group2_size);
				if (sim > 0)
					add_grp(i, grp_size, sim);
			//Complete Linkage
			} else if (1 == LINKAGE_MTD)
				add_grp(i, grp_size, max(sim1, sim2));
			//Single Linkage
			else if (2 == LINKAGE_MTD)
				if (min(sim1, sim2) > 0)
					add_grp(i, grp_size, min(sim1, sim2));
		}
	}

	//沿着顶向量遍历每个需要更新的列
	for (size_t i = 0; i < dim2_size; i++) {
		if (update_dim2[i]) {
			unit = get_unit(i, group1);
			double sim1 = (NULL != unit) ? unit->sim : 0;
			unit = get_unit(i, group2);
			double sim2 = (NULL != unit) ? unit->sim : 0;
			if (0 == LINKAGE_MTD) {
				double sim = (sim1 * group1_size + sim2 * group2_size) / 
					(group1_size + group2_size);
				if (sim > 0)
					add_grp(i, grp_size, sim);
			} else if (1 == LINKAGE_MTD)
				add_grp(i, grp_size, max(sim1, sim2));
			else if (2 == LINKAGE_MTD)
				if (min(sim1, sim2) > 0)
					add_grp(i, grp_size, min(sim1, sim2));
		}
	}

	size_t new_grp = grp.size();
	//”分组向量“扩容
	grp.resize(new_grp + 1);
	//扩容后的空间先存入group1内所有对象
	grp[new_grp].objs = grp[group1].objs;
	//再存入group2内所有对象（append）
	grp[new_grp].objs.insert(grp[new_grp].objs.end(), 
		grp[group2].objs.begin(), grp[group2].objs.end());
	if (0 == LINKAGE_MTD)
		grp[new_grp].sim = (grp[group1].sim * group1_size * (group1_size - 1)
			+ grp[group2].sim * group2_size * (group2_size - 1) + group1_size * group2_size * max_sim * 2) 
			/ ((group1_size + group2_size) * (group1_size + group2_size - 1));
	else if (1 == LINKAGE_MTD)
		grp[new_grp].sim = max(max(grp[group1].sim, grp[group2].sim), max_sim);
	else if (2 == LINKAGE_MTD)
		grp[new_grp].sim = min(min(grp[group1].sim, grp[group2].sim), max_sim);

	//清除group1和group2两组内的对象
	grp[group1].objs.clear();
	grp[group1].sim = 0;
	grp[group2].objs.clear();
	grp[group2].sim = 0;

	//定位(group1, group2, ...)这个点
	sp_unit* unit_max = get_unit(group1, group2);
	sp_unit* tmp = NULL;

	//更新链表的链接关系
	//删除(group1, ..., ...)一整行
	unit = vec_dim1[group1].head;
	while (NULL != unit) {
		//先不删除(group1,group2,sim)这个节点
		if (unit != unit_max) {
			//如果该节点上方没有节点（列链表头）
			if (NULL == unit->up) {
				//如果该节点下方也没有节点（列链表尾）
				if (NULL == unit->down)
					//该节点在列链表上是个孤立节点，直接把”顶向量“头记录清空
					vec_dim2[unit->group2].head = NULL;
				else {
					//”顶向量“的头记录更新
					vec_dim2[unit->group2].head = unit->down;
					//该节点的下方节点变成整个列链表的头节点
					unit->down->up = NULL;
				}
			} else {
				//该节点的上方节点，直接向下链接该节点的下方节点
				unit->up->down = unit->down;
				//如果该节点有下方节点
				if (NULL != unit->down)
					//下方节点直接链接上方节点
					unit->down->up = unit->up;
			}
			//把沿着横向的更新路径推进一个节点
			tmp = unit->right;
			//释放当前节点
			free(unit);
			unit = tmp;
		} else
			unit = unit->right;
	}

	//删除(group2, ..., ...)一整行，这行不一定存在需要先判断
	if (dim1_size > group2) {
		unit = vec_dim1[group2].head;
		while (NULL != unit) {
			if (NULL == unit->up) {
				if (NULL == unit->down)
					vec_dim2[unit->group2].head = NULL;
				else {
					vec_dim2[unit->group2].head = unit->down;
					unit->down->up = NULL;
				}
			} else {
				unit->up->down = unit->down;
				if (NULL != unit->down)
					unit->down->up = unit->up;
			}
			tmp = unit->right;
			free(unit);
			unit = tmp;
		}
	}

	//删除(..., group1, ...)一整列
	unit = vec_dim2[group1].head;
	while (NULL != unit) {
		if (NULL == unit->left) {
			if (NULL == unit->right)
				vec_dim1[unit->group1].head = NULL;
			else {
				vec_dim1[unit->group1].head = unit->right;
				unit->right->left = NULL;
			}
		} else {
			unit->left->right = unit->right;
			if (NULL != unit->right)
				unit->right->left = unit->left;
		}
		tmp = unit->down;
		free(unit);
		unit = tmp;
	}

	//删除(..., group2, ...)一整列
	unit = vec_dim2[group2].head;
	while (NULL != unit) {
		if (unit != unit_max) {
			if (NULL == unit->left) {
				if (NULL == unit->right)
					vec_dim1[unit->group1].head = NULL;
				else {
					vec_dim1[unit->group1].head = unit->right;
					unit->right->left = NULL;
				}
			} else {
				unit->left->right = unit->right;
				if (NULL != unit->right)
					unit->right->left = unit->left;
			}
			tmp = unit->down;
			free(unit);
			unit = tmp;
		} else
			unit = unit->down;
	}

	//最后删除(group1,group2,sim)节点
	free(unit_max);
	unit_max = NULL;
	
	//”左向量“中删除下标为group1,group2（如果存在）指针
	vec_dim1[group1].head = NULL;
	vec_dim1[group1].max = NULL;
	if (dim1_size > group2) {
		vec_dim1[group2].head = NULL;
		vec_dim1[group2].max = NULL;
	}
	//”顶向量“中删除下标为group1,group2的指针
	vec_dim2[group1].head = NULL;
	vec_dim2[group2].head = NULL;

	//遍历每一行链表，更新最大相似度的指向
	for (size_t i = 0; i < dim1_size; i++) {
		if (update_max[i]) {
			sp_unit* unit = vec_dim1[i].head;
			double max = 0;
			//从左到右遍历，找到最大值
			if (NULL == unit)
				vec_dim1[i].max = NULL;
			else {
				while (NULL != unit) {
					if (unit->sim > max) {
						max = unit->sim;
						vec_dim1[i].max = unit;
					}
					unit = unit->right;
				}
			}
		}
	}
	return grp.size();
}

//释放结构占用内存
void free_mem(void) {
	sp_unit* tmp = NULL;
	sp_unit* unit = NULL;
	//沿着”左向量“逐行释放链表空间，并清除“左向量”的指针
	for (size_t i = 0; i < vec_dim1.size(); i++) {
		unit = vec_dim1[i].head;
		while (NULL != unit) {
			tmp = unit->right;
			free(unit);
			unit = tmp;
		}
		vec_dim1[i].head = NULL;
		vec_dim1[i].max = NULL;
	}
	//清除“顶向量”的指针
	for (size_t i = 0; i < vec_dim2.size(); i++)
		vec_dim2[i].head = NULL;

	vec_dim1.clear();
	vec_dim2.clear();
}

//聚类主过程
void clust_grp(void) {
	clock_t s_time = clock();
	if (DEBUG >= 1)
		cout << "Initializing..." << endl;
	//将相似度矩阵中的数值填入结构
	init_grp();
	clock_t e_time = clock();
	if (DEBUG >= 1)
		cout << "Done. (" << (e_time - s_time) / CLOCKS_PER_SEC << "s)" << endl << endl;

	s_time = clock();
	if (DEBUG >= 1)
		cout << "Clustering..." << endl;

	//初始组（对象）数
	size_t grp_size = grp.size();
	//如果没有全部合并成1组，或剩余最大相似度小于阈值 MIN_SIM_CLUST
	while (grp_size > 1) {
		sp_unit* unit = get_max();
		double max_sim = (NULL != unit) ? unit->sim : 0;
		if (NULL == unit || max_sim < MIN_SIM_CLUST)
			break;
		else {
			size_t group1 = unit->group1;
			size_t group2 = unit->group2;
			//每次合并，实际组数减少1
			if (DEBUG >= 1 && 0 == grp_size % DEBUG_FREQ_CLUST)
				cout << "\t" << grp_size << " groups left with max_sim=" << max_sim << "." << endl;
			grp_size--;
			size_t new_grp = merge_grp(group1, group2, max_sim) - 1;
			if (DEBUG >= 2)
				cout << "\tGroup " << group1 << " and group " << group2 << 
					" are merged into group " << new_grp << " (sim=" << setprecision(12)
					<< max_sim << ")." << endl;
		}
	}
	//聚类结束，释放结构占用内存
	free_mem();

	e_time = clock();
	if (DEBUG >= 1)
		cout << "Done. (" << (e_time - s_time) / CLOCKS_PER_SEC << "s)" << endl << endl;
}
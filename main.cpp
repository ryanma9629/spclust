#include "main.h"
#include "spclust.h"
#include "readdb.h"

//聚类结果
vector<grp_s> grp;

int _tmain(int argc, _TCHAR* argv[])
{

	readdb();
	//获得聚类结果

	clust_grp();
	
	//打印聚类结果
/*	cout << "Final groups:" << endl;
	for (size_t i = 0; i < grp.size(); i++) {
		if (grp[i].objs.size() >= 2) {
			cout << "\tGroup " << i << " (sim=" << grp[i].sim << "): ";
			for (size_t j = 0; j < grp[i].objs.size(); j++)
				cout << grp[i].objs[j] << " ";
			cout << endl;
		}
	}	*/	

	return 0;
}
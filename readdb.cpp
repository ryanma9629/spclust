#include "readdb.h"

const string connStr = "dsn=mysql";
const string sqlQuery1 = "select count(*) from test.sim100";
const string sqlQuery2 = "select obj1, obj2, sim from test.sim100";

//相似性矩阵，稀疏形式
vector<sp_sim> sims;

void odbc_show_err(unsigned int handletype, const SQLHANDLE& handle)
{
	SQLCHAR sqlStat[1024];
	SQLCHAR message[1024];
	if (SQL_SUCCESS == SQLGetDiagRec(handletype, handle, 1, sqlStat, NULL, message, 1024, NULL))
		cerr << "Message: " << message << endl << "SQLSTATE: " << sqlStat << endl;
}

size_t readdb(void)
{
	SQLHANDLE sqlEnvHandle;    
	SQLHANDLE sqlConnHandle;
	SQLHANDLE sqlStatHandle;

	sp_sim tmp;
	unsigned int cnt;
 
	if (SQL_SUCCESS != SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &sqlEnvHandle))
		exit(1);
 	if (SQL_SUCCESS != SQLSetEnvAttr(sqlEnvHandle,SQL_ATTR_ODBC_VERSION, (SQLPOINTER)SQL_OV_ODBC3, 0)) 
		exit(1);
 	if(SQL_SUCCESS!=SQLAllocHandle(SQL_HANDLE_DBC, sqlEnvHandle, &sqlConnHandle))
		exit(1);
 
	SQLCHAR retConnStr[1024];

	switch (SQLDriverConnect (sqlConnHandle, NULL, (SQLCHAR* )connStr.c_str(), 
			SQL_NTS, retConnStr, 1024, NULL,SQL_DRIVER_NOPROMPT)) {
		case SQL_SUCCESS_WITH_INFO:
			odbc_show_err(SQL_HANDLE_DBC, sqlConnHandle);
			break;
		case SQL_INVALID_HANDLE:
		case SQL_ERROR:
			odbc_show_err(SQL_HANDLE_DBC, sqlConnHandle);
			exit(1);
		default:
			break;
	}
	if (SQL_SUCCESS != SQLAllocHandle(SQL_HANDLE_STMT, sqlConnHandle, &sqlStatHandle))
		exit(1);
	if (SQL_SUCCESS != SQLExecDirect(sqlStatHandle, 
		(SQLCHAR* )sqlQuery1.c_str(), SQL_NTS)) {
		odbc_show_err(SQL_HANDLE_STMT, sqlStatHandle);
		exit(1);
	} else {
		if (SQL_SUCCESS == SQLFetch(sqlStatHandle))
			SQLGetData(sqlStatHandle, 1, SQL_C_ULONG, &cnt, 0, NULL);
	}
	if (cnt > 0) {
		if (SQL_SUCCESS != SQLAllocHandle(SQL_HANDLE_STMT, sqlConnHandle, &sqlStatHandle))
			exit(1);
		if (SQL_SUCCESS != SQLExecDirect(sqlStatHandle, 
			(SQLCHAR* )sqlQuery2.c_str(), SQL_NTS)) {
			odbc_show_err(SQL_HANDLE_STMT, sqlStatHandle);
			exit(1);
		} 
		clock_t s_time = clock();
		if (DEBUG >= 1)
			cout << "Reading data from DB...(" << cnt << " records)" << endl;
		char obj1[32];
		char obj2[32];
		unsigned short sim;
		unsigned int i = 0;
		while (SQL_SUCCESS == SQLFetch(sqlStatHandle)) {
			SQLGetData(sqlStatHandle, 1, SQL_C_CHAR, obj1, 32, NULL);
			SQLGetData(sqlStatHandle, 2, SQL_C_CHAR, obj2, 32, NULL);
			SQLGetData(sqlStatHandle, 3, SQL_C_USHORT, &sim, 0, NULL);
			tmp.obj1 = obj1;
			tmp.obj2 = obj2;
			tmp.sim = sim;
			sims.push_back(tmp);
			if (DEBUG >= 1 && i > 0 && 0 == i % DEBUG_FREQ_READDB)
				cout << "\t" << i << " records read." << endl;
			i++;
		}
		clock_t e_time = clock();
		if (DEBUG >= 1)
			cout << "Done. (" << (e_time - s_time) / CLOCKS_PER_SEC << "s)" << endl << endl;
	} else {
		cerr << "No recoreds to read." << endl;
		exit(1);
	}

	SQLFreeHandle(SQL_HANDLE_STMT, sqlStatHandle );
	SQLDisconnect(sqlConnHandle);
	SQLFreeHandle(SQL_HANDLE_DBC, sqlConnHandle);
	SQLFreeHandle(SQL_HANDLE_ENV, sqlEnvHandle);

	return sims.size();
}
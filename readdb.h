#ifndef READDB_H
#define READDB_H

#include "main.h"

void odbc_show_err(unsigned int, const SQLHANDLE&);
size_t readdb(void);

#endif
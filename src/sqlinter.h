#ifndef SQLINTER_H
#define SQLINTER_H

#include <sqlite3.h>
#include <string>
#include <exception>
#include <iostream>
#include <map>
#include "grid.h"
using namespace std;

enum grid_tables { 
  tab_uknown, tab_general,
  tab_bus, tab_branch, tab_gen, tab_gencost,
  tab_bus_rts, tab_branch_rts, tab_gen_rts, tab_gencost_rts
};
static map<string, grid_tables> map_gridTables;


class sqlInter
{
 public:
  sqlInter();
  ~sqlInter() { sqlite3_close(db); }
  int openDb(string name);
  void printDb(string table);
  int loadDb(string table, grid & gr);
  int load(grid & gr);
  int loadRTS(grid & gr);

  string getStr(grid_tables gt);
  void parseStmt(grid_tables gt, sqlite3_stmt *stmt, grid & gr, bool load);

 private:
  const char* _name;
  sqlite3 *db;
  
};

#endif

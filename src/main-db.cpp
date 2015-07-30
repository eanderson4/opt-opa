#include <stdlib.h>

#include <ctime>
#include "sqlinter.h"
#include "jsoninter.h"
//#include "grid_rts.h"
#include "isjn.h"
#include "iopa.h"


using namespace std;

int main(int argc, char* argv[]){
  
  if(argc<2){
    cout<<"cmd: DB case/30.db \n"
	<<"\trun main-db for case30\n";
    return 1;
  }
  
  sqlInter db;
  grid * gr;
  
  string db_name;
  string str_rts ("rts");
  db_name = argv[1];
  if(db.openDb(db_name)){
    
    bool loaded=false;
    if (db_name.find(str_rts) != string::npos){
      //RTS GRID
      //      gr = new grid_rts;
      //     loaded=db.loadRTS(*gr);
    }
    else{
      gr = new grid;
      loaded=db.load(*gr);
    }
    
    if(loaded){
      gr->buildMap();
      gr->printNums(cout);
      
      cout<<*gr<<endl;
    }
  }

  return 1;
}

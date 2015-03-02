#include "sqlinter.h"

sqlInter::sqlInter() {
  map_gridTables["bus"]=tab_bus;
  map_gridTables["branch"]=tab_branch;
  map_gridTables["gen"]=tab_gen;
  map_gridTables["gencost"]=tab_gencost;
  map_gridTables["general"]=tab_general;

  map_gridTables["bus_rts"]=tab_bus_rts;
  map_gridTables["branch_rts"]=tab_branch_rts;
  map_gridTables["gen_rts"]=tab_gen_rts;
  map_gridTables["gencost_rts"]=tab_gencost_rts;


  cout<<"map_gridTables contains "
      << map_gridTables.size()
      << " entries." <<endl;
}

int sqlInter::openDb(string name) {
  _name=name.c_str();
  cout<<"Attepting to open "<<_name<<endl;
  int rc =sqlite3_open_v2( _name, &db, SQLITE_OPEN_READWRITE, NULL );
  if( rc != SQLITE_OK ){
    cout<< "sqlite3 error: "<<sqlite3_errmsg(db)<<endl;
    sqlite3_close(db);
    return 0;
  }
  else{
    cout<<"DB opened\n";
    return 1;
  }
}

void sqlInter::printDb(string table) {
  grid temp_gr;
  sqlite3_stmt *stmt;
  grid_tables gt = map_gridTables[table];
  string str = getStr(gt);
  cerr<<"SQL: "<<str<<endl;
  int result = sqlite3_prepare(db, str.c_str(), str.length()+1, &stmt, NULL);
  if (result != SQLITE_OK) {
    cerr<<"Failed to prepare database: "<<sqlite3_errstr(result)<<endl;
    return;
  }  
  do {
    result = sqlite3_step (stmt) ;
    if (result == SQLITE_ROW) { //read data
      parseStmt(gt,stmt,temp_gr, false);
    }
  } while (result == SQLITE_ROW) ;

}

int sqlInter::load(grid & gr) {
  
  cout<<"Load Regular"<<endl;
  if(!loadDb("bus",gr)) return 0;
  if(!loadDb("branch",gr)) return 0;
  if(!loadDb("gen",gr)) return 0;
  if(!loadDb("gencost",gr)) return 0;
  
  return 1;
}

int sqlInter::loadRTS(grid & gr) {

  cout<<"Load RTS"<<endl;
  if(!loadDb("bus_rts",gr)) return 0;
  if(!loadDb("branch_rts",gr)) return 0;
  if(!loadDb("gencost_rts",gr)) return 0;
  if(!loadDb("gen_rts",gr)) return 0;


  
  return 1;
}


int sqlInter::loadDb(string table,grid & gr) {
  sqlite3_stmt *stmt;
  grid_tables gt = map_gridTables[table];
  string str = getStr(gt);
  cout<<"SQL: "<<str<<endl;
  int result = sqlite3_prepare(db, str.c_str(), str.length()+1, &stmt, NULL);
  if (result != SQLITE_OK) {
    cerr<<"Failed to prepare database: "<<sqlite3_errstr(result)<<endl;
    return 0;
  }  
  while (sqlite3_step(stmt) == SQLITE_ROW){
      parseStmt(gt,stmt, gr, true);
  }
  return 1;
}




string sqlInter::getStr(grid_tables gt) {
  switch (gt){
  case tab_bus:
    return "select * from bus";
  case tab_branch:
    return "select * from branch";
  case tab_gen:
    return "select * from gen";
  case tab_gencost:
    return "select * from gencost";
  case tab_general:
    return "select * from general";

  case tab_bus_rts:
    return "select * from bus";
  case tab_branch_rts:
    return "select * from branch";
  case tab_gen_rts:
    return "select * from genparm";
  case tab_gencost_rts:
    return "select * from gencost";


  default:
    cerr<< "Table not known" << gt <<endl;
  }
  return ".tables";
}

void sqlInter::parseStmt(grid_tables gt, sqlite3_stmt *stmt, grid & gr, bool load){
  switch (gt) {
  case tab_bus:
    {
    int num=sqlite3_column_int(stmt,0);
    int type=sqlite3_column_int(stmt,1);
    double pd=sqlite3_column_double(stmt,2);
    double qd=sqlite3_column_double(stmt,3);
    double gs=sqlite3_column_double(stmt,4);
    double bs=sqlite3_column_int(stmt,5);
    int area=sqlite3_column_int(stmt,6);
    double vm=sqlite3_column_int(stmt,7);
    double va=sqlite3_column_int(stmt,8);
    double basekv=sqlite3_column_int(stmt,9);
    int zone=sqlite3_column_int(stmt,10);
    double vmax=sqlite3_column_int(stmt,11);
    double vmin=sqlite3_column_int(stmt,12);
    if(load)
      gr.addBus(bus(num,type,pd,qd,gs,bs,area,vm,va,basekv,zone,vmax,vmin));
    else
      cout<<num<<" "<<pd<<" "<<qd<<endl;
    break;
    }
  case tab_branch:
    {
    int num=sqlite3_column_int(stmt,0);
    int fbus=sqlite3_column_int(stmt,1);
    int tbus=sqlite3_column_int(stmt,2);
    double br_r=sqlite3_column_double(stmt,3);
    double br_x=sqlite3_column_double(stmt,4);
    double br_b=sqlite3_column_double(stmt,5);
    double rate_a=sqlite3_column_double(stmt,6);
    double rate_b=sqlite3_column_double(stmt,7);
    double rate_c=sqlite3_column_double(stmt,8);
    double tap=sqlite3_column_double(stmt,9);
    double shift=sqlite3_column_double(stmt,10);
    int br_status=sqlite3_column_int(stmt,11);
    double ang_min=sqlite3_column_double(stmt,12);
    double ang_max=sqlite3_column_double(stmt,13);
    if(load)
      gr.addBranch(branch(num,fbus,tbus,br_r,br_x,br_b,rate_a,rate_b,rate_c,tap,shift,br_status,ang_min,ang_max));
    else
      cout<<num<<" "<<fbus<<" "<<tbus<<" "<<br_x<<" "<<rate_a<<endl;
    break;
    }
  case tab_gen:
    {
    int num=sqlite3_column_int(stmt,0);
    int bus=sqlite3_column_int(stmt,1);
    double pg = sqlite3_column_double(stmt,2);
    double qg = sqlite3_column_double(stmt,3);
    double qmax=sqlite3_column_double(stmt,4);
    double qmin=sqlite3_column_double(stmt,5);
    double vg=sqlite3_column_double(stmt,6);
    double mbase=sqlite3_column_double(stmt,7);
    double status=sqlite3_column_double(stmt,8);
    double pmax=sqlite3_column_double(stmt,9);
    double pmin=sqlite3_column_double(stmt,10);
    if(load)
      gr.addGen(gen(num,bus,pg,qg,qmax,qmin,vg,mbase,status,pmax,pmin));
    else
      cout<<num<<" "<<bus<<" "<<pg<<" "<<pmax<<" "<<pmin<<endl;
    break;
    }
  case tab_gencost:
    {
    int num=sqlite3_column_int(stmt,0);
    int model=sqlite3_column_int(stmt,1);
    double startup=sqlite3_column_double(stmt,2);
    double shutdown=sqlite3_column_double(stmt,3);
    int ncost=sqlite3_column_int(stmt,4);
    double c2=sqlite3_column_double(stmt,5);
    double c1=sqlite3_column_double(stmt,6);
    double c0=sqlite3_column_double(stmt,7);
    if(load)
      gr.addGenCost(num,model,startup,shutdown,ncost,c2,c1,c0);
    else
      cout<<num<<" "<<model<<" "<<startup<<" "<<shutdown<<" "<<ncost<<" "<<c2<<" "<<c1<<" "<<c0<<endl;
    break;
    }


  case tab_bus_rts:
    {
    int num=sqlite3_column_int(stmt,0);
    string name((char *)sqlite3_column_text(stmt,1));
    int type=sqlite3_column_int(stmt,2);
    double pd=sqlite3_column_double(stmt,3);
    double qd=sqlite3_column_double(stmt,4);
    double gs=sqlite3_column_double(stmt,5);
    double bs=sqlite3_column_int(stmt,6);
    int area=sqlite3_column_int(stmt,7);
    double basekv=sqlite3_column_int(stmt,8);
    int zone=sqlite3_column_int(stmt,9);
    if(load)
      gr.addBus(bus(num,type,pd,qd,gs,bs,area,1,1,basekv,zone,1,1));
    else
      cout<<num<<" "<<pd<<" "<<qd<<endl;
    break;
    }
  case tab_branch_rts:
    {
    int num=sqlite3_column_int(stmt,0);
    string name((char *)sqlite3_column_text(stmt,1));
    int fbus=sqlite3_column_int(stmt,2);
    int tbus=sqlite3_column_int(stmt,3);
    double br_length=sqlite3_column_double(stmt,4);
    double perm_lam=sqlite3_column_double(stmt,5);
    double duration=sqlite3_column_double(stmt,6);
    double trans_lam=sqlite3_column_double(stmt,7);
    double br_r=sqlite3_column_double(stmt,8);
    double br_x=sqlite3_column_double(stmt,9);
    double br_b=sqlite3_column_double(stmt,10);
    double rate_a=sqlite3_column_double(stmt,11);
    double rate_b=sqlite3_column_double(stmt,12);
    double rate_c=sqlite3_column_double(stmt,13);
    double tap=sqlite3_column_double(stmt,14);
    if(load)
      gr.addBranch(branch(num,fbus,tbus,br_r,br_x,br_b,rate_a,rate_b,rate_c,tap,0,1,-360,360));
    else
      cout<<num<<" "<<fbus<<" "<<tbus<<" "<<br_x<<" "<<rate_a<<endl;
    break;
    }
  case tab_gen_rts:
    {
    int num=sqlite3_column_int(stmt,0);
    int bus=sqlite3_column_int(stmt,1);
    string gid((char *)sqlite3_column_text(stmt,2));
    int id = sqlite3_column_double(stmt,3);
    double pg = sqlite3_column_double(stmt,4);
    double qg=sqlite3_column_double(stmt,5);
    double qmax=sqlite3_column_double(stmt,6);
    double qmin=sqlite3_column_double(stmt,7);
    double vs=sqlite3_column_double(stmt,8);
    double c2=gr.getC2(gid);
    double c1=gr.getC1(gid);
    string sizess=gid;
    double pmax=atof(sizess.erase(0,1).c_str());
    if(load)
      gr.addGen(gen(num,bus,pg,qg,qmax,qmin,vs,1,1,pmax,0,c2,c1));
    else
      cout<<num<<" "<<bus<<" "<<pg<<" "<<pmax<<endl;
    break;
    }
  case tab_gencost_rts:
    {
    int num=sqlite3_column_int(stmt,0);
    string gid((char *)sqlite3_column_text(stmt,1));
    double size=sqlite3_column_double(stmt,2);
    string type((char *)sqlite3_column_text(stmt,3));
    double ramp=sqlite3_column_double(stmt,4);
    double c2=sqlite3_column_double(stmt,5);
    double c1=sqlite3_column_double(stmt,6);
    if(load)
      gr.addGenType(gentype(num,gid,size,type,ramp,c2,c1));
    else
      cout<<num<<" "<<size<<" "<<c2<<" "<<c1<<endl;
    break;
    }

  case tab_general:
    {
    int num=sqlite3_column_int(stmt,0);
    string par((char *)sqlite3_column_text(stmt,1));
    string val((char *)sqlite3_column_text(stmt,2));
    cerr<<num<<" "<<par<<" "<<val<<endl;
    break;
    }
  default:
    cerr << "Table not known" <<endl;
  }
}


  

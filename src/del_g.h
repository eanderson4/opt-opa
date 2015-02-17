#ifndef DELG_H
#define DELG_H

#include <vector>
#include "grid.h"
using namespace std;


class del_g
{
 public:
  del_g(): have_demand(false), have_topo(false), have_cap(false){ }
 del_g(grid * gr): _gr(gr), have_demand(true), have_topo(false), have_cap(false){ del_demand.resize(gr->numBuses(),0); baseTopo(gr); }
  ~del_g(){}
  
  void addDemand(int num,double demand);
  void setStatus(int num,bool status);
  void baseTopo(grid * gr);

  bool haveTopo(){ return have_topo; }
  bool haveDemand(){ return have_demand; }
  double getStatus(int num){ return del_topo[num]; }
  double getDemand(int num){ return del_demand[num]; }

  friend ostream& operator<<(ostream& os, const del_g& dg);

 private:
  
  grid * _gr;

  bool have_demand;
  bool have_topo;
  bool have_cap;

  vector< double > del_demand;
  vector< bool > del_topo;
  vector< double > del_cap;

};

#endif

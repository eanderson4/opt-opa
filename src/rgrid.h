#ifndef RGRID_H
#define RGRID_H

#include <ilcplex/ilocplex.h>
#include "rbase.h"
#include "grid.h"

class rgrid : public rbase {

 public:
  
 rgrid() : rbase (), _loadshed(0),have_gencost(false) { }
  ~rgrid(){  }
  
  double getObjective(){ return objective; }
  IloCplex::CplexStatus getStatus(){ return status; }

  void getSolveInfo(IloCplex * cplex, double rt);
  void setGenCost(double genCost);
  void outputInfo(ostream & out);

  void displayOperatingPos(grid * gr);
  void displayBranchPos(grid * gr);
  void displayGenPos(grid * gr);
  void setLoadShed(double ls){ _loadshed=ls; }


 private:

  double realTime;
  double time;
  double objective;
  double gap;
  double _budget;
  double _genCost;
  double _loadshed;

  int iterations;
  int rootalg;
  
  string rootname;

  bool isQO;
  bool have_gencost;

  IloCplex::CplexStatus status;
  IloCplex::Algorithm alg;
  IloCplex::Algorithm subalg;
  IloCplex::Algorithm test;


};

#endif

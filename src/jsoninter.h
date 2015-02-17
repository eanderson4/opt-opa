#ifndef JSONINTER_H
#define JSONINTER_H


#include "grid.h"
#include "rgrid.h"
#include "gridcalc.h"
using namespace std;


class jsonInter
{
 public:
  jsonInter(){}
  ~jsonInter() { }

  void print(grid * gr, rgrid * rg);
  void printDetail(ostream & out,grid * gr, rgrid * rg,double L, double parm, double SD, double o,double r,vec f, vec g, vec beta, vec sd, vec z, vec p);

 private:
  
};

#endif

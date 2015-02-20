#ifndef IOPA_H
#define IOPA_H


#include "gridcalc.h"


using namespace std;

class iopa {

 public:
 iopa( grid * gr, gridcalc * gc, double L, double p ) : _gr(gr), _gc(gc),_L(L),_p(p) { }

  void runTrials(ostream & out,vec z,int N,double num,running_stat<double> st_tot);

 private:
  
  grid * _gr;
  gridcalc * _gc;


  double _L;
  double _p;


};


#endif


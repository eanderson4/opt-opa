#ifndef IOPA_H
#define IOPA_H


#include "gridcalc.h"


using namespace std;

class iopa {

 public:
 iopa( grid * gr, gridcalc * gc,  vec z, double L, double p ) : _gr(gr), _gc(gc),_z(z),_L(L),_p(p) { }

  void runTrials(double num);

 private:
  
  grid * _gr;
  gridcalc * _gc;


  vec _z;
  double _L;
  double _p;


};


#endif


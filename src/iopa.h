#ifndef IOPA_H
#define IOPA_H


#include "ijn1.h"
#include "gridcalc.h"


using namespace std;

class iopa {

 public:
 iopa( grid * gr, gridcalc * gc, double L, double p, double Lr, double pr ) : _gr(gr), _gc(gc),_L(L),_p(p),_Lr(Lr),_pr(pr) { }

  void runTrials(ostream & out,vec z,int N,double num);
  void runTrials(ostream & out,ijn1 * n1,vec f, vec g, mat SIG,double num, double cost);

 private:
  
  grid * _gr;
  gridcalc * _gc;


  double _L;
  double _p;
  double _Lr;
  double _pr;


};


#endif


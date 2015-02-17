#ifndef IN1_H
#define IN1_H

#include "gridcalc.h"
#include "igrid.h"

using namespace std;


class in1 : public igrid {

 public:
 in1(grid * gr,mat SIGy, mat Hw, double m1): igrid(gr), _SIGy(SIGy), _Hw(Hw), _m1( m1 ) { setupOld(); }
 in1(grid * gr,gridcalc * gc,mat SIGy, mat Hw, double m1): igrid(gr), _gc(gc), _SIGy(SIGy), _Hw(Hw), _m1( m1 ) { setup(); }
  void setup();
  void setupOld();
  rgrid * solveModel( isolve * is=NULL);
  bool postN1(int n,vec yn, IloCplex * cplex);
  vec getN1(int n, vec y0, vec g);

  class itlimit: public exception 
  {
    virtual const char* what() const throw()
    {
      return "Stopped due to iteration limit";
    }
  }itlimit;


 private:
  gridcalc * _gc;

  vector<IloNumVarArray> _yplus;
  vector<IloRangeArray> _yup;
  vector<IloRangeArray> _ydown;

  int addedBounds;

  mat _C;
  mat _Cg;
  mat _Hb;
  mat _L;
  vec _var0;
  mat _var;
  vec _check;
  mat _in;
  mat _SIGy;
  mat _Hw;  
  double _m1;




};
#endif

#ifndef ISJN_H
#define ISJN_H

#include "isj.h"
#include "igrid.h"

using namespace std;


class isjn : public isj {

 public:
 isjn(grid * gr,gridcalc * gc, mat SIGm,vec indexM, double L, double p, double pc, double eps, vec epsN, double epsG): isj( gr, gc, SIGm, indexM, L, p, pc, eps, epsG ), _epsN(epsN) { setup(); }
  void setup();
  rgrid * solveModel( isolve * is=NULL);
  bool postN1(int n,vec yn, vec zn,vec beta, vec sdn, IloCplex * cplex, int iteration=0);
  vec getN1(int n, vec y0);
  vec getSDN(int n, vec y0, mat Cov);


  IloNumVarArray getZ(int n){ return _z[n]; }
  IloNumVarArray getYplus(int n){ return _yplus[n]; }
  IloNumVarArray getSD(int n){ return _sd[n]; }

  vec getEpsN(){ return _epsN; }
  mat getLmat(){ return _L; }
  mat getHb(){ return _Hb; }
  mat getSigPsi(){ return _sigpsi; }

  vec getCheck(){ return _check; }
  mat getAddCut(){ return _addCut; }
  mat getIn(){ return _in; }


 private:
  IloNumVarArray _risk;
  IloRangeArray _riskConstraint;
  vector<IloNumVarArray> _z;
  vector<IloNumVarArray> _yplus;
  vector<IloNumVarArray> _sd;
  vector<IloRangeArray> _yup;
  vector<IloRangeArray> _ydown;
  
  vec _epsN;

  mat _L;
  mat _Hb;

  vec _check;
  mat _addCut;
  mat _in;

  mat _sigpsi;

};
#endif

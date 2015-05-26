#ifndef IOJ_H
#define IOJ_H

#include "isjn.h"
#include "igrid.h"

using namespace std;


class ioj : public isjn {

 public:
 ioj(grid * gr,gridcalc * gc, mat SIGm,vec indexM, double L, double p, double pc, double eps, vec epsN, double epsG, double epsLS, vec epsLSN): isjn(grid * gr,gridcalc * gc, mat SIGm,vec indexM, double L, double p, double pc, double eps, vec epsN, double epsG), _epsLS(epsLS),_epsLSN(epsLSN) { setup(); }
  void setup();
  rgrid * solveModel( isolve * is=NULL);
  bool postN1(int n,vec yn, vec zn,vec beta, vec sdn, IloCplex * cplex, int iteration=0);

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

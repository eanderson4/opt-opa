#ifndef IOJ_H
#define IOJ_H

#include "isjn.h"
#include "igrid.h"

using namespace std;


class ioj : public isjn {

 public:
 ioj(grid * gr,gridcalc * gc, mat SIGm,vec indexM, double L, double p, double pc, double eps, vec epsN, double epsG, double epsLS, vec epsLSN, vec xdesign): isjn(gr,gc,SIGm,indexM,L,p,pc,eps, epsN,epsG), _epsLS(epsLS),_epsLSN(epsLSN), _xdes(xdesign) { setup(); }
  void setup();
  rgrid * solveModel( isolve * is=NULL);

  bool postLS(vec f, vec z, vec l,vec beta, vec SIGy,IloCplex * cplex, int iteration=0);
  bool postLSN1(int n,vec yn, vec zn,vec ln, vec beta, vec sdn, IloCplex * cplex, int iteration=0);


 private:
  IloNumVar _riskLS;
  IloRange _riskConstraintLS;
  
  IloNumVarArray _l0;
  IloRangeArray _l0eq;

  IloNumVarArray _riskLSN;
  IloRangeArray _riskConstraintLSN;

  vector<IloNumVarArray> _lN;
  vector<IloRangeArray> _lNeq;

  vector<IloNumVarArray> _z;
  vector<IloNumVarArray> _yplus;
  vector<IloNumVarArray> _sd;
  vector<IloRangeArray> _yup;
  vector<IloRangeArray> _ydown;

 
  double _epsLS;
  vec _epsLSN;
  
  vec _xdes;

  vec _addCut;
  mat _addCutN1;
  mat _in;

};
#endif

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

  bool postLS(vec f, vec z, vec beta, vec SIGy,IloCplex * cplex, int iteration=0);
  bool postLSN1(int n,vec yn, vec zn,vec beta, vec sdn, IloCplex * cplex, int iteration=0);


 private:
  IloNumVar _riskLS0;
  IloRange _riskLS0Constraint;
  IloNumVarArray _riskLSN;
  IloRangeArray _riskLSNConstraint;

  vector<IloNumVarArray> _l;

 
  double _epsLS;
  vec _epsLSN;
  
  vec _xdes;

};
#endif

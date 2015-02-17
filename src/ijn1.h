#ifndef IJN1_H
#define IJN1_H

#include "ijcc.h"
#include "igrid.h"

using namespace std;


class ijn1 : public ijcc {

 public:
 ijn1(grid * gr,mat SIGy, mat Hw, double L, double p, double pc, double eps, double epsN): ijcc(gr,SIGy,L,p,pc,eps), _Hw(Hw),_epsN(epsN) { setup(); }
 ijn1(grid * gr,mat SIGy, mat Hw, vec beta, double L, double p, double pc, double eps, double epsN): ijcc(gr,SIGy,beta,L,p,pc,eps), _Hw(Hw),_epsN(epsN) { setup(); }
  void setup();
  rgrid * solveModel( isolve * is=NULL);
  bool postN1(int n,vec f, vec g,vec z, IloCplex * cplex, int iteration=0);
  vec getN1(int n, vec y0, vec g);
  vec getCheck(){ return _check; }
  double getNumCuts(int n){ return sum(_addCut.row(n)); }
  double getCutsLine(int i){ return sum(_addCut.col(i)); }
  double getTotalCuts();
  mat getVar(){ return _var; }

  class nowinner: public exception 
  {
    virtual const char* what() const throw()
    {
      return "No islanded bus found!";
    }
  }nowin;


 private:
  vector<IloNumVarArray> _z;
  vector<IloNumVarArray> _yplus;
  vector<IloRangeArray> _yup;
   vector<IloRangeArray> _ydown;
  IloRangeArray _riskConstraint;
  
  mat _addCut;

  mat _Hw;
  mat _L;
  mat _C;
  mat _Cg;
  mat _Hb;
  mat _var;
  vec _var0;
  double _epsN;

  vec _check;
  mat _in;

};
#endif

#ifndef ISJ_H
#define ISJ_H


#include <ctime>
#include "igrid.h"
#include "gridcalc.h"

using namespace std;


class isj : public igrid {

 public:
 isj(grid * gr, gridcalc * gc, mat SIGm, vec indexM,double L, double p, double pc, double eps, double epsG=1): igrid(gr), _gc(gc), _SIG(SIGm), _indexM(indexM) , _L(L), _p(p), _pc(pc), _eps(eps), _epsG(epsG) { setup(); }
  void setup();
  rgrid * solveModel( isolve * is=NULL);
  void lineLimitStatus(bool status);
  bool postCC(vec f, vec z, vec beta, vec SIGy,IloCplex * cplex, int iteration=0);
  mat getSIGm(){ return _SIG; }
  mat getA(){ return _gc->getH(); }
  mat getCb(){ return _gc->getC(); }  //line admittance matrix
  mat getCm(){ return _Cm; }          //volatile inject admittance matrix
  mat getCg(){ return _gc->getCm(); } //generator admittance matrix
  vec getIndexG(){ return _indexG; }
  double getEps(){ return _eps; }
  double getL(){ return _L; }
  double getP(){ return _p; }
  double getPc(){ return _pc; }
  gridcalc * getGC(){ return _gc; }
  IloNumVarArray getZ(){ return _z; }
  IloNumVarArray getYplus(){ return _yplus; }
  IloNumVarArray getBetaVar(){ return _beta; }
  IloRangeArray getGenUp(){ return _genup; }

  void setEps(double eps);

  double getSigDelta(){ return _sig_delta;}
  vec getSig(){ return _sig; };
  mat getSigger(){ return _sigger; }
  vec getBeta(){ return _betaSolve; }
  vec getSD(){ return _sdSolve; }
  void setBetaSolve(vec beta){ _betaSolve=beta; }
  void setSDSolve(vec sd){ _sdSolve=sd; }


  
  class iterlimit: public exception 
  {
    virtual const char* what() const throw()
    {
      return "Hit ITERATION LIMIT";
    }
  }itlimit;

 private:
  gridcalc * _gc;

  IloNumVar _risk;
  IloRange _riskConstraint;
  IloRange _betaSum;
  IloRange _budgetConstrain;
  IloNumVarArray _nu;
  IloNumVarArray _z;
  IloNumVarArray _yplus;
  IloNumVarArray _sd;
  IloNumVarArray _beta;
  IloRangeArray _yup;
  IloRangeArray _ydown;
  IloRangeArray _genup;
  IloRangeArray _gendown;

  mat _SIG;
  vec _indexM;
  vec _indexG;
  mat _Cm;
  double _L;
  double _p;
  double _pc;
  double _eps;
  double _epsG;

  mat _A;

  double _sig_delta;
  vec _sig;
  mat _sigger;
  
  vec _addCut;

  vec _betaSolve;
  vec _sdSolve;
  
};
#endif

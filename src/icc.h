#ifndef ICC_H
#define ICC_H

#include "igrid.h"
#include "gridcalc.h"

using namespace std;


class icc : public igrid {

 public:
 icc(grid * gr, gridcalc * gc, mat SIGm, vec indexM, double epsL, double epsG=1): igrid(gr), _gc(gc), _SIG(SIGm), _indexM(indexM), _epsL(epsL), _epsG(epsG) { setup(); }
  void setup();
  rgrid * solveModel( isolve * is=NULL);
  void lineLimitStatus(bool status);
  bool postCC(vec f,vec p, vec beta, vec SIGy,IloCplex * cplex, int iteration=0);
  mat getSIGm(){ return _SIG; }
  mat getA(){ return _gc->getH(); }
  mat getCb(){ return _gc->getC(); }  //line admittance matrix
  mat getCm(){ return _Cm; }          //volatile inject admittance matrix
  mat getCg(){ return _gc->getCm(); } //generator admittance matrix
  vec getIndexG(){ return _indexG; }
  double getEpsL(){ return _epsL; }
  gridcalc * getGC(){ return _gc; }
  IloNumVarArray getYplus(){ return _yplus; }
  IloNumVarArray getBetaVar(){ return _beta; }
  IloRangeArray getGenUp(){ return _genup; }

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

  class nanerror: public exception 
  {
    virtual const char* what() const throw()
    {
      return "NaN Error";
    }
  }nanerr;

 private:
  gridcalc * _gc;

  IloRange _betaSum;
  IloRange _budgetConstrain;
  IloNumVarArray _yplus;
  IloNumVarArray _sd;
  IloNumVarArray _beta;
  IloRangeArray _yup;
  IloRangeArray _ydown;
  IloRangeArray _lineprob;
  IloRangeArray _genup;
  IloRangeArray _gendown;

  mat _SIG;
  vec _indexM;
  vec _indexG;
  mat _Cm;
  double _epsL;
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

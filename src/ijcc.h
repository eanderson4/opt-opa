#ifndef IJCC_H
#define IJCC_H

#include "igrid.h"
#include "gridcalc.h"

using namespace std;


class ijcc : public igrid {

 public:
 ijcc(grid * gr,mat SIGy, double L, double p, double pc, double eps): igrid(gr), _SIGy(SIGy), _L(L), _p(p), _pc(pc), _eps(eps) { setup(); }
 ijcc(grid * gr,mat SIGy,vec beta, double L, double p, double pc, double eps): igrid(gr), _SIGy(SIGy), _beta(beta), _L(L), _p(p), _pc(pc), _eps(eps) { setup2(); }
  void setup();
  void setup2();
  rgrid * solveModel( isolve * is=NULL);
  void lineLimitStatus(bool status);
  bool postCC(vec f, vec z,IloCplex * cplex, int iteration=0);
  mat getSig(){ return _SIGy; }
  double getEps(){ return _eps; }
  double getL(){ return _L; }
  double getP(){ return _p; }
  double getPc(){ return _pc; }
  double getNumBaseCuts(){ return sum(_addCut); }
  double getBaseCutsLine(int i){ return _addCut(i); }
  IloNumVarArray getZ(){ return _z; }
  IloNumVarArray getFplus(){ return _fplus; }

  //Risk functions
  
  class itlimit: public exception 
  {
    virtual const char* what() const throw()
    {
      return "Stopped due to iteration limit";
    }
  }itlimit;

 private:
  IloNumVarArray _z;
  IloNumVarArray _fplus;
  IloRange _riskConstraint;
  IloRangeArray _fup;
  IloRangeArray _fdown;

  vec _addCut;

  mat _SIGy;
  vec _beta;
  double _L;
  double _p;
  double _pc;
  double _eps;

  IloNumVarArray _bt;
  
};
#endif

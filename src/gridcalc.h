#ifndef GRIDCALC_H
#define GRIDCALC_H

#include <iostream>
#include <stdlib.h>
#include <ilcplex/ilocplex.h>
#include <math.h>
#include "grid.h"
#include "igrid.h"
#include "rgrid.h"
#include "rv.h"
#include "del_g.h"

#define ARMA_DONT_USE_WRAPPER
#include <armadillo>


using namespace arma;
using namespace std;

class gridcalc
{

  public:
  gridcalc(grid * gr);
  ~gridcalc() {}

  vec getDelG(del_g dg);
  vec getDelF(vec delg, vec slack);
  mat getC(){ return _C; };  //Admittance Matrix
  mat getCm();  //GenAdmittance Matrix
  vec getIndexG();  //GenAdmittance Matrix
  mat getH(){ return _H; };  //Linear Shift Factor
  mat getHw(vec slack);      //with arbitray slack distribution
  mat getL(mat Hw);          //Line outage distribution factor
  void makeL(mat Hw){ _L = getL(Hw); }
  mat getMadeL(){ if(have_L) return _L;  else return NULL; }
  vec getSlack(){ return _slack; }
  //  vec getSlackDist(){ return getCm()*_slack; }
  vec convert(IloNumArray na);
  vec risk(vec f,vec varf, double L, double p, double pc);
  vec lineprob(vec f, vec varf);
  vec dz(double f,double varf, double L, double p, double pc);
  vec getN1(int n, vec y0, vec g,mat Hw);
  vec getD();

  class nowinner: public exception 
  {
    virtual const char* what() const throw()
    {
      return "No islanded bus found!";
    }
  }nowin;

 private:
  grid * _gr;

  mat _H;
  mat _L;
  mat _C;
  vec _slack;

  bool have_L;
  
 

};


#endif

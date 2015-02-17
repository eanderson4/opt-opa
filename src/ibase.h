#ifndef IBASE_H
#define IBASE_H

#include <ilcplex/ilocplex.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "grid.h"
#include "rgrid.h"

using namespace std;

class ibase {

 public:

 ibase() : _env( IloEnv() ) {  }
 ~ibase() {
   _env.end();
 }

 int buildModel( grid * gr ) ;
 void getBaseResults(IloCplex * cplex, rgrid * rg);

 IloEnv getEnv(){ return _env; }
 IloModel * getModel(){ return &_mod; }
 IloNumVarArray getG(){ return _g; }
 IloNumVarArray getF(){ return _f; }
 IloRangeArray getNodalBalance(){ return _nodalBalance; }
 IloRangeArray getBranchFlow(){ return _branchFlow; }
 IloRangeArray getPhaseAngle(){ return _phaseAngle; }

 private:
  
  IloNumVarArray _f; //branch flow
  IloNumVarArray _g; //generator dispatch
  IloNumVarArray _theta; //phase angle of nodes

  IloRangeArray _nodalBalance; //energy conservation at nodes
  IloRangeArray _branchFlow; //branch flow is bounded
  IloRangeArray _phaseAngle; //phase angle difference is bounded

  IloEnv _env;
  IloModel _mod;

};

class islack {

 public:
  islack(){}
 islack(IloNumArray g_nom,IloNumArray slack):_g_nom(g_nom),_slack(slack){}

  void buildSlack(grid * gf, IloModel * mod, IloNumVarArray g);
  void fixMismatch(grid * gr,IloNumVarArray g);

  IloNumArray getGNom(){ return _g_nom; }
  IloNumArray getSlack(){ return _slack; }
  IloNum getMismatch(){ return _mismatch; }

 private:
  IloNumArray _g_nom;
  IloNumArray _slack;

  IloNum _mismatch;
  IloRange _totalDemand;

};

class icost {
  
 public:
  icost(){}

  int buildCost(grid * gr, IloModel * mod, IloNumVarArray g);
  int buildCost(grid * gr, IloModel * mod, IloNumVarArray g, IloNumVarArray beta, double sigD);
  int buildCostWithLoadShed(grid * gr, IloModel * mod, IloNumVarArray g,IloNumVarArray ls);
  void constrainCost(grid * gr, IloModel * mod, IloNumVarArray g, double budget);
  double getCost(grid * gr, IloNumArray g);

 private:  
  IloRange _budgetConstraint;

};

class ished{

 public:
  ished(){}
  int buildLoadShed(grid * gr, IloModel * mod, IloRangeArray nb);
  IloNumVarArray getLS(){ return _ls; }
  void getLoadShed(IloCplex * cplex, rgrid * rg);

 private:
  IloNumVarArray _ls;

};


class isolve {

 public:
 isolve() : timeLimit( 1000 ) {}

  void setCplexParams(IloCplex * cplex){
    cplex->setParam( IloCplex::RootAlg, rootsolver ); //Auto, Primal, Dual, Barrier
    cplex->setParam( IloCplex::NodeAlg, nodesolver ); 
    cplex->setParam( IloCplex::TiLim, timeLimit );
  }

  void setSolver(IloCplex::Algorithm rootalg, IloCplex::Algorithm nodealg){
    rootsolver=rootalg;
    nodesolver=nodealg;
  }
  void setTime(int tl){
    timeLimit=tl;
  }

 private:

  IloCplex::Algorithm rootsolver;
  IloCplex::Algorithm nodesolver;
  int timeLimit;

};


#endif

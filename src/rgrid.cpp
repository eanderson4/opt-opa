#include "rgrid.h"

void rgrid::getSolveInfo(IloCplex * cplex, double rt){
  realTime=rt;
  status = cplex->getCplexStatus();
  objective = cplex->getObjValue();
  alg = cplex->getAlgorithm();
  iterations = cplex->getNiterations();
  rootalg=cplex->getParam(IloCplex::RootAlg);
  if(rootalg==IloCplex::Primal) rootname="Primal";
  if(rootalg==IloCplex::Dual) rootname="Dual";
  if(rootalg==IloCplex::Barrier) rootname="Barrier";
  
  isQO = cplex->isQO();
}

void rgrid::setGenCost(double genCost){
  have_gencost=true;
  _genCost=genCost;
}

void rgrid::outputInfo(ostream & out){
  if(status){
    out<<"Time: "<<realTime<<", Objective: "<<objective<<", Algorithm: "<<alg<<","<<rootname<<", it: "<<iterations<<endl;
    if(_loadshed>0)
      out<<"    -ls-    Load Shed: "<<_loadshed<<endl;
    if(have_gencost){
      out<<"    -gen-   Gen Cost: "<<_genCost;
      if(isQO) out<<"    -quad obj-";
      out<<"\n";
    }
  }
  else{
    out<<"Problem didn't solve"<<endl;
  }
} 


void rgrid::displayOperatingPos(grid * gr){
  if(status){
    displayBranchPos(gr);
    displayGenPos(gr);
  }
  else{
    cerr<<"No solve info"<<endl;
  }
}

void rgrid::displayGenPos(grid * gr){

  int nGen=gr->numGens();
  double totalG;
  double totalCapacity;
  double excess;

  cout<<"Generator Position"<<endl;
  for(int i=0;i<nGen;i++){
    gen gi = gr->getGen(i);
    double pmin = gi.getPmin();
    double pmax = gi.getPmax();
    double g = getG()[i];
    totalG += g;
    totalCapacity += pmax;
    cout<<i<<": "<<g<<" in [ "<<pmin<<", "<<pmax<<"]"<<endl;
  }
  excess = totalCapacity - totalG;
  cout<<"Gen: "<<totalG<<", Cap: "<<totalCapacity<<", Excess: "<<excess<<endl;

}


void rgrid::displayBranchPos(grid * gr){

  int nBranch=gr->numBranches();
  double totalFlow;
  double totalCapacity;
  double excess;

  cout<<"Branch Position"<<endl;
  for(int i=0;i<nBranch;i++){
    branch bi = gr->getBranch(i);
    double U = bi.getRateA();
    double f = getF()[i];
    totalFlow += fabs(f);
    totalCapacity += U;
    cout<<i<<": "<<f<<" in [ "<<-U<<", "<<U<<"]"<<endl;
  }
  excess = totalCapacity - totalFlow;
  cout<<"Flow: "<<totalFlow<<", Cap: "<<totalCapacity<<", Excess: "<<excess<<endl;

}

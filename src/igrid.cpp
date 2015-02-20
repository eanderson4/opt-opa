#include "igrid.h"

rgrid *  igrid::solveModel( isolve * is){

  rgrid * rg = new rgrid();
  
  time_t tstart;
  tstart = clock();
    
  IloCplex cplex(*getModel());
  
  if(is!=NULL) is->setCplexParams(&cplex);
  
  if (cplex.solve()){
      float tot= float(clock() - tstart) / CLOCKS_PER_SEC;
      cout<<"\n - Solve info -"<<endl;
      cout<<"MODEL solved in "<<tot<<endl;
      rg->getSolveInfo(&cplex,tot);
      getBaseResults(&cplex, rg);
      if(have_loadshed) getIshed().getLoadShed(&cplex, rg);

      cout<<"STATUS: "<<rg->getStatus()<<endl;
      cout<<"OBJECTIVE: "<<cplex.getObjValue()<<"\n"<<endl;
      double genCost = getIcost().getCost(_gr,rg->getG());
      rg->setGenCost(genCost);
      
  }
  else{
    cerr<<"Not solved"<<endl;
    cerr<<cplex.getStatus()<<endl;
  }
  cplex.end();
  
  return rg;
}

rgrid *  igrid::solveModel(IloCplex * cplex){

  rgrid * rg = new rgrid();
  
  time_t tstart;
  tstart = clock();
  
  if (cplex->solve()){
      float tot= float(clock() - tstart) / CLOCKS_PER_SEC;
      cout<<"\n - Solve info -"<<endl;
      cout<<"MODEL solved in "<<tot<<endl;
      rg->getSolveInfo(cplex,tot);
      getBaseResults(cplex, rg);
      if(have_loadshed) getIshed().getLoadShed(cplex, rg);

      cout<<"STATUS: "<<rg->getStatus()<<endl;
      cout<<"OBJECTIVE: "<<cplex->getObjValue()<<"\n"<<endl;
      double genCost = getIcost().getCost(_gr,rg->getG());
      rg->setGenCost(genCost);
      
  }
  else{
    cerr<<"Not solved"<<endl;
    cerr<<cplex->getStatus()<<endl;
  }

  
  return rg;
}


void igrid::allowLoadShed(){
  have_loadshed = _ils.buildLoadShed(_gr,getModel(),getNodalBalance());
  if(have_loadshed) cout<<"Built load shed variables"<<endl;
}

int igrid::addCost(){ 
  if(have_loadshed){
    have_cost=_ic.buildCostWithLoadShed(_gr,getModel(),getG(),_ils.getLS());
  }
  else{
    have_cost=_ic.buildCost(_gr,getModel(),getG());
  }
  if(have_cost) cout<<"Built cost object"<<endl;
  return have_cost; 
}

int igrid::addCost(IloNumVarArray beta, double sigD){ 
  have_cost=_ic.buildCost(_gr,getModel(),getG(),beta,sigD);
  if(have_cost) cout<<"Built cost object with Arbitrary slack"<<endl;
  return have_cost; 
}

void igrid::addSlack(IloNumArray g_nom, IloNumArray slack){
  _islk = islack(g_nom,slack);
  
  _islk.buildSlack(_gr,getModel(),getG());

  have_slack=true;
}


void igrid::modGrid( del_g mod ){
  _mod=mod;

  //  make modification
  stringstream ss;
  int nB = _gr->numBranches();

  if(mod.haveTopo()){
    
    for(int i=0; i<nB; i++){
      if(mod.getStatus(i)==false){
	// line i is outaged
	cout<<"Line "<<i<<" - ";
	getF()[i].setBounds(0,0);
	getBranchFlow()[i].setBounds(-IloInfinity,IloInfinity );
	getPhaseAngle()[i].setBounds(-IloInfinity,IloInfinity );
	cout<<"f=0, phase angles relaxed"<<endl;
	ss.str("");
	ss<<"f"<<i<<"[0,0]";
	getF()[i].setName( ss.str().c_str() );
	cout<<getF()[i]<<endl;
	cout<<getBranchFlow()[i]<<endl;
	cout<<getPhaseAngle()[i]<<endl;
      }
    }
    
  }

  if(mod.haveDemand()){
    for(int i=0; i<_gr->numBuses(); i++){
      double del_demand = mod.getDemand(i);
      if(del_demand != 0){
	IloRange nb = getNodalBalance()[i];
	nb.setExpr(nb.getExpr() + del_demand);
      }
    }
  }

  //MODIFYING GRID TOO

  if(mod.haveTopo()){
    int nB=_gr->numBranches();
    for(int i=0; i<nB; i++){
      if(mod.getStatus(i)==false){
	// line i is outaged
	cout<<"Line "<<i<<" - ";
	_gr->getBranch(i).setStatus(0);	
      }
    }
  }
  if(mod.haveDemand()){
    cout<<"Modify Demand (add)"<<endl;
    int nB=_gr->numBuses();
    for(int i=0;i<nB; i++){
      if(mod.getDemand(i) != 0){
	//	cout<<i<<": "<<_gr->getBus(i).getP()<<" + "<<mod.getDemand(i)<<endl;
	_gr->addPd(i,mod.getDemand(i));
      }
    }
  }


  if(have_slack){
    cout<<"Fix Mismatch"<<endl;
    cout<<"Total Demand: "<<_gr->getTotalDemand()<<endl;
    _islk.fixMismatch(_gr, getG());
  }


}

void igrid::unmodGrid( del_g mod ){
  _mod=mod;

  //  make modification
  stringstream ss;
  int nB = _gr->numBranches();

  if(mod.haveTopo()){
    
    for(int i=0; i<nB; i++){
      if(mod.getStatus(i)==false){
	// line i was outaged repair
	cout<<"Line "<<i<<" - ";
	double U = _gr->getBranch(i).getRateA();
	getF()[i].setBounds(-U,U);
	getBranchFlow()[i].setBounds(0,0);
	getPhaseAngle()[i].setBounds(-360,360 );
	cout<<"f unbounded, phase angles tight"<<endl;
	ss.str("");
	ss<<"f"<<i<<"["<<-U<<","<<U<<"]";
	getF()[i].setName( ss.str().c_str() );
	cout<<getF()[i]<<endl;
	cout<<getBranchFlow()[i]<<endl;
	cout<<getPhaseAngle()[i]<<endl;
      }
    }
    
  }

  if(mod.haveDemand()){
    for(int i=0; i<_gr->numBuses(); i++){
      double del_demand = mod.getDemand(i);
      if(del_demand != 0){
	IloRange nb = getNodalBalance()[i];
	nb.setExpr(nb.getExpr() - del_demand);
      }
    }
  }

  //MODIFYING GRID TOO

  if(mod.haveTopo()){
    int nB=_gr->numBranches();
    for(int i=0; i<nB; i++){
      if(mod.getStatus(i)==false){
	// line i was outaged, repair
	cout<<"Line "<<i<<" - ";
	_gr->getBranch(i).setStatus(1);	
      }
    }
  }
  if(mod.haveDemand()){
    cout<<"Modify Demand (remove)"<<endl;
    int nB=_gr->numBuses();
    for(int i=0;i<nB; i++){
      if(mod.getDemand(i) != 0){
	//	cout<<i<<": "<<_gr->getBus(i).getP()<<" - "<<mod.getDemand(i)<<endl;
	_gr->addPd(i,-mod.getDemand(i));
      }
    }
  }


  if(have_slack){
    cout<<"fix Mismatch"<<endl;
    cout<<_gr->getTotalDemand()<<endl;
    _islk.fixMismatch(_gr, getG());
  }


}




void igrid::removeDemand(del_g mod){
  int nB=_gr->numBuses();
  //IMPLEMENTATION
  for(int i=0; i<nB; i++){
    double del_demand = mod.getDemand(i);
    if(del_demand != 0){
      IloRange nb = getNodalBalance()[i];
      nb.setExpr(nb.getExpr() - del_demand);
    }
  }
  
  //GRID
  //  cout<<"Modify Demand (remove)"<<endl;
  for(int i=0;i<nB; i++){
    if(mod.getDemand(i) != 0){
      //      cout<<i<<": "<<_gr->getBus(i).getP()<<" + "<<mod.getDemand(i)<<endl;
      _gr->addPd(i,-mod.getDemand(i));
    }
  }
  
  
}

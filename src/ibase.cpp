#include "ibase.h"

int ibase::buildModel(grid * gr ) {
  stringstream ss;
    int nk = gr->numBranches();
    int nb = gr->numBuses();
    int ng = gr->numGens();

    _mod = IloModel(_env);

    _f = IloNumVarArray(_env,nk,-IloInfinity,IloInfinity);
    _g = IloNumVarArray(_env,ng,0,IloInfinity);
    _theta = IloNumVarArray(_env,nb,-IloInfinity,IloInfinity);

    _nodalBalance = IloRangeArray(_env, nb, 0, 0);
    _branchFlow = IloRangeArray(_env, nk, 0, 0);
    _phaseAngle = IloRangeArray(_env, nk, 0, 0);


    for(int i=0; i<nb; i++){
      bus bi = gr->getBus(i);
      int type = bi.getType();
      double pd =bi.getP();
      double gs = bi.getGs();

      _nodalBalance[i].setExpr( IloExpr(_env,pd + gs));
      ss.str("");
      ss<< "nb" << i;
      _nodalBalance[i].setName( ss.str().c_str() );

      ss.str("");
      ss<<"theta"<<i;
      if (type==3) {
	_theta[i].setBounds(0,0);
	ss<<"[0,0]";
      }
      else {
	_theta[i].setBounds(-360,360);
	ss<<"[-360,360]";
      }

      _theta[i].setName( ss.str().c_str() );

      	//_nodalBalance[i].setBounds(-IloInfinity,IloInfinity);
    }

    for(int k=0; k<nk; k++){
      branch bk = gr->getBranch(k);
      int from=gr->getFromBus(k);
      int to=gr->getToBus(k);
      double x = bk.getX();
      double tap = bk.getTap();
      double b;
      if (tap == 0) b = 1/x;
      else b = 1/(x*tap);
      double shift = bk.getShift();
      double limit = bk.getRateA();

      _branchFlow[k].setExpr( IloExpr( -_f[k] + b * (_theta[from] - _theta[to] - shift) ));
      ss.str("");
      ss<<"bf"<<k;
      _branchFlow[k].setName( ss.str().c_str() );

      _f[k].setBounds(-limit,limit);
      ss.str("");
      ss<<"f"<<k; //<<"["<<-limit<<","<<limit<<"]";
      _f[k].setName( ss.str().c_str() );


      _nodalBalance[from].setExpr(_nodalBalance[from].getExpr() + _f[k]);
      _nodalBalance[to].setExpr(_nodalBalance[to].getExpr() - _f[k]);

      double angleMin = bk.getAngMin();
      double angleMax = bk.getAngMax();
      _phaseAngle[k].setExpr( IloExpr( _theta[from] - _theta[to]));
      _phaseAngle[k].setBounds(angleMin,angleMax);     
      ss.str("");
      ss<<"pa"<<k;
      _phaseAngle[k].setName( ss.str().c_str() );
    }


    for(int j=0; j<ng; j++){
      gen gj = gr->getGen(j);
      //      double pd = gj.getP();
      double pmin = gj.getPmin();
      double pmax = gj.getPmax();
      double status = gj.getStatus();
      int bus = gr->getBusNum(gj.getBus());
      _nodalBalance[bus].setExpr(_nodalBalance[bus].getExpr() - _g[j]);
      if( status <= 0) _g[j].setBounds(0,0);  //gen out of service
            else _g[j].setBounds(0,pmax);  //DONT RESPECT LOWER LIMIT
      //            else _g[j].setBounds(pmin,pmax);  //gen in service respect gen limits
      ss.str("");
      ss<<"g"<<j;
      _g[j].setName( ss.str().c_str() );
      //      cout<<_nodalBalance[bus]<<endl;
    }

    _mod.add(_f);
    _mod.add(_g);
    _mod.add(_theta);
    _mod.add(_nodalBalance);
    _mod.add(_branchFlow);
    _mod.add(_phaseAngle);

    return 1;
}


void ibase::getBaseResults(IloCplex * cplex, rgrid * rg){

  cplex->getValues(rg->getF(),_f);
  cplex->getValues(rg->getG(),_g);
  cplex->getValues(rg->getTheta(), _theta);
  
}

void islack::buildSlack(grid * gr, IloModel * mod, IloNumVarArray g){
  
  IloEnv env=mod->getEnv();
  int nG=gr->numGens();

  _mismatch = IloSum(_g_nom) - gr->getTotalDemand();

  cout<<"Building fixed slack model\n"
      <<"Current mismatch: "<<_mismatch<<endl;
  

  for(int i=0;i<nG;i++){
    IloNum gset=_g_nom[i]-_slack[i]*_mismatch;
    g[i].setBounds(gset,gset);
  }
  double td=gr->getTotalDemand();
  _totalDemand = IloRange(env,td,td,"totaldemand");
  _totalDemand.setExpr( IloSum(g) );
  mod->add(_totalDemand);
  


}

void islack::fixMismatch(grid * gr, IloNumVarArray g){
  int nG=gr->numGens();
  double td=gr->getTotalDemand();
  cout<<_g_nom<<endl;
  _mismatch = IloSum(_g_nom) - td;

  cout<<"mismatch: "<<_mismatch<<endl;
  

  for(int i=0;i<nG;i++){
    IloNum gset=_g_nom[i]-_slack[i]*_mismatch;
    g[i].setBounds(gset,gset);
  }
  _totalDemand.setBounds(td,td);
  cout<<_totalDemand<<endl;

}


int icost::buildCost(grid * gr, IloModel * mod, IloNumVarArray g){
  int ng = gr->numGens();
  IloEnv env = mod->getEnv();
  IloExpr exp(env);
  bool quadratic=false;
  for(int i =0; i<ng; i++){
    gen gi = gr->getGen(i);
    double c2 = gi.getC2();
    double c1 = gi.getC1();
    double c0 = gi.getC0();
    exp += g[i]*g[i]*c2 + c1*g[i] + c0;
    if(c2>0) quadratic=true;
  }
  mod->add(IloMinimize(env,exp));

  if (quadratic) cout<<"The objective is quadratic"<<endl;
  return 1;
}
int icost::buildCost(grid * gr, IloModel * mod, IloNumVarArray g, IloNumVarArray beta, double sigD){
  int ng = gr->numGens();
  IloEnv env = mod->getEnv();
  IloExpr exp(env);
  bool quadratic=false;
  for(int i =0; i<ng; i++){
    gen gi = gr->getGen(i);
    double c2 = gi.getC2();
    double c1 = gi.getC1();
    double c0 = gi.getC0();
    exp += (g[i]*g[i] + beta[i]*beta[i]*sigD)*c2 + c1*g[i] + c0;
    if(c2>0) quadratic=true;
  }
  mod->add(IloMinimize(env,exp));
  cout<<"Arbitrary Slack Model"<<endl;
  if (quadratic) cout<<"The objective is quadratic"<<endl;
  return 1;
}

int icost::buildCostWithLoadShed(grid * gr, IloModel * mod, IloNumVarArray g, IloNumVarArray ls){
  double lsPenalty=gr->getLoadShedPenalty();

  int ng = gr->numGens();
  IloEnv env = mod->getEnv();
  IloExpr exp(env);
  for(int i =0; i<ng; i++){
    gen gi = gr->getGen(i);
    double c2 = gi.getC2();
    double c1 = gi.getC1();
    double c0 = gi.getC0();
    exp += g[i]*g[i]*c2 + c1*g[i] + c0;
  }
  int nb = gr->numBuses();
  for(int i=0; i<nb; i++){
    exp += lsPenalty*ls[i];
  }
  mod->add(IloMinimize(env,exp));
  exp.end();
  return 1;
}

void icost::constrainCost(grid * gr, IloModel * mod, IloNumVarArray g, double budget){
  int ng = gr->numGens();
  IloEnv env = mod->getEnv();
  IloExpr exp(env);
  for(int i =0; i<ng; i++){
    gen gi = gr->getGen(i);
    double c2 = gi.getC2();
    double c1 = gi.getC1();
    double c0 = gi.getC0();
    exp += g[i]*g[i]*c2 + c1*g[i] + c0;
  }
  _budgetConstraint = IloRange(env,0,budget);
  _budgetConstraint.setExpr(exp);
  mod->add(_budgetConstraint);
  exp.end();
}

double icost::getCost(grid * gr, IloNumArray g){
  int ng = gr->numGens();
  double cost;

  cost=0;
  for(int i =0; i<ng; i++){
    gen gi = gr->getGen(i);
    double c2 = gi.getC2();
    double c1 = gi.getC1();
    double c0 = gi.getC0();
    cost=cost+c2*g[i]*g[i] + c1*g[i]+c0;
  }
  return cost;

}


int ished::buildLoadShed(grid * gr, IloModel * mod, IloRangeArray nb){

  int checknb = gr->numBuses();  
  if(nb.getSize() != checknb){
    cerr<<"Dimension error on number of buses"<<endl;
    return 0;
  }

  IloEnv env = mod->getEnv();
  _ls = IloNumVarArray(env, checknb, 0, IloInfinity);

  for(int i=0; i<checknb; i++){
    bus bi = gr->getBus(i);
    double pd =bi.getP();
    double gs = bi.getGs();
    
    nb[i].setExpr( nb[i].getExpr() - _ls[i] );
    //    _ls[i].setBounds(0,(pd+gs));
  }

   
  return 1;
}

void ished::getLoadShed(IloCplex * cplex, rgrid * rg){
  cplex->getValues(rg->getLS(),_ls);
  double loadshed = IloSum(rg->getLS());
  rg->setLoadShed(loadshed);
}

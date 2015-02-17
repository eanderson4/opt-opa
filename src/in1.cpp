#include "in1.h"

rgrid *  in1::solveModel( isolve * is){

  rgrid * rg = new rgrid();

  time_t tstart;
  tstart = clock();
    
  int Nl = getGrid()->numBranches();
  gridcalc gc(getGrid());
  
  IloCplex cplex(*getModel());
  
  if(is!=NULL) is->setCplexParams(&cplex);
  int n=0;
  double large=0;

  //  cerr<<"n\tyn\tzn\tan\tbn\tun\tvn\n";
  cout<<"Check: "<<sum(_check)<<" / "<<getGrid()->numBranches()<<endl;

  if (cplex.solve()){
    bool systemfail=true;
    while(systemfail){
      n++; if (n>500) throw itlimit;  
      cout<<" ----- Iteration - "<<n<<" ----- "<<endl;

      IloNumArray fsolve(getEnv(),Nl);
      cplex.getValues(fsolve,getF());
      vec f0=gc.convert(fsolve);

      IloNumArray gsolve(getEnv(),Nl);
      cplex.getValues(gsolve,getG());
      vec g0=gc.convert(gsolve);
      systemfail=false;
      
      if(n==2) large=4000;
      for(int i=0;i<Nl;i++){
	if(_check(i) && i<(30 + large)){
	  cout<<"Contingency: "<<i<<endl;
	  vec fn = getN1(i,f0,g0);
	  if(!fn.is_finite()) fn.print("fn: ");
	  cout<<"Post eval"<<endl;
	  systemfail += postN1(i,fn,&cplex);
	}
	else cout<<"dont check Contingency "<<i<<endl;
      }

      
      if(!systemfail && n>=2) break;
      if(!cplex.solve()) break;

    }

    //Risk constraint satisfied, record solution
    float total= float(clock() - tstart) / CLOCKS_PER_SEC;  
    cout<<"\n - Solve info -"<<endl;
    cout<<"MODEL solved in "<<total<<endl;
    rg->getSolveInfo(&cplex,total);
    getBaseResults(&cplex, rg);
    if(isLoadShed()) getIshed().getLoadShed(&cplex, rg);
    
    cout<<"STATUS: "<<rg->getStatus()<<endl;
    cout<<"OBJECTIVE: "<<cplex.getObjValue()<<"\n"<<endl;
    double genCost = getIcost().getCost(getGrid(),rg->getG());
    rg->setGenCost(genCost);
    
    cout<<"Iterations: "<<n<<endl;
    cout<<"Cuts: "<<addedBounds<<endl;
    cout<<"Time: "<<total<<endl;
        
  }

  else{
    cerr<<"Not solved"<<endl;
    cerr<<cplex.getStatus()<<endl;
  }
  cplex.end();
  
  return rg;
}

bool in1::postN1(int n, vec yn, IloCplex * cplex){
  if(_check(n)!=1) return false;

  stringstream ss;
  int Nl = getGrid()->numBranches();

  bool fail=false;
  for(int i=0;i<Nl;i++){
    double U=getGrid()->getBranch(i).getRateA();
    double Up=_m1*U;
    if(yn(i)>Up){
      cout<<"Failed"<<endl;
      fail=true;
      if(sum(_in.row(n))==0){
	_yplus[n] = IloNumVarArray(getEnv(),Nl,0,IloInfinity);	  
	_yup[n] = IloRangeArray(getEnv(),Nl,0,IloInfinity);	  
	_ydown[n] = IloRangeArray(getEnv(),Nl,0,IloInfinity);	  
      }
      _in(n,i)=1;
      addedBounds++;

      _yup[n][i].setExpr( _yplus[n][i] - (getF()[i] + _L(i,n)*(getF()[n])) );
      _ydown[n][i].setExpr( _yplus[n][i] + (getF()[i] + _L(i,n)*getF()[n]) );

      _yplus[n][i].setBounds(0, Up);
      ss.str("");
      ss<<"y"<<n<<","<<i<<"[0,"<<Up<<"]";
      _yplus[n][i].setName( ss.str().c_str() );
      
      getModel()->add(_yplus[n][i]);
      getModel()->add(_yup[n][i]);
      getModel()->add(_ydown[n][i]);
      
    }
  }

  return fail;
}




void in1::setup(){
  cout<<"Setup N-1 standard formulation"<<endl;
  stringstream ss;

  
  IloEnv env = getEnv();

  int Nl = getGrid()->numBranches();

  _check = vec(Nl,fill::ones);

  _C = _gc->getC();
  _Cg = _gc->getCm();
  cout<<"Get Branch Sensitivities"<<endl;
  _Hb=_gc->getH()*trans(_C);
  cout<<"Build L"<<endl;
  _L=_gc->getL(_Hw);

  _in = mat(Nl,Nl,fill::zeros);
  addedBounds=0;

  cout<<"Here"<<endl;;
  mat Sig = _SIGy;
  _var0 = Sig.diag();
  _var=mat(Nl,Nl);
  
  for(int n=0;n<Nl;n++){  //Line i outage
    for(int j=0;j<Nl;j++){
      if(_Hb(n,n)<=1-.000001 || _Hb(n,n)>=1+.0000001){
	_var(n,j)=Sig(j,j) + 2*_L(j,n)*Sig(n,j) + pow(_L(j,n),2)*Sig(n,n);
      }
      else{
	//	_var(n,j)=Sig(j,j) + 2*_Hb(j,n)*Sig(n,j) + pow(_Hb(j,n),2)*Sig(n,n);
	_check(n)=0;
      }
    }
  }
  cout<<sum(_check)<<" / "<<getGrid()->numBranches()<<endl;
  cout<<"there"<<endl;


  for(int n=0;n<Nl;n++){
    _yplus.push_back(IloNumVarArray(env,1,0,IloInfinity));
    _yup.push_back(IloRangeArray(env, 1,0,IloInfinity));
    _ydown.push_back(IloRangeArray(env,1,0,IloInfinity));

  }
  cout<<"DONT BUILDING"<<endl;
}




void in1::setupOld(){
  cout<<"Setup N-1 standard formulation"<<endl;
  stringstream ss;

  
  IloEnv env = getEnv();

  int Nl = getGrid()->numBranches();
  gridcalc gc(getGrid());

  _check = vec(Nl,fill::ones);

  _C = gc.getC();
  _Cg = gc.getCm();
  cout<<"Get Branch Sensitivities"<<endl;
  _Hb=_Hw*trans(_C);
  cout<<"Build L"<<endl;
  _L=gc.getL(_Hw);

  _in = mat(Nl,Nl,fill::zeros);
  addedBounds=0;

  cout<<"Here"<<endl;;
  mat Sig = _SIGy;
  _var0 = Sig.diag();
  _var=mat(Nl,Nl);
  
  for(int n=0;n<Nl;n++){  //Line i outage
    for(int j=0;j<Nl;j++){
      if(_Hb(n,n)<=1-.000001 || _Hb(n,n)>=1+.0000001){
	_var(n,j)=Sig(j,j) + 2*_L(j,n)*Sig(n,j) + pow(_L(j,n),2)*Sig(n,n);
      }
      else{
	//	_var(n,j)=Sig(j,j) + 2*_Hb(j,n)*Sig(n,j) + pow(_Hb(j,n),2)*Sig(n,n);
	_check(n)=0;
      }
    }
  }
  cout<<sum(_check)<<" / "<<getGrid()->numBranches()<<endl;
  cout<<"there"<<endl;


  for(int n=0;n<Nl;n++){
    _yplus.push_back(IloNumVarArray(env,1,0,IloInfinity));
    _yup.push_back(IloRangeArray(env, 1,0,IloInfinity));
    _ydown.push_back(IloRangeArray(env,1,0,IloInfinity));

  }
  cout<<"DONT BUILDING"<<endl;
}


vec in1::getN1(int n, vec y0, vec g){
  grid * gr = getGrid();
  vec yn;

  if(_Hb(n,n)<=1-.000001 || _Hb(n,n)>=1+.0000001){
    //    cout<<"Line "<<n<<endl;
    yn = y0[n]*_L.col(n)+y0;
  }
  else{
    cout<<"PROBLEM"<<endl;
  }
  
  return yn;

}

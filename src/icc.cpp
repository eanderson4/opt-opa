#include "icc.h"

rgrid *  icc::solveModel( isolve * is){

  rgrid * rg = new rgrid();

  double tol = pow(10,-4);
  time_t tstart;
  tstart = clock();
    
  int Nl = getGrid()->numBranches();
  int Ng = getGrid()->numGens();
  int Nm = _Cm.n_cols;
  
  IloCplex cplex(*getModel());
  
  if(is!=NULL) is->setCplexParams(&cplex);
  int n=0;  

  if (cplex.solve()){
    bool systemfail=true;

    mat Am = getA()*getCm();
    vec ones(Nm,1,fill::ones);
    mat ker = ones.t()*_SIG*ones;
    mat AmSone = Am*_SIG*ones;
    mat Ag = getA()*getCg();
    while(systemfail){
      n++; if (n>100) throw itlimit;
      if(cplex.getStatus() == IloAlgorithm::Infeasible){
	cout<<"Infeasible"<<endl;
	break;
      }
      cout<<"Iteration: "<<n<<endl;
      cout<<"Checking Line Probabilities Constraint"<<endl;

      IloNumArray betasolve(getEnv(),Ng);
      IloNumArray ysolve(getEnv(),Nl);
      IloNumArray sdsolve(getEnv(),Nl);

      cplex.getValues(ysolve,getF());
      cplex.getValues(betasolve,_beta);
      cplex.getValues(sdsolve,_sd);

      vec y=_gc->convert(ysolve);      
      vec beta=_gc->convert(betasolve);      
      vec sds=_gc->convert(sdsolve);      
      
      //      beta.t().print("beta: ");
      //      sds.t().print("sd: ");

      mat term3=Ag*beta;
      mat SIGy = term3*ker*term3.t()-term3*AmSone.t()-AmSone*term3.t()+_sigger;
      //      (term3*ker*term3.t()).print("one: ");
      //      (term3*AmSone.t()).print("two: ");
      //      _sigger.print("sig: ");

      vec p = _gc->lineprob(y,SIGy.diag());

      //      p.t().print("lineprob: ");

      double maxp=p.max();

      for(int i=0;i<Nl;i++){
	if(p(i)>0){
	  double U=getGrid()->getBranch(i).getRateA();
	  cout<<i<<":\t"<<abs(y(i))<<"\t"<<U<<"\t+/- "<<sqrt(SIGy(i,i))<<"\t"<<p(i)<<endl;
	}
      }

      if(maxp <= _epsL+5*tol){
	systemfail = false;
      }
      else{    
	systemfail = true;
	postCC(y,p,beta,SIGy.diag(),&cplex);
      }

      if(systemfail) {
	cout<<"Re-Solving LP"<<endl;
	cplex.solve();
	cout<<"LP Solved, perform risk analysis"<<endl;
      }
      else{
	_betaSolve = beta;
	_sdSolve = SIGy.diag();
      }
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
    
    cout<<"\n\nRisk constraint satisfied\n\n"<<endl;
    cout<<"Iterations: "<<n<<endl;
    cout<<"Time: "<<total<<endl;
    //    _addCut.t().print("cuts: ");
  }

  else{
    cerr<<"Not solved"<<endl;
    cerr<<cplex.getStatus()<<endl;
  }
  cplex.end();
  
  return rg;

}      

bool icc::postCC(vec y, vec p, vec beta, vec SIGy,IloCplex * cplex, int iteration){
  stringstream ss;
  //define tolerance for line risk > 0
  double tol = pow(10,-4);
  int Nl = getGrid()->numBranches();

  
  ranvar rv;
  double nu = rv.PHIInverse(1-_epsL);
  vec pi = getA()*getCg()*beta;
  cout<<"CUTTING ----------"<<endl;
  for(int i=0; i<Nl; i++){
    if (p(i)>_epsL+tol){
      cout<<"Line "<<i<<" cuts"<<endl;
      double U=getGrid()->getBranch(i).getRateA();
      if(_addCut(i)==0){
	cout<<"Initialize Cutting Variables for Line "<<i<<endl;
	_lineprob[i].setExpr( U - _yplus[i] - _sd[i]*nu );
	//	cout<<_lineprob[i]<<endl;
	_yup[i].setExpr( _yplus[i] - getF()[i] );
	_ydown[i].setExpr( _yplus[i] + getF()[i] );
	_yplus[i].setBounds(0,U);	  
	getModel()->add(_yup[i]);
	getModel()->add(_ydown[i]);
	getModel()->add(_lineprob[i]);
	ss.str("");
	ss<<"yplus"<<i<<"[0,"<<U<<"]";
	_yplus[i].setName( ss.str().c_str() );
      }
      _addCut(i)=_addCut(i)+1;
      
      
      //Add cuts to describe standard deviation of branch flow
      IloRange cut_sd(getEnv(),-IloInfinity,0);
      
      int Ng = getGrid()->numGens();
      
      double pi_i = pi(i);
      double sd_i = sqrt(SIGy(i));
      
      cut_sd.setExpr( sd_i - _sd[i] );
      double term = (pi_i*_sig_delta - _sig(i))/sd_i;
      //	cout<<"\n";
      for( int j=0; j<Ng;j++){
	//double pf_j = term;
	double pf_j = _A(i,_indexG(j))*term;
	cut_sd.setExpr( cut_sd.getExpr() + pf_j*(_beta[j]-beta(j)) );
      }
      //      	cout<<cut_sd<<endl;
      getModel()->add(cut_sd);
      
      //	cout<<"\n\n";
    }
  }
  return true;
}     




void icc::lineLimitStatus(bool status){
  int Nl = getGrid()->numBranches();
  
  ranvar rv;
  //  double Ueps = rv.ginv(_eps,_L,_p,_pc);

  for(int i=0;i<Nl;i++){
    double U = getGrid()->getBranch(i).getRateA();
    if(!status){
      //      getF()[i].setBounds(-U*Ueps,U*Ueps);
      getF()[i].setBounds(-IloInfinity,IloInfinity);
    }
    else{
      double U = getGrid()->getBranch(i).getRateA();
      getF()[i].setBounds(-U,U);
    }
  }


}


void icc::setup(){
  cout<<"--cc-- ===Setup=== slack chance constraint"<<endl;
  stringstream ss;

  ranvar rv;

  IloEnv env = getEnv();
  int Nl = getGrid()->numBranches();
  int Nb = getGrid()->numBuses();
  int Ng = getGrid()->numGens();
  int Nm = _indexM.n_elem;

  mat Cm(Nb,Nm,fill::zeros);

  for(int i=0;i<Nm;i++){
    Cm(_indexM(i),i)=1;
  }
  _Cm = Cm;

  _indexG = _gc->getIndexG();

  _addCut = vec(Nl,fill::zeros);
    _A = getA();
  //Calculate branch variance terms --------------------
  cout<<"Build Static Variance Terms: "<<endl;

  mat Ag = _A*getCg();
  mat Ak = _A*getCm();
  vec ones2(Nm,1,fill::ones);
  _sig_delta = accu(_SIG);
  _sig = Ak*_SIG*ones2;
  _sigger = Ak*_SIG*Ak.t();
      
  if(!_sigger.is_finite()){
    Ag.print("Ag: ");
    Ak.print("Ak: ");
    getCm().print("Cm: ");
    cin>>Nm;
    throw nanerr;
  }

  //Build CPLEX model ----------------------------------
  cout<<"Build CPLEX model for base: "<<endl;

  _yplus = IloNumVarArray(env,Nl,0,IloInfinity);
  _sd = IloNumVarArray(env,Nl,0,IloInfinity);
  _beta = IloNumVarArray(env,Ng,0,1);

  _yup = IloRangeArray(env, Nl,0,IloInfinity);
  _ydown = IloRangeArray(env,Nl,0,IloInfinity);
  _lineprob = IloRangeArray(env,Nl,0,IloInfinity);
  _genup = IloRangeArray(env, Ng,0,IloInfinity);
  _gendown = IloRangeArray(env,Ng,0,IloInfinity);

  //  _riskConstraint.setExpr( IloSum(_z) );
  _betaSum = IloRange(env,1,1,"betasum");
  _betaSum.setExpr( IloSum(_beta) );

  for(int i=0;i<Nl;i++){
    ss.str("");
    ss<<"sd"<<i<<"[0,inf]";
    _sd[i].setName( ss.str().c_str() );
  }
  cout<<"GenInfo (epsG="<<_epsG<<")"<<endl;

  for(int j=0;j<Ng;j++){
    ss.str("");
    ss<<"beta"<<j<<"[0,1]";
    _beta[j].setName( ss.str().c_str() );

    double pmax = getGrid()->getGen(j).getPmax();
    //    double pmin = getGrid()->getGen(j).getPmin();
    double eta;
    if(_epsG == 1) eta=0;
    else eta = rv.PHIInverse(1-_epsG);

    cout<<"gen "<<j<<": "<<pmax<<endl;
    _genup[j].setExpr(getG()[j] - _beta[j]*sqrt(_sig_delta)*eta );
    _genup[j].setBounds(0,pmax);
    //    _gendown[j].setExpr( getG()[j] - _beta[j]*sqrt(_sig_delta)*eta );
    //    _gendown[j].setBounds(pmin,IloInfinity);

    cout<<"eta*sqrt(TV): "<<eta*sqrt(_sig_delta)<<endl;
    cout<<_genup[j]<<endl;
    getModel()->add(_genup[j]);
  }
  
  addCost(_beta,_sig_delta);

  getModel()->add(_sd);

  getModel()->add(_beta);
  //  getModel()->add(_genup);
  //  getModel()->add(_gendown);
  getModel()->add(_betaSum);


}

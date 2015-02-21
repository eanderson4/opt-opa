#include "isjn.h"

rgrid *  isjn::solveModel( isolve * is){

  rgrid * rg = new rgrid();
  
  time_t tstart;
  tstart = clock();

  int Nl = getGrid()->numBranches();
  int Ng = getGrid()->numGens();
  int Nm = getCm().n_cols;
  
  IloCplex cplex(*getModel());
  
  if(is!=NULL) is->setCplexParams(&cplex);
  int n=0;


  cout<<"Check: "<<sum(_check)<<" / "<<getGrid()->numBranches()<<endl;
  
  vec check2(Nl,fill::zeros);
  bool fullcheck=false;


  if (cplex.solve()){
    bool systemfail=true;
    mat Am = getA()*getCm();
    vec ones(Nm,1,fill::ones);
    vec onesL(Nl,1,fill::ones);
    mat ker = ones.t()*getSIGm()*ones;
    mat AmSone = Am*getSIGm()*ones;
    mat Ag = getA()*getCg();
    while(systemfail){
      n++; if (n>100) throw itlimit;  
      cout<<" --SJN-- --- Iteration - "<<n<<" ----- "<<endl;

      IloNumArray xsolve(getEnv(),Ng);
      IloNumArray betasolve(getEnv(),Ng);
      IloNumArray ysolve(getEnv(),Nl);

      cplex.getValues(xsolve,getG());
      cplex.getValues(ysolve,getF());
      cplex.getValues(betasolve,getBetaVar());

      vec x=getGC()->convert(xsolve);      
      vec y=getGC()->convert(ysolve);      
      vec beta=getGC()->convert(betasolve);      

      mat term3=Ag*beta;
      mat SIGy = term3*ker*term3.t()-term3*AmSone.t()-AmSone*term3.t()+getSigger();
      //      mat term = getA()*(getCg()*beta*ones.t() - getCm());    
      //      mat SIGy = term*getSIGm()*term.t();
      vec z=getGC()->risk(y,SIGy.diag(),getL(),getP(),getPc());


      x.t().print("x: ");
      y.t().print("y: ");
      beta.t().print("beta: ");
      //      pi.t().print("pi: ");
      z.t().print("z: ");
      cout<<"r: "<<accu(z)<<endl;

      cout<<"Base system"<<endl;
      systemfail = postCC(y,z,beta,SIGy.diag(),&cplex);
      
      cout<<getGenUp()[1]<<endl;
      cout<<x(1)<<"\t"<<beta(1)<<endl;

      //      vec ycheck = onesL*y + _L*y;
      

      for(int i=0;i<Nl;i++){
	if((_check(i) && check2(i)) || (_check(i) && fullcheck)){
	  cout<<"Contingency: "<<i<<endl;
	  vec yn = getN1(i,y);
	  //	  vec sdn = getSDN(i,y,SIGy);
	  vec sdn(Nl,fill::zeros);
  
	  for(int e=0;e<Nl;e++){
	    sdn(e) = SIGy(e,e) + 2*_L(e,i)*SIGy(e,i)+ _L(e,i)*_L(e,i)*SIGy(i,i);
	    if(sdn(e)<0.000001) sdn(e)=0;
	  }

	  /*	  double sdT=getSigDelta();
	  vec sig=getSig();

	  vec psi_en(Nl,fill::zeros);
	  for(int e=0;e<Nl;e++){
	    psi_en(e) = pi(e) + _L(e,i)*pi(i);
	  }
	  vec sd_en(Nl,fill::zeros);
	  for(int e=0;e<Nl;e++){
	    double term= psi_en(e)*psi_en(e)*sdT - 2*psi_en(e)*(sig(e) + _L(e,i)*sig(i)) + _sigpsi(e,i) ;
	    if(term<0 && term>=-.0000001) term=0;
	    sd_en(e) = sqrt( term);
	  }
	  vec sigpsi = _sigpsi.col(i);
	  vec error_sd = sdn - square(sd_en); 
	  */

	  vec zn=getGC()->risk(yn,sdn,getL(),getP(),getPc());
	  if(!yn.is_finite()) yn.print("yn: ");
	  cout<<"Post eval"<<endl;
	  check2(i) =postN1(i,yn,zn,beta,sdn,&cplex,n);
	  systemfail += check2(i);
	}
	else cout<<"dC "<<i<<"\t";
      }
      //      if(n<=5) systemfail=true;
      if(!systemfail && fullcheck){
	cout<<"Solved!"<<endl;
	setBetaSolve(beta);
	setSDSolve(SIGy.diag());
	break;
      }
      if(systemfail && fullcheck) fullcheck=false;
      if(!systemfail && !fullcheck) {
	fullcheck=true;
	systemfail=true;
      }


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
    
    cout<<"\n\nRisk constraint satisfied\n\n"<<endl;
    cout<<"Iterations: "<<n<<endl;
    cout<<"Time: "<<total<<endl;
    //    _addCut.print("cuts: ");        
  }

  else{
    cerr<<"Not solved"<<endl;
    cerr<<cplex.getStatus()<<endl;
  }
  cplex.end();
  
  return rg;
}


void isjn::setup(){
  cout<<"--sjn--  ===Setup=== N-1 joint chance constraint"<<endl;
  stringstream ss;
  
  IloEnv env = getEnv();

  int Nl = getGrid()->numBranches();
  gridcalc gc(getGrid());
  _addCut = mat(Nl,Nl,fill::zeros);
  _check = vec(Nl,fill::ones);

  _Hb=getA()*trans(getCb());
  cout<<"Build L"<<endl;
  _L=gc.getL(getA());

  mat sigger = getSigger();
  mat sigpsi(Nl,Nl);

  for(int e=0;e<Nl;e++){
    for(int n=0;n<Nl;n++){
      double term = 0;
      term = sigger(e,e)+2*_L(e,n)*sigger(e,n)+ _L(e,n)*_L(e,n)*sigger(n,n);
      sigpsi(e,n)=term;
    }
  }
 
  _sigpsi = sigpsi;

  for(int n=0;n<Nl;n++){  //Line i outage
    for(int j=0;j<Nl;j++){
      if(_Hb(n,n)>=1-.000001 && _Hb(n,n)<=1+.0000001){
	_check(n)=0;
      }
    }
  }

  cout<<sum(_check)<<" / "<<getGrid()->numBranches()<<endl;
  cout<<"there"<<endl;

  _in = mat(Nl,Nl,fill::zeros);

  _riskConstraint = IloRangeArray(env,Nl,0,IloInfinity);
  for(int n=0;n<Nl;n++){
    _z.push_back(IloNumVarArray(env,1,0,IloInfinity));
    _yplus.push_back(IloNumVarArray(env,1,0,IloInfinity));
    _sd.push_back(IloNumVarArray(env,1,0,IloInfinity));
    _yup.push_back(IloRangeArray(env, 1,0,IloInfinity));
    _ydown.push_back(IloRangeArray(env,1,0,IloInfinity));

  }
  getModel()->add(_riskConstraint);  
  cout<<"DONT BUILDING"<<endl;
  
}

bool isjn::postN1(int n, vec yn, vec zn, vec beta, vec sdn, IloCplex * cplex, int iteration){
  
  if(_check(n)!=1) return false;
  //define tolerance for line risk > 0
  stringstream ss;

  double tol = 1*pow(10,-3);
  int Nl = getGrid()->numBranches();

  mat A=getA();
  vec indexG = getIndexG();

  ranvar rv;
  double L = getL();
  double p = getP();
  double pc = getPc();
  double Ueps = rv.ginv(_epsN(n),L,p,pc);
  double rn = sum(zn);
  cout<<"Risk: "<<rn<<endl;
  if(rn<=_epsN(n)+tol) return false;  // SYSTEM SUCCEEDED
  else{
    cout<<"CUTTING ----------"<<endl;
    for(int i=0; i<Nl; i++){
      if (zn(i)>tol){
	if(sum(_addCut.col(n))==0){
	  _z[n] = IloNumVarArray(getEnv(),Nl,0,_epsN(n));	  
	  _yplus[n] = IloNumVarArray(getEnv(),Nl,0,IloInfinity);    
	  _sd[n] = IloNumVarArray(getEnv(),Nl,0,IloInfinity);	  
	  _yup[n] = IloRangeArray(getEnv(),Nl,0,IloInfinity);	  
	  _ydown[n] = IloRangeArray(getEnv(),Nl,0,IloInfinity);	  
	  _riskConstraint[n].setBounds(0,_epsN(n));
	}
	if(_addCut(i,n)==0){
	  _riskConstraint[n].setExpr( _riskConstraint[n].getExpr() + _z[n][i] );
	  _yup[n][i].setExpr( _yplus[n][i] - (getF()[i] + _L(i,n)*(getF()[n])) );
	  _ydown[n][i].setExpr( _yplus[n][i] + (getF()[i] + _L(i,n)*getF()[n]) );
	  double U = getGrid()->getBranch(i).getRateA();
	  _yplus[n][i].setBounds(0,U*Ueps);	  
	  getModel()->add(_yup[n][i]);
	  getModel()->add(_ydown[n][i]);
	  ss.str("");
	  ss<<"yplus"<<n<<","<<i<<"[0,"<<U*Ueps<<"]";
	  _yplus[n][i].setName( ss.str().c_str() );
	  ss.str("");
	  ss<<"sd"<<n<<","<<i<<"[0,inf]";
	  _sd[n][i].setName( ss.str().c_str() );
	  ss.str("");
	  ss<<"z"<<n<<","<<i<<"[0,"<<_epsN(n)<<"]";
	  _z[n][i].setName( ss.str().c_str() );
	  
	}
	
	//Add cuts for each line with positive risk
	double y_i = abs(yn(i));
	double sd_i = sqrt(sdn(i));
	double U=getGrid()->getBranch(i).getRateA();
	double dmu=rv.deriveMu(L,p,pc,y_i/U,sd_i/U);
	double dsigma=rv.deriveSigma(L,p,pc,y_i/U,sd_i/U);
	cout<<"dmu: "<<dmu<<", dsigma: "<<dsigma<<endl;
	IloRange cut(getEnv(),-IloInfinity,0);
	cut.setExpr( dmu/U*(_yplus[n][i] - y_i) + dsigma/U*(_sd[n][i] - sd_i) + zn(i) - _z[n][i]);
	cout<<cut<<endl;
	getModel()->add(cut);
	_addCut(i,n)=_addCut(i,n)+1;
	
	//Add cuts to describe standard deviation of branch flow
	IloRange cut_sd(getEnv(),-IloInfinity,0);

	int Ng = getGrid()->numGens();
	double pi_i = dot(getA().row(i),getCg()*beta);
	double pi_n = dot(getA().row(n),getCg()*beta);
	double psi = pi_i + _L(i,n)*pi_n;

	cut_sd.setExpr( sd_i - _sd[n][i] );
	double term=0;
	//	if(sd_i>0)
	  term=(psi * getSigDelta() - (getSig()(i) + _L(i,n)*getSig()(n)))/sd_i;
	cout<<"\n";
	for( int j=0; j<Ng;j++){
	  double pf_j = (A(i,indexG(j)) + _L(i,n)*A(n,indexG(j)))*term;
	  cut_sd.setExpr( cut_sd.getExpr() + pf_j*(getBetaVar()[j]-beta(j)) );
	}
	cout<<"\n";
	cout<<"sd_i: "<<sd_i<<" - "<<sqrt(sdn(i))<<endl;
	cout<<"pi_i: "<<pi_i<<endl;

	cout<<cut_sd<<endl;
	getModel()->add(cut_sd);

	cout<<"RC: ";
	cout<<_riskConstraint[n]<<endl;
	cout<<"\n\n";
       	    
      }
    }
    return true;
  }     
  return true;
}


vec isjn::getN1(int n, vec y0){
  return  (y0[n]*_L.col(n)+y0);
}

vec isjn::getSDN(int n, vec y0, mat Cov){
  int Nl = getGrid()->numBranches();
  vec sd(Nl,fill::zeros);
  
  for(int e=0;e<Nl;e++){
    sd(e) = Cov(e,e) + 2*_L(e,n)*Cov(e,n)+ _L(e,n)*_L(e,n)*Cov(n,n);
  }

  return sd;
}


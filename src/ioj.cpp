#include "ioj.h"

rgrid *  ioj::solveModel( isolve * is){

  rgrid * rg = new rgrid();
  
  time_t tstart;
  tstart = clock();

  int Nl = getGrid()->numBranches();
  int Ng = getGrid()->numGens();
  int Nm = getCm().n_cols;
  
  IloCplex cplex(*getModel());
  
  if(is!=NULL) is->setCplexParams(&cplex);
  int n=0;

  //  _xdes.t().print("design weighting: ");
  cout<<"Design weighting (first 10)"<<endl;
  for(int ci=0;ci<10;ci++) cout<<_xdes(ci)<<"\t";
  cout<<endl;
  cout<<"Check: "<<sum(getCheck())<<" / "<<getGrid()->numBranches()<<endl;
  
  vec check2(Nl,fill::zeros);
  bool fullcheck=false;


  if (cplex.solve()){
    bool systemfail=true;
    int numfailed=0;
    mat Am = getA()*getCm();
    vec ones(Nm,1,fill::ones);
    vec onesL(Nl,1,fill::ones);
    vec zprev(Nl,1,fill::zeros);
    mat ker = ones.t()*getSIGm()*ones;
    mat AmSone = Am*getSIGm()*ones;
    mat Ag = getA()*getCg();
    mat Lmat = getLmat();
    while(systemfail){
      n++; if (n>100) { cout<<"IT LIMIT!!"<<endl; throw itlimit;  }
      cout<<endl;
      cout<<" --OJ-- --- Iteration - "<<n<<" ----- "<<endl;

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
      vec l=_xdes % z;

      //      x.t().print("x: ");
      //      y.t().print("y: ");
      //      beta.t().print("beta: ");
      //      pi.t().print("pi: ");
      //      z.t().print("z: ");
      //      l.t().print("l: ");
      cout<<"r: "<<accu(z)<<endl;
      cout<<"l_tot: "<<accu(l)<<endl;


      cout<<"Base system"<<endl;
      systemfail = postLS(y,z,zprev,l,beta,SIGy.diag(),&cplex);
      zprev=z;
      //cout<<getGenUp()[1]<<endl;
      //      cout<<x(1)<<"\t"<<beta(1)<<endl;

      //      vec ycheck = onesL*y + Lmat*y;
      
      numfailed=0;
      for(int i=0;i<Nl;i++){
	if((getCheck()(i) && check2(i)) || (getCheck()(i) && fullcheck)){
	  cout<<"Contingency: "<<i<<endl;
	  vec yn = getN1(i,y);
	  //	  vec sdn = getSDN(i,y,SIGy);
	  vec sdn(Nl,fill::zeros);
  
	  for(int e=0;e<Nl;e++){
	    sdn(e) = SIGy(e,e) + 2*Lmat(e,i)*SIGy(e,i)+ Lmat(e,i)*Lmat(e,i)*SIGy(i,i);
	    if(sdn(e)<0.0000001) sdn(e)=0;
	  }

	  /*	  double sdT=getSigDelta();
	  vec sig=getSig();

	  vec psi_en(Nl,fill::zeros);
	  for(int e=0;e<Nl;e++){
	    psi_en(e) = pi(e) + Lmat(e,i)*pi(i);
	  }
	  vec sd_en(Nl,fill::zeros);
	  for(int e=0;e<Nl;e++){
	    double term= psi_en(e)*psi_en(e)*sdT - 2*psi_en(e)*(sig(e) + Lmat(e,i)*sig(i)) + _sigpsi(e,i) ;
	    if(term<0 && term>=-.0000001) term=0;
	    sd_en(e) = sqrt( term);
	  }
	  vec sigpsi = _sigpsi.col(i);
	  vec error_sd = sdn - square(sd_en); 
	  */

	  vec zn=getGC()->risk(yn,sdn,getL(),getP(),getPc());
	  vec ln=_xdes % zn;
	  //	  if(!yn.is_finite()) yn.print("yn: ");
	  cout<<"Post eval"<<endl;
	  check2(i) =postLSN1(i,yn,zn,ln,beta,sdn,&cplex,n);
	  systemfail += check2(i);
	  numfailed += check2(i);
	  if(numfailed>3*n) {
	    cout<<"\n\nOver "<<3*n<<" contingencies have failed, resolve\n"<<endl;
	    fullcheck=false;
	  }
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


    cout<<"RC: "<<_riskConstraintLS<<endl;
    cout<<"_l0: "<<_l0<<endl;
    cout<<"z: "<<getZ()<<endl;
    
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


void ioj::setup(){
  cout<<"--oj--  ===Setup=== N-1 joint chance constraint with LOADSHED"<<endl;
  stringstream ss;
  
  IloEnv env = getEnv();

  int Nl = getGrid()->numBranches();
  gridcalc gc(getGrid());

  cout<<"HeRE"<<endl;

  _addCut = vec(Nl,fill::zeros);
  _addCutN1 = mat(Nl,Nl,fill::zeros);
  _in = mat(Nl,Nl,fill::zeros);

  _riskConstraintLS = IloRange(env,0,_epsLS);
  _l0 = IloNumVarArray(env,Nl,0,_epsLS);
  _l0eq = IloRangeArray(env,Nl,0,IloInfinity);
  _riskConstraintLSN = IloRangeArray(env,Nl,0,IloInfinity);
  for(int n=0;n<Nl;n++){
    _lN.push_back(IloNumVarArray(env,1,0,IloInfinity));
    _lNeq.push_back(IloRangeArray(env,1,0,IloInfinity));
        _z.push_back(IloNumVarArray(env,1,0,IloInfinity));
        _yplus.push_back(IloNumVarArray(env,1,0,IloInfinity));
        _sd.push_back(IloNumVarArray(env,1,0,IloInfinity));
        _yup.push_back(IloRangeArray(env, 1,0,IloInfinity));
        _ydown.push_back(IloRangeArray(env,1,0,IloInfinity));
    getF()[n].setBounds(-IloInfinity,IloInfinity);
    getZ()[n].setBounds(0,IloInfinity);
  }
  getModel()->add(_riskConstraintLS);  
  getModel()->add(_riskConstraintLSN);  
  //  getModel()->add(_riskConstraintLSN);  
  

  getRiskConstraint().setBounds(0,IloInfinity);
  //  getRiskConstraintN1().setBounds(0,IloInfinity);

  cout<<"DONT BUILDING"<<endl;
  
}

bool ioj::postLS(vec y, vec z,vec zprev,vec l, vec beta, vec SIGy,IloCplex * cplex, int iteration){
  stringstream ss;
  //define tolerance for line risk > 0
  double tol = 5*pow(10,-4);
  int Nl = getGrid()->numBranches();

  double account=0;
  double accountLS=0;
  
  mat A=getA();  
  ranvar rv;
  double eps = getEps();
  double L = getL();
  double p = getP();
  double pc = getPc();
  double Ueps = rv.ginv(eps,L,p,pc);
  double r = sum(z);
  double rl = sum(l);
    cout<<"Risk: "<<r<<endl;
    cout<<"Risk LS: "<<rl<<", Desired: "<<_epsLS<<endl;
  if(rl<=_epsLS+tol) return false;
  else{
    cout<<"CUTTING ----------"<<endl;
    vec pi = getA()*getCg()*beta;
    for(int i=0; i<Nl; i++){
      //      cout<<i<<endl;
      if (z(i)>.0000001){
      if (zprev(i) != z(i)){
       	account += z(i);
       	accountLS += l(i);
	cout<<"Line "<<i<<" cuts"<<endl;
	cout<<"z: "<<z(i)<<", l: "<<l(i)<<endl;
	if(_addCut(i)==0){
	  cout<<"Initialize Cutting Variables for Line "<<i<<endl;
	  _riskConstraintLS.setExpr( _riskConstraintLS.getExpr() + _l0[i]);
	  getYUp()[i].setExpr( getYplus()[i] - getF()[i] );
	  getYDown()[i].setExpr( getYplus()[i] + getF()[i] );
	  //	  cout<<"Set up N-1 Flow Variables \t"<<endl;
	  double U = getGrid()->getBranch(i).getRateA();
	  //	  cout<<"Get U\t";
	  //	  getYplus[i].setBounds(0,U*Ueps);	  
	  getYplus()[i].setBounds(0,IloInfinity);	  
	  getModel()->add(getYUp()[i]);
	  getModel()->add(getYDown()[i]);
	  ss.str("");
	  //	  ss<<"yplus"<<i<<"[0,"<<U*Ueps<<"]";
	  ss<<"yplus"<<i<<"[0,inf]";
	  getYplus()[i].setName( ss.str().c_str() );
	  //	  cout<<"Add and name \t"<<endl;
	  _l0eq[i].setExpr( _l0[i] - _xdes(i)*getZ()[i] );
	  getModel()->add(_l0eq[i]);
	  ss.str("");
	  ss<<"l"<<i<<"[0,"<<_epsLS<<"]";
	  _l0[i].setName( ss.str().c_str() );
	  //	  cout<<"LS risk constraint"<<endl;

	}
	//Add cuts for each line with positive risk
	//	cout<<"Calculated "<<endl;
	double y_i = abs(y(i));
	//	cout<<"\ty"<<endl;
	double sd_i = sqrt(SIGy(i));
	//	cout<<"\tsd"<<endl;
	//	cout<<"z_"<<i<<": "<<z(i)<<", y_"<<i<<": "<<y_i<<endl;
	double U=getGrid()->getBranch(i).getRateA();
	//	cout<<"\tU"<<endl;
	double dmu=rv.deriveMu(getL(),getP(),getPc(),y_i/U,sd_i/U);
	//	cout<<"Got dmu"<<endl;
	double dsigma=rv.deriveSigma(getL(),getP(),getPc(),y_i/U,sd_i/U);
	//	cout<<"Got dsigma\t"<<endl;
	IloRange cut(getEnv(),-IloInfinity,0);
	cut.setExpr( dmu/U*(getYplus()[i] - y_i) + dsigma/U*(getSDVar()[i] - sd_i) + z(i) - getZ()[i]);
	//	cout<<cut<<endl;
	getModel()->add(cut);
	//	cout<<"Added Cut\t"<<endl;
	_addCut(i)=_addCut(i)+1;
	
	//Add cuts to describe standard deviation of branch flow
	IloRange cut_sd(getEnv(),-IloInfinity,0);

	int Ng = getGrid()->numGens();

	//	double pi_i = pi(i);


	cut_sd.setExpr( sd_i - getSDVar()[i] );
	double term = (pi(i)*getSigDelta() - getSig()(i))/sd_i;
	//	cout<<"\n";
	for( int j=0; j<Ng;j++){
	  //double pf_j = term;
	  //	  cout<<j<<",";
	  //	  cout<<"A: "<<getA()(i,getIndexG()(j));
	  //	  cout<<"in ("<<getGrid()->getGen(j).getPmin()<<",";
	  //	  cout<<getGrid()->getGen(j).getPmax()<<")"<<endl;
	  double pf_j = A(i,getIndexG()(j))*term;
	  cut_sd.setExpr( cut_sd.getExpr() + pf_j*(getBetaVar()[j]-beta(j)) );
	}
	
	
	//	cout<<cut_sd<<endl;
	getModel()->add(cut_sd);
	//	cout<<"Added SD Cut"<<endl;

	//	cout<<"\n\n";
      }
      else{ cout<<"z("<<i<<") has not changed, no cut"<<endl; }
      }
    }
    cout<<"Accounted: "<<account<<endl;
    cout<<"Accounted LS: "<<accountLS<<endl;
    return true;
  }     

  return true;
}



bool ioj::postLSN1(int n, vec yn, vec zn, vec l,vec beta, vec sdn, IloCplex * cplex, int iteration){
  
  if(getCheck()(n)!=1) return false;
  //define tolerance for line risk > 0
  stringstream ss;

  double tol = 1*pow(10,-3);
  int Nl = getGrid()->numBranches();

  mat A=getA();
  mat Lmat=getLmat();
  vec indexG = getIndexG();

  double account=0;
  double accountLS=0;

  ranvar rv;
  double L = getL();
  double p = getP();
  double pc = getPc();
  double Ueps = rv.ginv(getEpsN()(n),L,p,pc);
  double rn = sum(zn);
  double rln = sum(l);
  
    cout<<"Risk: "<<rn<<endl;
  cout<<"Risk LS: "<<rln<<endl;
  
  if(rln<=_epsLSN(n)+tol+.1) return false;  // SYSTEM SUCCEEDED
  else{
    cout<<"CUTTING ----------"<<endl;
    //    l.t().print("LS: ");
    for(int i=0; i<Nl; i++){
      if (zn(i)>tol/100){
       	account += zn(i);
       	accountLS += l(i);
	cout<<"line: "<<i<<" z:"<<zn(i)<<" l:"<<l(i)<<endl;
	if(sum(_addCutN1.col(n))==0){
	  cout<<"Initialize Cutting Variables for contingency "<<n<<endl;
	  _lN[n] = IloNumVarArray(getEnv(),Nl,0,_epsLSN(n));	  
	  _lNeq[n] = IloRangeArray(getEnv(),Nl,0,0);	  
	  _z[n] = IloNumVarArray(getEnv(),Nl,0,IloInfinity);	  
	  _yplus[n] = IloNumVarArray(getEnv(),Nl,0,IloInfinity);    
	  _sd[n] = IloNumVarArray(getEnv(),Nl,0,IloInfinity);	  
	  _yup[n] = IloRangeArray(getEnv(),Nl,0,IloInfinity);	  
	  _ydown[n] = IloRangeArray(getEnv(),Nl,0,IloInfinity);	  
	  _riskConstraintLSN[n].setBounds(0,_epsLSN(n));
	}
	if(_addCutN1(i,n)==0){
	  cout<<"Initialize Cutting Variables for Line "<<i<<" in contingency "<<n<<endl;
	  _lNeq[n][i].setExpr( _lN[n][i] - _xdes(i)*_z[n][i] );
	  _riskConstraintLSN[n].setExpr( _riskConstraintLSN[n].getExpr() + _lN[n][i] );
	  _yup[n][i].setExpr( _yplus[n][i] - (getF()[i] + Lmat(i,n)*(getF()[n])) );
	  _ydown[n][i].setExpr( _yplus[n][i] + (getF()[i] + Lmat(i,n)*getF()[n]) );
	  double U = getGrid()->getBranch(i).getRateA();
	  _yplus[n][i].setBounds(0,IloInfinity);	  
	  getModel()->add(_yup[n][i]);
	  getModel()->add(_ydown[n][i]);
	  getModel()->add(_lNeq[n][i]);
	  ss.str("");
	  ss<<"yplus"<<n<<","<<i<<"[0,inf]";
	  _yplus[n][i].setName( ss.str().c_str() );
	  ss.str("");
	  ss<<"sd"<<n<<","<<i<<"[0,inf]";
	  _sd[n][i].setName( ss.str().c_str() );
	  ss.str("");
	  ss<<"z"<<n<<","<<i<<"[0,inf]";
	  _z[n][i].setName( ss.str().c_str() );
	  ss.str("");
	  ss<<"lN"<<n<<","<<i<<"[0,"<<_epsLSN(n)<<"]";
	  _lN[n][i].setName( ss.str().c_str() );

	}
	
	//Add cuts for each line with positive risk
	double y_i = abs(yn(i));
	double sd_i = sqrt(sdn(i));
	double U=getGrid()->getBranch(i).getRateA();
	double dmu=rv.deriveMu(L,p,pc,y_i/U,sd_i/U);
	double dsigma=rv.deriveSigma(L,p,pc,y_i/U,sd_i/U);
	//	cout<<"dmu: "<<dmu<<", dsigma: "<<dsigma<<endl;
	IloRange cut(getEnv(),-IloInfinity,0);
	cut.setExpr( dmu/U*(_yplus[n][i] - y_i) + dsigma/U*(_sd[n][i] - sd_i) + zn(i) - _z[n][i]);
	//	cout<<cut<<endl;
	getModel()->add(cut);
	_addCutN1(i,n)=_addCutN1(i,n)+1;
	
	//Add cuts to describe standard deviation of branch flow
	IloRange cut_sd(getEnv(),-IloInfinity,0);

	int Ng = getGrid()->numGens();
	double pi_i = dot(getA().row(i),getCg()*beta);
	double pi_n = dot(getA().row(n),getCg()*beta);
	double psi = pi_i + Lmat(i,n)*pi_n;

	cut_sd.setExpr( sd_i - _sd[n][i] );
	double term=0;
	//	if(sd_i>0)
	  term=(psi * getSigDelta() - (getSig()(i) + Lmat(i,n)*getSig()(n)))/sd_i;
	  //	cout<<"\n";
	  if(sd_i>0){
	    for( int j=0; j<Ng;j++){
	      double pf_j = (A(i,indexG(j)) + Lmat(i,n)*A(n,indexG(j)))*term;
	      /*	  if (std::isinf(pf_j)){
			  cout<<"A(i,g): "<<A(i,indexG(j))<<endl;
			  cout<<"Lmat: "<<Lmat(i,n)<<endl;
			  cout<<"A(n,g): "<<A(n,indexG(j))<<endl;
			  cout<<"term: "<<term<<endl;
			  cout<<"psi: "<<psi<<endl;
			  cout<<"sd_i: "<<sd_i<<endl;
			  }*/
	      cut_sd.setExpr( cut_sd.getExpr() + pf_j*(getBetaVar()[j]-beta(j)) );
	    }
	    getModel()->add(cut_sd);
	  }
	    //	cout<<"\n";
	    //	cout<<"sd_i: "<<sd_i<<" - "<<sqrt(sdn(i))<<endl;
	    //	cout<<"pi_i: "<<pi_i<<endl;
	    
	//	cout<<"cut: "<<cut<<endl;
	//	cout<<"cut_sd: "<<cut_sd<<endl;


	//	cout<<"RC: ";
	//	cout<<_riskConstraintLSN[n]<<endl;
	//	cout<<"Lneq: "<<_lNeq[n][i]<<endl;
	//	cout<<"\n\n";
       	    
      }
    }

    cout<<"Accounted: "<<account<<endl;
    cout<<"Accounted LS: "<<accountLS<<endl;
    return true;
  }     


  return true;
}


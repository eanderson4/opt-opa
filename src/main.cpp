#include <stdlib.h>

#include <ctime>
#include "sqlinter.h"
#include "jsoninter.h"

#include "icc.h"
#include "isj.h"
#include "isjn.h"
#include "ijn1.h"
#include "in1.h"
#include "iopa.h"

using namespace std;

int main(int argc, char* argv[]){

  if(argc<=12){
    cout<<"cmd: pow case/30.db <m0> <m1> <e0> <e1> <L> <p> <B> <T> <L2> <p2>\n"
	<<"\trun main for case30\n"
	<<"\t<m0> base capacity multiplier\n"
	<<"\t<m1> contingency capacity multiplier\n"
	<<"\t<e0> base risk\n"
	<<"\t<e1> contingency risk\n"
	<<"\t<pl> line fail prob\n"
	<<"\t<pg> generator fail prob\n"
	<<"\t<L> no risk intercept\n"
	<<"\t<p> probability of failure at nominal\n"
	<<"\t<B> Variance budget\n"
	<<"\t<T> num trials\n"
	<<"\t<L2> no risk intercept (OPA)\n"
	<<"\t<p2> probability of failure at nominal (OPA)\n"<<endl;
    return 1;
  }

  double m0=atof(argv[2]);
  double m1=atof(argv[3]);
  double mg=atof(argv[4]);
  double eps=atof(argv[5]);
  double epsN=atof(argv[6]);
  double epsL=atof(argv[7]);
  double epsG=atof(argv[8]);
  double L=atof(argv[9]);
  double p=atof(argv[10]);
  double B=atof(argv[11]);
  double pc=.85;
  int T = atoi(argv[12]);

  double opaL,opap;
  if(argc>13){
    opaL=atof(argv[13]);
    opap=atof(argv[14]);
  }
  else {
    opaL=L;
    opap=.5;
  }
   

  ///  int sn=atoi(argv[7]); //standard deviation test

  sqlInter db;
  grid * gr = new grid;

  string db_name;
  
  db_name = argv[1];
  db.openDb(db_name);
  db.load(*gr);
  gr->buildMap();
  gr->printNums(cout);

  cout<<*gr<<endl;
  


  int Nl = gr->numBranches();
  for(int i=0;i<Nl;i++){
    branch bi = gr->getBranch(i);
    double U = bi.getRateA();
    gr->setCapacity(i,m0*U);
  }

  int Ng = gr->numGens();
  for(int j=0;j<Ng;j++){
    gen gj = gr->getGen(j);
    double U = gj.getPmax();
    gr->setGenPmax(j,mg*U);
  } 

  int Nb = gr->numBuses();
  vector<int> nodes;
  vector<int> demand;
  for(int i=0;i<Nb;i++){
    bus bn = gr->getBus(i);
    double p = bn.getP();
    
    if(p>1){ 
      nodes.push_back(i);
      demand.push_back(p);
    }  
  }

  int Nm=nodes.size();
  vec indexM(Nm,fill::zeros);


  mat Cm(Nb,Nm,fill::zeros);
  vec mu(Nm,fill::zeros);
  mat SIG(Nm,Nm,fill::zeros);
  for(int i=0;i<Nm;i++){
    Cm(nodes[i],i)=1;
    indexM(i)=nodes[i];
    SIG(i,i)=pow(.05*demand[i],2)*B;
  }
  SIG.diag().t().print("StDv Volatile Injects: ");
  double TV = accu(SIG);

  //Generate second demand scenario
  arma_rng::set_seed_random();
  mat as = chol(SIG);

  cout<<*gr<<endl;

  arma_rng::set_seed_random(); 
  gridcalc gc(gr);  
    

  mat Cg=gc.getCm();
  vec alpha(Ng,fill::zeros);

  for(int i=0;i<Ng;i++){
    gen gi=gr->getGen(i);
    double status = gi.getStatus();
    if (status>0)   alpha(i)=1;
  }
  int numSlack = accu(alpha);
  alpha = alpha/numSlack;
  alpha.t().print("alpha: ");

  vec slack=Cg*alpha;

  double randcost=0;
  for(int i=0;i<Ng;i++){
    double c2=gr->getGen(i).getC2();
    double genrandcost=c2*pow(alpha(i),2)*TV;
    randcost=randcost+genrandcost;
  }

  mat Hw = gc.getHw(slack);
  mat A = gc.getH();
  mat Lo = gc.getL(A);  


  mat SIGy(Nl,Nl,fill::zeros);
  SIGy = Hw*Cm*SIG*(Hw*Cm).t();
  vec sd = SIGy.diag();
  


  vec eN(Nl, fill::ones);
  eN = eN*epsN;
  

  //Set Probability info
  ranvar rv;


  igrid ig(gr);
  ig.addCost();
  isolve is;
  is.setSolver(IloCplex::Dual,IloCplex::Dual);
  //  rgrid * rbase = ig.solveModel(&is);

  double TG;
  
  vec ones(Nm,1,fill::ones);  
  
  vec ztemp(Nl,fill::zeros);
  //    iopa opa(gr,&gc,opaL,opap,L,p);
  //    opa.runTrials(cout,ztemp,15,10);
    
  //    double test;
  //    cin>>test;
  

    ijn1 n1(gr,  SIGy,Hw,L,p,pc,eps,epsN);
  n1.addCost();
  vec check = n1.getCheck();

    isj nom(gr, &gc, SIG, indexM, L, p, pc, 1);
    nom.lineLimitStatus(true);



      try{
	//	igrid nom(gr);
	//	nom.addCost();
	//	in1 nom1(gr, SIGy, Hw,m1);
	rgrid * rnom = nom.solveModel(&is);

	

	
	double o0=rnom->getObjective();
	vec f0=gc.convert(rnom->getF());
	vec g0=gc.convert(rnom->getG());
	vec beta0=nom.getBeta();
	mat term0 = A*(Cg*beta0*ones.t() - Cm);    
	mat SIGy0 = term0*SIG*term0.t();
	vec sd0 = SIGy0.diag();
	vec z0=gc.risk(f0,sd0,L,p,pc);
 	double r0 = sum(z0);
	vec p0=gc.lineprob(f0,SIGy.diag());
	IloCplex::CplexStatus s0=rnom->getStatus();
	TG = accu(g0);


	vec fU0(Nl);
	vec sdU0(Nl);
	for(int i=0;i<Nl;i++){
	  double U = gr->getBranch(i).getRateA();
	  fU0(i) = abs(f0(i))/U;
	  sdU0(i) = sd0(i)/U;
	}
	
	running_stat<double> stats_r0;
	for(int i=0;i<Nl;i++){
	  if(check(i)==1){
	    vec f0n = n1.getN1(i,f0,g0);
	    vec sdn(Nl);
	    for(int e=0;e<Nl;e++){
	      double U = gr->getBranch(i).getRateA();
	      sdn(e) = SIGy0(e,e) + 2*Lo(e,i)*SIGy0(e,i)+ Lo(e,i)*Lo(e,i)*SIGy0(i,i);
	      if(sdn(e)<0 && sdn(e)>=-.0000001) sdn(e)=0;
	      //	      else sdn(e)=sdn(e)/U;
	    }
	    vec z0n=gc.risk(f0n,sdn,L,p,pc);
	    double r0n = sum(z0n);
	    stats_r0(r0n);
	  }
	}
	
	cout<<"OPF"<<"\t"<<s0<<endl;
	cout<<"C0: "<<o0<<endl;
	cout<<"r0 - "<<r0<<endl;
	cout << "count = " << stats_r0.count() << endl;
	cout << "mean = " << stats_r0.mean() << endl;
	cout << "stdv  = " << stats_r0.stddev()  << endl;
	cout << "min  = " << stats_r0.min()  << endl;
	cout << "max  = " << stats_r0.max()  << endl;
	cout<<endl;
	cerr<<"opf\t"<<o0<<"\t"<<r0<<"\t"<<stats_r0.mean()<<"\t"<<stats_r0.max()<<endl;

	in1 nom1(gr, SIGy, Hw,m1);
	nom1.addCost();
	rgrid * rnom1 = nom1.solveModel(&is);
	
	double o1=rnom1->getObjective();
	vec f1=gc.convert(rnom1->getF());
	vec g1=gc.convert(rnom1->getG());
	vec beta1=nom.getBeta();
	mat term1 = A*(Cg*beta1*ones.t() - Cm);    
	mat SIGy1 = term1*SIG*term1.t();
	vec sd1 = SIGy1.diag();
	vec z1=gc.risk(f1,sd1,L,p,pc);
	double r1 = sum(z1);
	vec pnom1=gc.lineprob(f1,SIGy.diag());
	IloCplex::CplexStatus s1=rnom1->getStatus();
	TG = accu(g1);


	vec fU1(Nl);
	vec sdU1(Nl);
	for(int i=1;i<Nl;i++){
	  double U = gr->getBranch(i).getRateA();
	  fU1(i) = abs(f1(i))/U;
	  sdU1(i) = sd1(i)/U;
	}
	
	running_stat<double> stats_r1;
	for(int i=0;i<Nl;i++){
	  if(check(i)==1){
	    vec f1n = n1.getN1(i,f1,g1);
	    vec sdn(Nl);
	    for(int e=0;e<Nl;e++){
	      double U = gr->getBranch(i).getRateA();
	      sdn(e) = SIGy1(e,e) + 2*Lo(e,i)*SIGy1(e,i)+ Lo(e,i)*Lo(e,i)*SIGy1(i,i);
	      if(sdn(e)<0 && sdn(e)>=-.0000001) sdn(e)=0;
	      //	      else sdn(e)=sdn(e)/U;
	    }
	    vec z1n=gc.risk(f1n,sdn,L,p,pc);
	    double r1n = sum(z1n);
	    stats_r1(r1n);
	  }
	}
	
	cout<<"OPF n1"<<"\t"<<s1<<endl;
	cout<<"C1: "<<o1<<endl;
	cout<<"r1 - "<<r1<<endl;
	cout << "count = " << stats_r1.count() << endl;
	cout << "mean = " << stats_r1.mean() << endl;
	cout << "stdv  = " << stats_r1.stddev()  << endl;
	cout << "min  = " << stats_r1.min()  << endl;
	cout << "max  = " << stats_r1.max()  << endl;
	cout<<endl;
	cerr<<"opfn1\t"<<o1<<"\t"<<r1<<"\t"<<stats_r1.mean()<<"\t"<<stats_r1.max()<<endl;	


	

	icc cc(gr, &gc, SIG, indexM, epsL, epsG);
  	rgrid * rcc = cc.solveModel(&is);
	
	double o=rcc->getObjective();
	vec f=gc.convert(rcc->getF());
	vec g=gc.convert(rcc->getG());
	vec beta=cc.getBeta();
	mat term = A*(Cg*beta*ones.t() - Cm);    
	mat SIGycc = term*SIG*term.t();
	vec sd = SIGycc.diag();
	vec z=gc.risk(f,sd,L,p,pc);
	double r = sum(z);
	vec p1=gc.lineprob(f,sd);
	IloCplex::CplexStatus s=rcc->getStatus();


	vec fU(Nl);
	vec sdU(Nl);
	for(int i=0;i<Nl;i++){
	  double U = gr->getBranch(i).getRateA();
	  fU(i) = abs(f(i))/U;
	  sdU(i) = sd(i)/U;
	}
	
	running_stat<double> stats_risk;
	for(int i=0;i<Nl;i++){
	  if(check(i)==1){
	    vec fn = n1.getN1(i,f,g);
	    vec sdn(Nl);
	    for(int e=0;e<Nl;e++){
	      double U = gr->getBranch(i).getRateA();
	      sdn(e) = SIGycc(e,e) + 2*Lo(e,i)*SIGycc(e,i)+ Lo(e,i)*Lo(e,i)*SIGycc(i,i);
	      if(sdn(e)<0 && sdn(e)>=-.0000001) sdn(e)=0;
	      //	      else sdn(e)=sdn(e)/U;
	    }
	    vec zn=gc.risk(fn,sdn,L,p,pc);
	    double rn = sum(zn);
	    stats_risk(rn);
	  }
	}
	
	cout<<"CC"<<"\t"<<s<<endl;
	cout<<"C: "<<o<<endl;
	cout<<"r - "<<r<<endl;
	cout << "count = " << stats_risk.count() << endl;
	cout << "mean = " << stats_risk.mean() << endl;
	cout << "stdv  = " << stats_risk.stddev()  << endl;
	cout << "min  = " << stats_risk.min()  << endl;
	cout << "max  = " << stats_risk.max()  << endl;
	cout<<endl;
	cerr<<"cc\t"<<o<<"\t"<<r<<"\t"<<stats_risk.mean()<<"\t"<<stats_risk.max()<<endl;	

	

	isj sj(gr, &gc, SIG, indexM, L, p, pc, eps);
	rgrid * rsj = sj.solveModel(&is);

	
	double o3=rsj->getObjective();
	vec f3=gc.convert(rsj->getF());
	vec g3=gc.convert(rsj->getG());
	vec beta3=sj.getBeta();
	mat term3 = A*(Cg*beta3*ones.t() - Cm);    
	mat SIGy3 = term3*SIG*term3.t();
	vec sd3 = SIGy3.diag();
	vec z3=gc.risk(f3,sd3,L,p,pc);
	double r3 = sum(z3);
	vec p3=gc.lineprob(f3,sd3);
	IloCplex::CplexStatus s3=rsj->getStatus();



	vec fU3(Nl);
	vec sdU3(Nl);
	for(int i=0;i<Nl;i++){
	  double U = gr->getBranch(i).getRateA();
	  fU3(i) = abs(f3(i))/U;
	  sdU3(i) = sd3(i)/U;
	}
	
	running_stat<double> stats_r3;
	for(int i=0;i<Nl;i++){
	  if(check(i)==1){
	    vec f3n = n1.getN1(i,f3,g3);
	    vec sdn(Nl);
	    for(int e=0;e<Nl;e++){
	      double U = gr->getBranch(i).getRateA();
	      sdn(e) = SIGy3(e,e) + 2*Lo(e,i)*SIGy3(e,i)+ Lo(e,i)*Lo(e,i)*SIGy3(i,i);
	      if(sdn(e)<0 && sdn(e)>=-.0000001) sdn(e)=0;
	      //	      else sdn(e)=sdn(e)/U;
	    }
	    vec z3n=gc.risk(f3n,sdn,L,p,pc);
	    double r3n = sum(z3n);
	    stats_r3(r3n);
	  }
	}
	
	cout<<"SJ"<<"\t"<<s3<<endl;
	cout<<"C3: "<<o3<<endl;
	cout<<"r3 - "<<r3<<endl;
	cout << "count = " << stats_r3.count() << endl;
	cout << "mean = " << stats_r3.mean() << endl;
	cout << "stdv  = " << stats_r3.stddev()  << endl;
	cout << "min  = " << stats_r3.min()  << endl;
	cout << "max  = " << stats_r3.max()  << endl;
	cout<<endl;
	cerr<<"jcc\t"<<o3<<"\t"<<r3<<"\t"<<stats_r3.mean()<<"\t"<<stats_r3.max()<<endl;

	isjn sjn(gr, &gc, SIG, indexM, L, p, pc, eps, eN, .5);
	rgrid * rsjn = sjn.solveModel(&is);

	cout<<"HERE"<<endl;
	
	double o4=rsjn->getObjective();
	vec f4=gc.convert(rsjn->getF());
	vec g4=gc.convert(rsjn->getG());
	vec beta4=sjn.getBeta();
	mat term4 = A*(Cg*beta4*ones.t() - Cm);    
	mat SIGy4 = term4*SIG*term4.t();
	vec sd4 = SIGy4.diag();
	vec z4=gc.risk(f4,sd4,L,p,pc);
	double r4 = sum(z4);
	vec p4=gc.lineprob(f4,sd4);
	IloCplex::CplexStatus s4=rsjn->getStatus();

	vec fU4(Nl);
	vec sdU4(Nl);
	for(int i=0;i<Nl;i++){
	  double U = gr->getBranch(i).getRateA();
	  fU4(i) = abs(f4(i))/U;
	  sdU4(i) = sd4(i)/U;
	}
	
	running_stat<double> stats_r4;
	for(int i=0;i<Nl;i++){
	  if(check(i)==1){
	    vec f4n = n1.getN1(i,f4,g4);
	    vec sdn(Nl);
	    for(int e=0;e<Nl;e++){
	      double U = gr->getBranch(i).getRateA();
	      sdn(e) = SIGy4(e,e) + 2*Lo(e,i)*SIGy4(e,i)+ Lo(e,i)*Lo(e,i)*SIGy4(i,i);
	      if(sdn(e)<0 && sdn(e)>=-.0000001) sdn(e)=0;
	      //	      else sdn(e)=sdn(e)/U;
	    }
	    vec z4n=gc.risk(f4n,sdn,L,p,pc);
	    double r4n = sum(z4n);
	    stats_r4(r4n);
	  }
	}
	
	cout<<"SJ N1"<<"\t"<<s4<<endl;
	cout<<"C4: "<<o4<<endl;
	cout<<"r4 - "<<r4<<endl;
	cout << "count = " << stats_r4.count() << endl;
	cout << "mean = " << stats_r4.mean() << endl;
	cout << "stdv  = " << stats_r4.stddev()  << endl;
	cout << "min  = " << stats_r4.min()  << endl;
	cout << "max  = " << stats_r4.max()  << endl;
	cout<<endl;
	cerr<<"jccn4\t"<<o4<<"\t"<<r4<<"\t"<<stats_r4.mean()<<"\t"<<stats_r4.max()<<endl;	

    

    
    
    cout<<"2nd stage Comparison"<<endl;
    //    cout<<*gr<<endl;
    SIG.diag().t().print("Cov.diag m: ");
    alpha.t().print("slack: ");
    

    cout<<"Parameters: "<<endl;
    cout<<"m0\tm1\tmg\teps\tepsN\tpL\tpG\tL\tp\tB\tT"<<endl;
    cout<<m0<<"\t"<<m1<<"\t"<<mg<<"\t"<<eps<<"\t"<<epsN<<"\t"<<epsL<<"\t"<<epsG<<"\t"<<L<<"\t"<<p<<"\t"<<B<<"\t"<<T<<endl;

    cout<<"Total Gen: "<<TG<<endl;
    cout<<"Total Variance: "<<TV<<" ( "<<sqrt(TV)<<" )"<<endl;
    cout<<"Total Random Cost: "<<randcost<<endl;
    cout<<"U eps: "<<rv.ginv(eps,L,p,pc)<<endl;

    cout<<"\n\n";

    cout<<" --- \t --- RISK 2nd--- \t ---\n"<<endl;

	cout<<"OPF"<<"\t"<<s0<<endl;
	cout<<"C0: "<<o0<<endl;
	cout<<"r0 - "<<r0<<endl;
	cout << "count = " << stats_r0.count() << endl;
	cout << "mean = " << stats_r0.mean() << endl;
	cout << "stdv  = " << stats_r0.stddev()  << endl;
	cout << "min  = " << stats_r0.min()  << endl;
	cout << "max  = " << stats_r0.max()  << endl;
	cout<<endl;



	cout<<"OPF n1"<<"\t"<<s1<<endl;
	cout<<"C1: "<<o1<<endl;
	cout<<"r1 - "<<r1<<endl;
	cout << "count = " << stats_r1.count() << endl;
	cout << "mean = " << stats_r1.mean() << endl;
	cout << "stdv  = " << stats_r1.stddev()  << endl;
	cout << "min  = " << stats_r1.min()  << endl;
	cout << "max  = " << stats_r1.max()  << endl;
	cout<<endl;



	cout<<"CC"<<"\t"<<s<<endl;
	cout<<"C: "<<o<<endl;
	cout<<"r - "<<r<<endl;
	cout << "count = " << stats_risk.count() << endl;
	cout << "mean = " << stats_risk.mean() << endl;
	cout << "stdv  = " << stats_risk.stddev()  << endl;
	cout << "min  = " << stats_risk.min()  << endl;
	cout << "max  = " << stats_risk.max()  << endl;
	cout<<endl;


	cout<<"SJ"<<"\t"<<s3<<endl;
	cout<<"C3: "<<o3<<endl;
	cout<<"r3 - "<<r3<<endl;
	cout << "count = " << stats_r3.count() << endl;
	cout << "mean = " << stats_r3.mean() << endl;
	cout << "stdv  = " << stats_r3.stddev()  << endl;
	cout << "min  = " << stats_r3.min()  << endl;
	cout << "max  = " << stats_r3.max()  << endl;
	cout<<endl;


	cout<<"SJ N1"<<"\t"<<s4<<endl;
	cout<<"C4: "<<o4<<endl;
	cout<<"r4 - "<<r4<<endl;
	cout << "count = " << stats_r4.count() << endl;
	cout << "mean = " << stats_r4.mean() << endl;
	cout << "stdv  = " << stats_r4.stddev()  << endl;
	cout << "min  = " << stats_r4.min()  << endl;
	cout << "max  = " << stats_r4.max()  << endl;
	cout<<endl;

	cerr<<"\n";

	cout<<"\n\n";


    cout<<" --- \t --- RISK --- \t ---\n"<<endl;
    cout<<"OPF"<<endl;
    cout << "risk: "<<endl;
    cout << r0 <<endl;
    cout<<endl;

    cout<<"OPF n1"<<endl;
    cout << "risk: "<<endl;
    cout << r1 <<endl;
    cout<<endl;

    cout<<"CC"<<endl;
    cout << "risk: "<<endl;
    cout << r <<endl;
    cout<<endl;


    cout<<"SJ"<<endl;
    cout << "risk: "<<endl;
    cout << r3 <<endl;
    cout<<endl;
    cout<<"SJ N1"<<endl;
    cout << "risk: "<<endl;
    cout << r4 <<endl;
    cout<<endl;

    cout<<"\n\n";
    cout<<" --- \t --- COST --- \t ---\n"<<endl;
    cout<<"OPF"<<endl;
    cout << "costs: "<<endl;
    cout << o0 <<endl;
    cout<<endl;
    cout<<"OPF n1"<<endl;
    cout << "costs: "<<endl;
    cout << o1 <<endl;
    cout<<endl;


    cout<<"CC"<<endl;
    cout << "costs: "<<endl;
    cout << o <<endl;
    cout<<endl;


    cout<<"SJ"<<endl;
    cout << "costs: "<<endl;
    cout << o3 <<endl;
    cout<<endl;

    cout<<"SJ N1"<<endl;
    cout << "costs: "<<endl;
    cout << o4 <<endl;
    cout<<endl;




    /*
    ofstream myopf( "opf.json" );
    ofstream mycc( "cc.json" );
    ofstream myjcc( "jcc.json" );

    jsonInter ji;
    //    ji.print(gr,rsj);
    ji.printDetail(myopf,gr,rnom,L,p,sqrt(TV),o0,r0,f0,g0,beta0,sd0,z0,p0);
    ji.printDetail(mycc,gr,rcc,L,p,sqrt(TV),o,r,f,g,beta,sd,z,p1);
    ji.printDetail(myjcc,gr,rsjn,L,p,sqrt(TV),o4,r4,f4,g4,beta4,sd4,z4,p4);
    //    ji.printDetail(cerr,gr,rsj,L,p,o4,r4,f4,g4,beta4,sd4,z4,p4);

    myopf.close();
    mycc.close();
    myjcc.close();
    */

    if(argc>13){
    check.t().print("check: ");
        
    ofstream myopf( "opf.out" );
    ofstream myopfn1( "opfn1.out" );
    ofstream mycc( "cc.out" );
    ofstream myjcc( "jcc.out" );
    ofstream myjccn1( "jccn1.out" );
    iopa opa(gr,&gc,opaL,opap,L,p);

    opa.runTrials(myopf, &n1, f0, g0, SIGy0,T,o0);
    opa.runTrials(myopfn1, &n1, f1, g1, SIGy1,T,o1);
    opa.runTrials(mycc, &n1, f, g, SIGycc,T,o);
    opa.runTrials(myjcc, &n1, f3, g3, SIGy3,T,o3);
    opa.runTrials(myjccn1, &n1, f4, g4, SIGy4,T,o4);

    /*    for(int n=0;n<Nl;n++){
      if(check(n)){
	vec f1n = n1.getN1(n,f1,g1);
	vec z1n=gc.risk(f1n,sd1,L,p,pc);
	z1n(n)=1;
	z1n.t().print("Risk: ");	
	//	f1n.print("f35: ");
	opa.runTrials(myopa, z1n,n,T);
      }
    }
    */
    myopf.close();
    myopfn1.close();
    mycc.close();
    myjcc.close();
    myjccn1.close();

    ofstream mycomp( "comp.out" );

    mycomp<<" --- \t --- RISK 2nd--- \t ---\n"<<endl;

	mycomp<<"OPF"<<"\t"<<s0<<endl;
	mycomp<<"C0: "<<o0<<endl;
	mycomp<<"r0 - "<<r0<<endl;
	mycomp << "count = " << stats_r0.count() << endl;
	mycomp << "mean = " << stats_r0.mean() << endl;
	mycomp << "stdv  = " << stats_r0.stddev()  << endl;
	mycomp << "min  = " << stats_r0.min()  << endl;
	mycomp << "max  = " << stats_r0.max()  << endl;
	mycomp<<endl;


	mycomp<<"OPF n1"<<"\t"<<s1<<endl;
	mycomp<<"C1: "<<o1<<endl;
	mycomp<<"r1 - "<<r1<<endl;
	mycomp << "count = " << stats_r1.count() << endl;
	mycomp << "mean = " << stats_r1.mean() << endl;
	mycomp << "stdv  = " << stats_r1.stddev()  << endl;
	mycomp << "min  = " << stats_r1.min()  << endl;
	mycomp << "max  = " << stats_r1.max()  << endl;
	mycomp<<endl;


	mycomp<<"CC"<<"\t"<<s<<endl;
	mycomp<<"C: "<<o<<endl;
	mycomp<<"r - "<<r<<endl;
	mycomp << "count = " << stats_risk.count() << endl;
	mycomp << "mean = " << stats_risk.mean() << endl;
	mycomp << "stdv  = " << stats_risk.stddev()  << endl;
	mycomp << "min  = " << stats_risk.min()  << endl;
	mycomp << "max  = " << stats_risk.max()  << endl;
	mycomp<<endl;

	mycomp<<"SJ"<<"\t"<<s3<<endl;
	mycomp<<"C3: "<<o3<<endl;
	mycomp<<"r3 - "<<r3<<endl;
	mycomp << "count = " << stats_r3.count() << endl;
	mycomp << "mean = " << stats_r3.mean() << endl;
	mycomp << "stdv  = " << stats_r3.stddev()  << endl;
	mycomp << "min  = " << stats_r3.min()  << endl;
	mycomp << "max  = " << stats_r3.max()  << endl;
	mycomp<<endl;

	mycomp<<"SJ N1"<<"\t"<<s4<<endl;
	mycomp<<"C4: "<<o4<<endl;
	mycomp<<"r4 - "<<r4<<endl;
	mycomp << "count = " << stats_r4.count() << endl;
	mycomp << "mean = " << stats_r4.mean() << endl;
	mycomp << "stdv  = " << stats_r4.stddev()  << endl;
	mycomp << "min  = " << stats_r4.min()  << endl;
	mycomp << "max  = " << stats_r4.max()  << endl;
	mycomp<<endl;

	mycomp<<"\n\n";
    
	mycomp.close();

    }
      }
      catch(IloException& e){
	cerr<<"Concert exception: "<<e<<endl; 
      }
      catch(exception& e){
	cerr<<"Exception: "<<e.what()<<endl; 
      }
      catch(...){
	cerr<<"Unknown Error "<<endl;
      }



  return 0;



}

   

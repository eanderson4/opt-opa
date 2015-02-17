#include <stdlib.h>

#include <ctime>
#include "sqlinter.h"
#include "jsoninter.h"

#include "icc.h"
#include "isj.h"
#include "isjn.h"
#include "ijn1.h"
#include "in1.h"

using namespace std;

int main(int argc, char* argv[]){

  if(argc<=12){
    cout<<"cmd: pow case/30.db <m0> <m1> <e0> <e1> <L> <p> <B> <T>\n"
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
	<<"\t<T> num trials"<<endl;
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
  
  
  

  

    ijn1 n1(gr,  SIGy,Hw,L,p,pc,eps,epsN);
  n1.addCost();
  vec check = n1.getCheck();

    isj nom(gr, &gc, SIG, indexM, L, p, pc, 1);
    nom.lineLimitStatus(true);




      




      try{
	//	igrid nom(gr);
	//	nom.addCost();

	rgrid * rnom = nom.solveModel(&is);

	

	
	double o0=rnom->getObjective();
	vec f0=gc.convert(rnom->getF());
	vec g0=gc.convert(rnom->getG());
	vec beta0=nom.getBeta();
	vec z0=gc.risk(f0,SIGy.diag(),L,p,pc);
	double r0 = sum(z0);
	vec p0=gc.lineprob(f0,SIGy.diag());
	IloCplex::CplexStatus s0=rnom->getStatus();
	TG = accu(g0);
	vec sd0 = nom.getSD();

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
	    vec z0n=gc.risk(f0n,SIGy.diag(),L,p,pc);
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
	


	

	icc cc(gr, &gc, SIG, indexM, epsL, epsG);
  	rgrid * rcc = cc.solveModel(&is);
	
	double o=rcc->getObjective();
	vec f=gc.convert(rcc->getF());
	vec g=gc.convert(rcc->getG());
	vec beta=cc.getBeta();
	vec sd=cc.getSD();
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
	    vec zn=gc.risk(fn,sd,L,p,pc);
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
	

	

	isj sj(gr, &gc, SIG, indexM, L, p, pc, eps);
	rgrid * rsj = sj.solveModel(&is);

	
	double o4=rsj->getObjective();
	vec f4=gc.convert(rsj->getF());
	vec g4=gc.convert(rsj->getG());
	vec beta4=sj.getBeta();
	vec sd4=sj.getSD();
	vec z4=gc.risk(f4,sd4,L,p,pc);
	double r4 = sum(z4);
	vec p4=gc.lineprob(f4,sd4);
	IloCplex::CplexStatus s4=rsj->getStatus();

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
	    vec z4n=gc.risk(f4n,sd4,L,p,pc);
	    double r4n = sum(z4n);
	    stats_r4(r4n);
	  }
	}
	
	cout<<"SJ"<<"\t"<<s4<<endl;
	cout<<"C4: "<<o4<<endl;
	cout<<"r4 - "<<r4<<endl;
	cout << "count = " << stats_r4.count() << endl;
	cout << "mean = " << stats_r4.mean() << endl;
	cout << "stdv  = " << stats_r4.stddev()  << endl;
	cout << "min  = " << stats_r4.min()  << endl;
	cout << "max  = " << stats_r4.max()  << endl;
	cout<<endl;
	

    

    
    
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
    cout<<" --- \t --- RISK --- \t ---\n"<<endl;
    cout<<"OPF"<<endl;
    cout << "risk: "<<endl;
    cout << r0 <<endl;
    cout<<endl;

    cout<<"CC"<<endl;
    cout << "risk: "<<endl;
    cout << r <<endl;
    cout<<endl;


    cout<<"SJ"<<endl;
    cout << "risk: "<<endl;
    cout << r4 <<endl;
    cout<<endl;

    cout<<"\n\n";
    cout<<" --- \t --- COST --- \t ---\n"<<endl;
    cout<<"OPF"<<endl;
    cout << "costs: "<<endl;
    cout << o0 <<endl;
    cout<<endl;

    cout<<"CC"<<endl;
    cout << "costs: "<<endl;
    cout << o <<endl;
    cout<<endl;


    cout<<"SJ"<<endl;
    cout << "costs: "<<endl;
    cout << o4 <<endl;
    cout<<endl;

    ofstream myopf( "opf.json" );
    ofstream mycc( "cc.json" );
    ofstream myjcc( "jcc.json" );

    jsonInter ji;
    //    ji.print(gr,rsj);
    ji.printDetail(myopf,gr,rnom,L,p,sqrt(TV),o0,r0,f0,g0,beta0,sd0,z0,p0);
    ji.printDetail(mycc,gr,rcc,L,p,sqrt(TV),o,r,f,g,beta,sd,z,p1);
    ji.printDetail(myjcc,gr,rsj,L,p,sqrt(TV),o4,r4,f4,g4,beta4,sd4,z4,p4);
    //    ji.printDetail(cerr,gr,rsj,L,p,o4,r4,f4,g4,beta4,sd4,z4,p4);

    myopf.close();
    mycc.close();
    myjcc.close();

    //    f4.print("f4: ");
    //    sqrt(sd4).print("sd4: ");
    //    z4.print("z4: ");

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

   

#include "iopa.h"


void iopa::runTrials(ostream & out,vec z,int N,double num){
  int Nl = z.n_elem; 
  std::clock_t start;
  double duration;

  start = std::clock();  
   
  igrid ig(_gr);
  ig.allowLoadShed();
  ig.addCost();

  
  isolve is;
  is.setSolver(IloCplex::Dual,IloCplex::Dual);

  
  IloCplex cplex(*ig.getModel());
  is.setCplexParams(&cplex);


  rgrid * rbase = ig.solveModel(&cplex);
  vec f0=_gc->convert(rbase->getF());
  vec g0=_gc->convert(rbase->getG());
  double d0=sum(g0);

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"Duration: "<<duration<<endl;
  running_stat<double> stats_duration;
  stats_duration(duration);
  is.setTime(1);
  is.setCplexParams(&cplex);

  
  running_stat<double> stats_ls;
  running_stat<double> stats_stages;
  IloNumArray flow(ig.getEnv());
  IloNumArray gen(ig.getEnv());

  for(int n=0;n<num;n++){
    del_g mod(_gr); 
    
    //    z.t().print("Initial Risk: ");
    for(int i=0;i<Nl;i++){
      double randu =((double)rand()/(double)RAND_MAX);
      if(z(i) >= randu)      mod.setStatus(i,false);
    }
    
    //    cout<<"Initial Failures: "<<endl;
    //    cout<<mod<<endl;
    

    /////////LOOP this shit!
    bool fail=true;
    bool record=true;
    int stages=0;
    double d;
    while(fail){
      stages++;
      fail=false;
      record=true;
      ig.modGrid(mod);
      
      start = std::clock();  
      
      if(cplex.solve()){
	
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//	cout<<"Duration: "<<duration<<endl;
	stats_duration(duration);
	//	is.setTime(stats_duration.max()*100);
	//	is.setCplexParams(&cplex); 

	cplex.getValues(flow,ig.getF());
	cplex.getValues(gen,ig.getG());
	
	vec f=_gc->convert(flow);
	vec g=_gc->convert(gen);
	d=sum(g);
	
	vec z =_gc->risk(f,_L,_p);
	
	//	f.t().print("ff: ");
	//	g.t().print("gf: ");
	//	z.t().print("zf: ");
	
	for(int i=0; i < Nl; i++){
	  if(z(i)>0){
	    double randu =((double)rand()/(double)RAND_MAX);
	    if (randu <= z(i)){
	      //fail line i
	      mod.setStatus(i,false);
	      fail=true;
	      cout<<"Line "<<i<<" failed"<<endl;
	    }
	  }    
	  
	}
	
      }
      else{ 
	//cplex failed to solve, catalog and move on
	cerr<<"Not solved"<<endl;
	cerr<<cplex.getStatus()<<endl;
	cerr<<cplex.getCplexStatus()<<endl;
	cerr<<mod<<endl;
	record=false;
      }
    }
     
    //    f0.t().print("f0: ");
    //    g0.t().print("g0: ");
    
    if(record){
      //      cout<<"Stages: "<<stages<<endl;
      //      cout<<mod<<endl;
      //      cout<<"d0: "<<d0<<endl;
      //      cout<<"df: "<<d<<endl;
      
      stats_ls(d0-d);
      stats_stages(stages);   
    }
    
    ig.unmodGrid(mod);
  }
  

  //output inmportatnt informationation!
  out<<"\n"<<endl;
  out<<"Contingency: "<<N<<endl;
  out<<"\n"<<endl;
  out<<"Load Shed: "<<endl;
  out << "count = " << stats_ls.count() << endl;
  out << "mean = " << stats_ls.mean() << endl;
  out << "stdv  = " << stats_ls.stddev()  << endl;
  out << "ster  = " << stats_ls.stddev()/sqrt(stats_ls.count())  << endl;
  out << "min  = " << stats_ls.min()  << endl;
  out << "max  = " << stats_ls.max()  << endl;
  out<<"\n"<<endl;
  out<<"Stages: "<<endl;
  out << "count = " << stats_stages.count() << endl;
  out << "mean = " << stats_stages.mean() << endl;
  out << "stdv  = " << stats_stages.stddev()  << endl;
  out << "min  = " << stats_stages.min()  << endl;
  out << "max  = " << stats_stages.max()  << endl;
  out<<"\n"<<endl;
  


  cplex.end();  


}


void iopa::runTrials(ostream & out,ijn1 * n1,vec f, vec g, mat SIG,double num,double cost){
  int Nl = f.n_elem; 
  std::clock_t start;
  double duration;

  start = std::clock();  
   
  igrid ig(_gr);
  ig.allowLoadShed();
  ig.addCost();

  out<<"N\tc\tr\tT\tLS\tSD\tSE\tmin\tmax"<<endl;

  
  isolve is;
  //  is.setSolver(IloCplex::Dual,IloCplex::Dual);

  
  IloCplex cplex(*ig.getModel());
  //  is.setCplexParams(&cplex);
  //  cplex.setParam( IloCplex::RootAlg, 0 ); 

  rgrid * rbase = ig.solveModel(&cplex);
  vec f0=_gc->convert(rbase->getF());
  vec g0=_gc->convert(rbase->getG());
  double d0=sum(g0);


  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"Duration: "<<duration<<endl;
  running_stat<double> stats_duration;
  stats_duration(duration);
  //  is.setTime(10);
  //  is.setCplexParams(&cplex);

  vec check = n1->getCheck();
  mat Lo = n1->getLo();
  vec z(Nl,fill::zeros);

  running_stat<double> statsT_ls;
  running_stat<double> statsT_stages;
  
  vec sd0=SIG.diag();
  for(int e=0;e<Nl;e++){
    double U = _gr->getBranch(e).getRateA();
    if(sd0(e)<0.00000001) sd0(e)=0;
  }
  z=_gc->risk(f,sd0,_Lr,_pr,.85);	
  double r0 = sum(z);

    for(int n=0;n<Nl;n++){
      if(check(n)){
	vec fn = n1->getN1(n,f,g);
	vec sdn(Nl);
	for(int e=0;e<Nl;e++){
	  double U = _gr->getBranch(e).getRateA();
	  sdn(e) = SIG(e,e) + 2*Lo(e,n)*SIG(e,n)+ Lo(e,n)*Lo(e,n)*SIG(n,n);
	  if(sdn(e)<0 && sdn(e)>=-.0000001) sdn(e)=0;

	}
	z=_gc->risk(fn,sdn,_Lr,_pr,.85);
	double r = sum(z);
	z(n)=1;
	//	z.t().print("Risk: ");	




  
  running_stat<double> stats_ls;
  running_stat<double> stats_stages;
  IloNumArray flow(ig.getEnv());
  IloNumArray gen(ig.getEnv());

  for(int n=0;n<num;n++){
    del_g mod(_gr); 
    
    //    z.t().print("Initial Risk: ");
    for(int i=0;i<Nl;i++){
      double randu =((double)rand()/(double)RAND_MAX);
      if(z(i) >= randu)      mod.setStatus(i,false);
    }
    
    //    cout<<"Initial Failures: "<<endl;
    //    cout<<mod<<endl;
    

    /////////LOOP this shit!
    bool fail=true;
    bool record=true;
    int stages=0;
    double d;
    while(fail){
      stages++;
      fail=false;
      record=true;
      ig.modGrid(mod);
      
      start = std::clock();  
      
      if(cplex.solve() || cplex.getCplexStatus() == IloCplex::NumBest){
	
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//	cout<<"Duration: "<<duration<<endl;
	stats_duration(duration);
	//	is.setTime(stats_duration.max()*100);
	//	is.setCplexParams(&cplex); 

	cplex.getValues(flow,ig.getF());
	cplex.getValues(gen,ig.getG());
	
	vec f=_gc->convert(flow);
	vec g=_gc->convert(gen);
	d=sum(g);
	
	vec z =_gc->risk(f,_L,_p);
	
	//	f.t().print("ff: ");
	//	g.t().print("gf: ");
	//	z.t().print("zf: ");
	
	for(int i=0; i < Nl; i++){
	  if(z(i)>0){
	    double randu =((double)rand()/(double)RAND_MAX);
	    if (randu <= z(i)){
	      //fail line i
	      mod.setStatus(i,false);
	      fail=true;
	      cout<<"Line "<<i<<" failed"<<endl;
	    }
	  }    
	  
	}
	
      }
      else{ 
	//cplex failed to solve, catalog and move on
	cerr<<"Not solved"<<endl;
	cerr<<cplex.getStatus()<<endl;
	cerr<<cplex.getCplexStatus()<<endl;
	cerr<<mod<<endl;
	record=false;

      }
    }
     
    //    f0.t().print("f0: ");
    //    g0.t().print("g0: ");
    
    if(record){
      //      cout<<"Stages: "<<stages<<endl;
      //      cout<<mod<<endl;
      //      cout<<"d0: "<<d0<<endl;
      //      cout<<"df: "<<d<<endl;
      double loadshed;
      if(d0-d < .000001) loadshed=0;
      else loadshed=d0-d;
      stats_ls(loadshed);
      stats_stages(stages);   
      statsT_ls(loadshed);
      statsT_stages(stages);   

      



    }
    
    ig.unmodGrid(mod);
  }
  

  //output inmportatnt informationation!
  /*  out<<"\n"<<endl;
  out<<"Contingency: "<<n<<endl;
  out<<"\n"<<endl;
  out<<"Load Shed: "<<endl;
  out << "count = " << stats_ls.count() << endl;
  out << "mean = " << stats_ls.mean() << endl;
  out << "stdv  = " << stats_ls.stddev()  << endl;
  out << "ster  = " << stats_ls.stddev()/sqrt(stats_ls.count())  << endl;
  out << "min  = " << stats_ls.min()  << endl;
  out << "max  = " << stats_ls.max()  << endl;
  out<<"\n"<<endl;
  out<<"Stages: "<<endl;
  out << "count = " << stats_stages.count() << endl;
  out << "mean = " << stats_stages.mean() << endl;
  out << "stdv  = " << stats_stages.stddev()  << endl;
  out << "min  = " << stats_stages.min()  << endl;
  out << "max  = " << stats_stages.max()  << endl;
  out<<"\n"<<endl;
  */
  out<<n<<"\t"<<cost<<"\t"<<r<<"\t"<<stats_ls.count()<<"\t"<<stats_ls.mean()<<"\t"<<stats_ls.stddev()<<"\t"<<stats_ls.stddev()/sqrt(stats_ls.count())<<"\t"<<stats_ls.min()<<"\t"<<stats_ls.max()<<endl;


      }
    }
    /*
    out<<"Total: "<<endl;
  out<<"\n"<<endl;
  out<<"Load Shed: "<<endl;
  out << "count = " << statsT_ls.count() << endl;
  out << "mean = " << statsT_ls.mean() << endl;
  out << "stdv  = " << statsT_ls.stddev()  << endl;
  out << "ster  = " << statsT_ls.stddev()/sqrt(statsT_ls.count())  << endl;
  out << "min  = " << statsT_ls.min()  << endl;
  out << "max  = " << statsT_ls.max()  << endl;
  out<<"\n"<<endl;
  out<<"Stages: "<<endl;
  out << "count = " << statsT_stages.count() << endl;
  out << "mean = " << statsT_stages.mean() << endl;
  out << "stdv  = " << statsT_stages.stddev()  << endl;
  out << "min  = " << statsT_stages.min()  << endl;
  out << "max  = " << statsT_stages.max()  << endl;
  out<<"\n"<<endl;
    */
    cerr<<statsT_ls.count()<<"\t"<<cost<<"\t"<<r0<<"\t"<<statsT_ls.mean()<<"\t"<<statsT_ls.stddev()<<"\t"<<statsT_ls.stddev()/sqrt(statsT_ls.count())<<"\t"<<statsT_ls.min()<<"\t"<<statsT_ls.max()<<endl;

  cplex.end();  

  rbase->displayOperatingPos(_gr);

}





void iopa::runTrials(ostream & out, ostream & out2, ijn1 * n1,vec f, vec g, mat SIG,double num,double cost){
  int Nl = f.n_elem; 
  std::clock_t start;
  double duration;

  start = std::clock();  
   
  igrid ig(_gr);
  ig.allowLoadShed();
  ig.addCost();

  out<<"N\tc\tr\tT\tLS\tSD\tSE\tmin\tmax"<<endl;

  
  isolve is;
  //  is.setSolver(IloCplex::Dual,IloCplex::Dual);

  
  IloCplex cplex(*ig.getModel());
  //  is.setCplexParams(&cplex);
  //  cplex.setParam( IloCplex::RootAlg, 0 ); 

  rgrid * rbase = ig.solveModel(&cplex);
  vec f0=_gc->convert(rbase->getF());
  vec g0=_gc->convert(rbase->getG());
  double d0=sum(g0);


  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"Duration: "<<duration<<endl;
  running_stat<double> stats_duration;
  stats_duration(duration);
  //  is.setTime(10);
  //  is.setCplexParams(&cplex);

  vec check = n1->getCheck();
  mat Lo = n1->getLo();
  vec z(Nl,fill::zeros);
  vec z0(Nl,fill::zeros);

  running_stat<double> statsT_ls;
  running_stat<double> statsT_stages;
  
  vec sd0=SIG.diag();
  for(int e=0;e<Nl;e++){
    double U = _gr->getBranch(e).getRateA();
    if(sd0(e)<0.00000001) sd0(e)=0;
  }
  z=_gc->risk(f,sd0,_Lr,_pr,.85);	
  double r0 = sum(z);

    for(int n=0;n<Nl;n++){
      if(check(n)){
	vec fn = n1->getN1(n,f,g);
	vec sdn(Nl);
	for(int e=0;e<Nl;e++){
	  double U = _gr->getBranch(e).getRateA();
	  sdn(e) = SIG(e,e) + 2*Lo(e,n)*SIG(e,n)+ Lo(e,n)*Lo(e,n)*SIG(n,n);
	  if(sdn(e)<0 && sdn(e)>=-.0000001) sdn(e)=0;

	}
	z=_gc->risk(fn,sdn,_Lr,_pr,.85);
	double r = sum(z);
	z(n)=1;
	//	z.t().print("Risk: ");	
	z0=z;



  
  running_stat<double> stats_ls;
  running_stat<double> stats_stages;
  IloNumArray flow(ig.getEnv());
  IloNumArray gen(ig.getEnv());

  for(int n=0;n<num;n++){
    del_g mod(_gr); 
    
    //    z.t().print("Initial Risk: ");
    for(int i=0;i<Nl;i++){
      double randu =((double)rand()/(double)RAND_MAX);
      if(z(i) >= randu)      mod.setStatus(i,false);
    }
    
    //    cout<<"Initial Failures: "<<endl;
    //    cout<<mod<<endl;
    

    /////////LOOP this shit!
    bool fail=true;
    bool record=true;
    int stages=0;
    double d;
    while(fail){
      stages++;
      fail=false;
      record=true;
      ig.modGrid(mod);
      
      start = std::clock();  
      
      if(cplex.solve() || cplex.getCplexStatus() == IloCplex::NumBest){
	
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//	cout<<"Duration: "<<duration<<endl;
	stats_duration(duration);
	//	is.setTime(stats_duration.max()*100);
	//	is.setCplexParams(&cplex); 

	cplex.getValues(flow,ig.getF());
	cplex.getValues(gen,ig.getG());
	
	vec f=_gc->convert(flow);
	vec g=_gc->convert(gen);
	d=sum(g);
	
	vec z =_gc->risk(f,_L,_p);
	
	//	f.t().print("ff: ");
	//	g.t().print("gf: ");
	//	z.t().print("zf: ");
	
	for(int i=0; i < Nl; i++){
	  if(z(i)>0){
	    double randu =((double)rand()/(double)RAND_MAX);
	    if (randu <= z(i)){
	      //fail line i
	      mod.setStatus(i,false);
	      fail=true;
	      cout<<"Line "<<i<<" failed"<<endl;
	    }
	  }    
	  
	}
	
      }
      else{ 
	//cplex failed to solve, catalog and move on
	cerr<<"Not solved"<<endl;
	cerr<<cplex.getStatus()<<endl;
	cerr<<cplex.getCplexStatus()<<endl;
	cerr<<mod<<endl;
	record=false;

      }
    }
     
    //    f0.t().print("f0: ");
    //    g0.t().print("g0: ");
    
    if(record){
      //      cout<<"Stages: "<<stages<<endl;
      //      cout<<mod<<endl;
      //      cout<<"d0: "<<d0<<endl;
      //      cout<<"df: "<<d<<endl;
      double loadshed;
      if(d0-d < .000001) loadshed=0;
      else loadshed=d0-d;
      stats_ls(loadshed);
      stats_stages(stages);   
      statsT_ls(loadshed);
      statsT_stages(stages);   


    }
    
    ig.unmodGrid(mod);
  }
  

  //output inmportatnt informationation!
  /*  out<<"\n"<<endl;
  out<<"Contingency: "<<n<<endl;
  out<<"\n"<<endl;
  out<<"Load Shed: "<<endl;
  out << "count = " << stats_ls.count() << endl;
  out << "mean = " << stats_ls.mean() << endl;
  out << "stdv  = " << stats_ls.stddev()  << endl;
  out << "ster  = " << stats_ls.stddev()/sqrt(stats_ls.count())  << endl;
  out << "min  = " << stats_ls.min()  << endl;
  out << "max  = " << stats_ls.max()  << endl;
  out<<"\n"<<endl;
  out<<"Stages: "<<endl;
  out << "count = " << stats_stages.count() << endl;
  out << "mean = " << stats_stages.mean() << endl;
  out << "stdv  = " << stats_stages.stddev()  << endl;
  out << "min  = " << stats_stages.min()  << endl;
  out << "max  = " << stats_stages.max()  << endl;
  out<<"\n"<<endl;
  */
  out2<<n<<" "<<r<<" "<<stats_ls.count()<<" "<<stats_ls.mean()<<" "<<stats_ls.stddev();
      for(int i=0;i<Nl;i++){
	if(z(i)>.0005)
	  out2<<" "<<i<<","<<z(i);
      }
      out2<<endl;



  out<<n<<"\t"<<cost<<"\t"<<r<<"\t"<<stats_ls.count()<<"\t"<<stats_ls.mean()<<"\t"<<stats_ls.stddev()<<"\t"<<stats_ls.stddev()/sqrt(stats_ls.count())<<"\t"<<stats_ls.min()<<"\t"<<stats_ls.max()<<endl;


      }
    }
    /*
    out<<"Total: "<<endl;
  out<<"\n"<<endl;
  out<<"Load Shed: "<<endl;
  out << "count = " << statsT_ls.count() << endl;
  out << "mean = " << statsT_ls.mean() << endl;
  out << "stdv  = " << statsT_ls.stddev()  << endl;
  out << "ster  = " << statsT_ls.stddev()/sqrt(statsT_ls.count())  << endl;
  out << "min  = " << statsT_ls.min()  << endl;
  out << "max  = " << statsT_ls.max()  << endl;
  out<<"\n"<<endl;
  out<<"Stages: "<<endl;
  out << "count = " << statsT_stages.count() << endl;
  out << "mean = " << statsT_stages.mean() << endl;
  out << "stdv  = " << statsT_stages.stddev()  << endl;
  out << "min  = " << statsT_stages.min()  << endl;
  out << "max  = " << statsT_stages.max()  << endl;
  out<<"\n"<<endl;
    */
    cerr<<statsT_ls.count()<<"\t"<<cost<<"\t"<<r0<<"\t"<<statsT_ls.mean()<<"\t"<<statsT_ls.stddev()<<"\t"<<statsT_ls.stddev()/sqrt(statsT_ls.count())<<"\t"<<statsT_ls.min()<<"\t"<<statsT_ls.max()<<endl;

  cplex.end();  

  rbase->displayOperatingPos(_gr);

}


double iopa::runTrials(ostream & out, ostream & out2, ostream & mycomp, ijn1 * n1,vec f, vec g, mat SIG,vec xdes,double num,double cost,int Nstart){
  int Nl = f.n_elem; 
  std::clock_t start;
  double duration;

  start = std::clock();  
   
  igrid ig(_gr);
  ig.allowLoadShed();
  ig.addCost();

  out<<"N\tc\tr\trl\tT\tLS\tSD\tSE\tmin\tmax"<<endl;

  
  isolve is;
  //  is.setSolver(IloCplex::Dual,IloCplex::Dual);

  
  IloCplex cplex(*ig.getModel());
  //  is.setCplexParams(&cplex);
  //  cplex.setParam( IloCplex::RootAlg, 0 ); 

  rgrid * rbase = ig.solveModel(&cplex);
  vec f0=_gc->convert(rbase->getF());
  vec g0=_gc->convert(rbase->getG());
  double d0=sum(g0);


  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"Duration: "<<duration<<endl;
  running_stat<double> stats_duration;
  stats_duration(duration);
  //  is.setTime(10);
  //  is.setCplexParams(&cplex);

  vec check = n1->getCheck();
  mat Lo = n1->getLo();
  vec z(Nl,fill::zeros);
  vec z0(Nl,fill::zeros);



  running_stat<double> statsT_ls;
  running_stat<double> statsT_stages;
  
  vec sd0=SIG.diag();
  for(int e=0;e<Nl;e++){
    double U = _gr->getBranch(e).getRateA();
    if(sd0(e)<0.00000001) sd0(e)=0;
  }
  z=_gc->risk(f,sd0,_Lr,_pr,.85);	
  double r0 = sum(z);
  vec l0 = xdes % z;

    for(int n=Nstart;n<Nl;n++){
      if(check(n)){
	vec fn = n1->getN1(n,f,g);
	vec sdn(Nl);
	for(int e=0;e<Nl;e++){
	  double U = _gr->getBranch(e).getRateA();
	  sdn(e) = SIG(e,e) + 2*Lo(e,n)*SIG(e,n)+ Lo(e,n)*Lo(e,n)*SIG(n,n);
	  if(sdn(e)<0 && sdn(e)>=-.0000001) sdn(e)=0;

	}
	z=_gc->risk(fn,sdn,_Lr,_pr,.85);
	double r = sum(z);

	z(n)=1;
	vec  l = xdes.t()*z;
	//	z.t().print("Risk: ");	
	z0=z;




  
  running_stat<double> stats_ls;
  running_stat<double> stats_stages;
  IloNumArray flow(ig.getEnv());
  IloNumArray gen(ig.getEnv());

  for(int n=0;n<num;n++){
    del_g mod(_gr); 
    
    //    z.t().print("Initial Risk: ");
    for(int i=0;i<Nl;i++){
      double randu =((double)rand()/(double)RAND_MAX);
      if(z(i) >= randu)      mod.setStatus(i,false);
    }
    
    //    cout<<"Initial Failures: "<<endl;
    //    cout<<mod<<endl;
    

    /////////LOOP this shit!
    bool fail=true;
    bool record=true;
    int stages=0;
    double d;
    while(fail){
      stages++;
      fail=false;
      record=true;
      ig.modGrid(mod);
      
      start = std::clock();  
      
      if(cplex.solve() || cplex.getCplexStatus() == IloCplex::NumBest){
	
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//	cout<<"Duration: "<<duration<<endl;
	stats_duration(duration);
	//	is.setTime(stats_duration.max()*100);
	//	is.setCplexParams(&cplex); 

	cplex.getValues(flow,ig.getF());
	cplex.getValues(gen,ig.getG());
	
	vec f=_gc->convert(flow);
	vec g=_gc->convert(gen);
	d=sum(g);
	
	vec z =_gc->risk(f,_L,_p);
	

	//	f.t().print("ff: ");
	//	g.t().print("gf: ");
	//	z.t().print("zf: ");
	
	for(int i=0; i < Nl; i++){
	  if(z(i)>0){
	    double randu =((double)rand()/(double)RAND_MAX);
	    if (randu <= z(i)){
	      //fail line i
	      mod.setStatus(i,false);
	      fail=true;
	      cout<<"Line "<<i<<" failed"<<endl;
	    }
	  }    
	  
	}
	
      }
      else{ 
	//cplex failed to solve, catalog and move on
	cerr<<"Not solved"<<endl;
	cerr<<cplex.getStatus()<<endl;
	cerr<<cplex.getCplexStatus()<<endl;
	cerr<<mod<<endl;
	record=false;

      }
    }
     
    //    f0.t().print("f0: ");
    //    g0.t().print("g0: ");
    
    if(record){
      //      cout<<"Stages: "<<stages<<endl;
      //      cout<<mod<<endl;
      //      cout<<"d0: "<<d0<<endl;
      //      cout<<"df: "<<d<<endl;
      double loadshed;
      if(d0-d < .000001) loadshed=0;
      else loadshed=d0-d;
      stats_ls(loadshed);
      stats_stages(stages);   
      statsT_ls(loadshed);
      statsT_stages(stages);   


    }
    
    ig.unmodGrid(mod);
  }
  

  //output inmportatnt informationation!
  /*  out<<"\n"<<endl;
  out<<"Contingency: "<<n<<endl;
  out<<"\n"<<endl;
  out<<"Load Shed: "<<endl;
  out << "count = " << stats_ls.count() << endl;
  out << "mean = " << stats_ls.mean() << endl;
  out << "stdv  = " << stats_ls.stddev()  << endl;
  out << "ster  = " << stats_ls.stddev()/sqrt(stats_ls.count())  << endl;
  out << "min  = " << stats_ls.min()  << endl;
  out << "max  = " << stats_ls.max()  << endl;
  out<<"\n"<<endl;
  out<<"Stages: "<<endl;
  out << "count = " << stats_stages.count() << endl;
  out << "mean = " << stats_stages.mean() << endl;
  out << "stdv  = " << stats_stages.stddev()  << endl;
  out << "min  = " << stats_stages.min()  << endl;
  out << "max  = " << stats_stages.max()  << endl;
  out<<"\n"<<endl;
  */
  out2<<n<<" "<<r<<" "<<stats_ls.count()<<" "<<stats_ls.mean()<<" "<<stats_ls.stddev();
      for(int i=0;i<Nl;i++){
	if(z(i)>.0005)
	  out2<<" "<<i<<","<<z(i);
      }
      out2<<endl;


      double rl = accu(l);
      out<<n<<"\t"<<cost<<"\t"<<r<<"\t"<<rl<<"\t"<<stats_ls.count()<<"\t"<<stats_ls.mean()<<"\t"<<stats_ls.stddev()<<"\t"<<stats_ls.stddev()/sqrt(stats_ls.count())<<"\t"<<stats_ls.min()<<"\t"<<stats_ls.max()<<endl;


      }
    }
    /*
    out<<"Total: "<<endl;
  out<<"\n"<<endl;
  out<<"Load Shed: "<<endl;
  out << "count = " << statsT_ls.count() << endl;
  out << "mean = " << statsT_ls.mean() << endl;
  out << "stdv  = " << statsT_ls.stddev()  << endl;
  out << "ster  = " << statsT_ls.stddev()/sqrt(statsT_ls.count())  << endl;
  out << "min  = " << statsT_ls.min()  << endl;
  out << "max  = " << statsT_ls.max()  << endl;
  out<<"\n"<<endl;
  out<<"Stages: "<<endl;
  out << "count = " << statsT_stages.count() << endl;
  out << "mean = " << statsT_stages.mean() << endl;
  out << "stdv  = " << statsT_stages.stddev()  << endl;
  out << "min  = " << statsT_stages.min()  << endl;
  out << "max  = " << statsT_stages.max()  << endl;
  out<<"\n"<<endl;
    */
    cerr<<"Load Shed: "<<endl;
        cerr<<statsT_ls.count()<<"\t"<<cost<<"\t"<<r0<<"\t"<<statsT_ls.mean()<<"\t"<<statsT_ls.stddev()<<"\t"<<statsT_ls.stddev()/sqrt(statsT_ls.count())<<"\t"<<statsT_ls.min()<<"\t"<<statsT_ls.max()<<endl;
	cerr<<"Stages: "<<endl;
    cerr<<statsT_stages.count()<<"\t"<<cost<<"\t"<<r0<<"\t"<<statsT_stages.mean()<<"\t"<<statsT_stages.stddev()<<"\t"<<statsT_stages.stddev()/sqrt(statsT_stages.count())<<"\t"<<statsT_stages.min()<<"\t"<<statsT_stages.max()<<endl;


    mycomp<<"Load Shed: (count cost risk ls_mean ls_stddv ls_stderr ls_min ls_max)"<<endl;
        mycomp<<statsT_ls.count()<<"\t"<<cost<<"\t"<<r0<<"\t"<<statsT_ls.mean()<<"\t"<<statsT_ls.stddev()<<"\t"<<statsT_ls.stddev()/sqrt(statsT_ls.count())<<"\t"<<statsT_ls.min()<<"\t"<<statsT_ls.max()<<endl;
	mycomp<<"Stages: (count cost risk S_mean S_stddv S_stderr S_min S_max)"<<endl;
    mycomp<<statsT_stages.count()<<"\t"<<cost<<"\t"<<r0<<"\t"<<statsT_stages.mean()<<"\t"<<statsT_stages.stddev()<<"\t"<<statsT_stages.stddev()/sqrt(statsT_stages.count())<<"\t"<<statsT_stages.min()<<"\t"<<statsT_stages.max()<<endl;


  cplex.end();  
  

  return statsT_ls.mean();
  //  rbase->displayOperatingPos(_gr);

}




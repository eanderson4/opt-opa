#include "iopa.h"


void iopa::runTrials(ostream & out,vec z,int N,double num, running_stat<double> st_tot){
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

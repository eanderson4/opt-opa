#include "iopa.h"


void iopa::runTrials(double num){
  int Nl = _z.n_elem; 
  
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
  
  running_stat<double> stats_ls;
  running_stat<double> stats_stages;

  for(int n=0;n<num;n++){
    del_g mod(_gr); 
    
    _z.t().print("Initial Risk: ");
    for(int i=0;i<Nl;i++){
      double randu =((double)rand()/(double)RAND_MAX);
      if(_z(i) >= randu)      mod.setStatus(i,false);
    }
    
    cout<<"Initial Failures: "<<endl;
    cout<<mod<<endl;
    

    /////////LOOP this shit!
    bool fail=true;
    int stages=0;
    double d;
    while(fail){
      stages++;
      fail=false;
      ig.modGrid(mod);
      
      rgrid * rg = ig.solveModel(&cplex);
      
      vec f=_gc->convert(rg->getF());
      vec g=_gc->convert(rg->getG());
      d=sum(g);
      
      vec z =_gc->risk(f,_L,_p);
      
      f.t().print("ff: ");
      g.t().print("gf: ");
      z.t().print("zf: ");
      
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
    
    
    f0.t().print("f0: ");
    g0.t().print("g0: ");
    
    
    cout<<"Stages: "<<stages<<endl;
    cout<<mod<<endl;
    cout<<"d0: "<<d0<<endl;
    cout<<"df: "<<d<<endl;
    
    stats_ls(d0-d);
    stats_stages(stages);

    ig.unmodGrid(mod);
  }
  

  //output inmportatnt informationation!
  cout<<"\n"<<endl;
  cout<<"Load Shed: "<<endl;
  cout << "count = " << stats_ls.count() << endl;
  cout << "mean = " << stats_ls.mean() << endl;
  cout << "stdv  = " << stats_ls.stddev()  << endl;
  cout << "ster  = " << stats_ls.stddev()/sqrt(stats_ls.count())  << endl;
  cout << "min  = " << stats_ls.min()  << endl;
  cout << "max  = " << stats_ls.max()  << endl;
  cout<<"\n"<<endl;
  cout<<"Stages: "<<endl;
  cout << "count = " << stats_stages.count() << endl;
  cout << "mean = " << stats_stages.mean() << endl;
  cout << "stdv  = " << stats_stages.stddev()  << endl;
  cout << "min  = " << stats_stages.min()  << endl;
  cout << "max  = " << stats_stages.max()  << endl;
  cout<<"\n"<<endl;
  


  cplex.end();  


}

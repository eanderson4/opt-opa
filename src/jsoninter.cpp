#include "jsoninter.h"

void jsonInter::print(grid * gr, rgrid * rg) {
  cout<<"Print to JSON"<<endl;

  int nB=gr->numBuses();
  int nR=gr->numBranches();
  int nG=gr->numGens();

  //Outputing resolts to JSON format
  cerr<<"{\n"
      <<"\t\"dataset\": {  \n"
      <<"\t\t\"nodes\": [ \n";

  //Buses
  for(int i =0; i<nB; i++){
    bus b = gr->getBus(i);
    double ni=0; //net inject
    cerr<<"\t\t\t { \"name\": \""<<b.getNum()<<"\",\n"
	<<"\t\t\t   \"index\": "<<i<<",\n"
	<<"\t\t\t   \"demand\": "<<b.getP()<<",\n";
    ni=ni+b.getP();

    for(int j =0; j<nG; j++){
      gen g =gr->getGen(j); 
      int bus = g.getBus();
      if(bus==b.getNum()){
	double p = rg->getG()[j];
	if( p > 0.001 ){
	  cerr <<"\t\t\t   \"p\": "<<p<<",\n";
	  ni=ni-p;
	}
	else{
	  cerr <<"\t\t\t   \"p\": "<<0<<",\n";
	}
      }
    }
    cerr <<"\t\t\t   \"ni\": "<<ni<<" ";
    if(i!=nB-1) cerr<<"},"<<endl;
    else  cerr<<"} \n"
	      <<"\t\t],"<<endl;

  }

  //Branches
  cerr<<"\t\t \"edges\": [ \n";
  for(int i =0; i<nR; i++){
    branch b = gr->getBranch(i);
    int from= gr->getBusNum(b.getFrom());
    int to= gr->getBusNum(b.getTo());
    double X=gr->getBranch(i).getX();
    double cap=gr->getBranch(i).getRateA();
    double flow=double(rg->getF()[i]);
    cerr<<"\t\t\t { \"source\": "<<from<<", \"target\": "<<to<<","<<endl;
    cerr<<"\t\t\t   \"index\": "<<i<<",\n";
    cerr<<"\t\t\t   \"X\": "<<X<<",\n";
    cerr<<"\t\t\t   \"flow\": "<<flow<<",\n";;
    cerr<<"\t\t\t   \"cap\": "<<cap;
    if (i!=nR-1) cerr<<" },"<<endl;
    else cerr<<" }"<<endl;
  }
  cerr<<"\t\t ],\n"
      <<"\t\t \"gens\": [ \n";
  for(int i =0; i<nG; i++){
    gen g =gr->getGen(i); 
    int bus = g.getBus();
    double p = g.getP();
    double pmax = g.getPmax();
    double pmin = g.getPmin();
    double c2 = g.getC2();
    double c1 = g.getC1();
    double c0 = g.getC0();
    double cost = p*p*c2 + p*c1 + c0;
    if( p > 0.001 ){
      cerr<<"\t\t\t { \"bus\": "<<bus<<",\n"
	  <<"\t\t\t   \"p\": "<<p<<",\n"
	  <<"\t\t\t   \"pmax\": "<<pmax<<",\n"
	  <<"\t\t\t   \"pmin\": "<<pmin<<",\n"
	  <<"\t\t\t   \"cost\": "<<cost<<",\n"
	  <<"\t\t\t   \"c2\": "<<c2<<",\n"
	  <<"\t\t\t   \"c1\": "<<c1<<",\n"
	  <<"\t\t\t   \"c0\": "<<c0;
      if (i!=nG-1) cerr<<" },"<<endl;
      else cerr<<" }"<<endl;
    }
  }
  cerr<<"\t\t ]\n"
      <<"\t} \n"
      <<"}"<<endl;
      
}


void jsonInter::printDetail(ostream & out, grid * gr,rgrid * rg,double L, double parm, double SD, double o, double r, vec f, vec g, vec beta, vec sd, vec z, vec p) {
  cout<<"Print to JSON"<<endl;

  int nB=gr->numBuses();
  int nR=gr->numBranches();
  int nG=gr->numGens();
  
  vec fU(nR);
  vec sdU(nR);
  for(int i=0;i<nR;i++){
    double U = gr->getBranch(i).getRateA();
    fU(i) = abs(f(i))/U;
    sdU(i) = sd(i)/U;
  }
  

  //Outputing resolts to JSON format
  out<<"{\n"
      <<"\t\"dataset\": {  \n"
      <<"\t\t\"cost\": "<<o<<", \n"
      <<"\t\t\"risk\": "<<r<<", \n"
      <<"\t\t\"L\": "<<L<<", \n"
      <<"\t\t\"p\": "<<parm<<", \n"
      <<"\t\t\"SD\": "<<SD<<", \n"
      <<"\t\t\"nodes\": [ \n";

  //Buses
  for(int i =0; i<nB; i++){
    bus b = gr->getBus(i);
    double ni=0; //net inject
    out<<"\t\t\t { \"name\": \""<<b.getNum()<<"\",\n"
	<<"\t\t\t   \"index\": "<<i<<",\n"
	<<"\t\t\t   \"demand\": "<<b.getP()<<",\n";
    ni=ni+b.getP();
    
    bool generator=false;
    for(int j =0; j<nG; j++){
      gen g =gr->getGen(j); 
      int bus = g.getBus();
      if(bus==b.getNum()){
	generator=true;
	double p = rg->getG()[j];
	if( p > 0.001 ){
	  out <<"\t\t\t   \"p\": "<<p<<",\n";
	  ni=ni-p;
	}
	else{
	  out <<"\t\t\t \"p\": "<<0<<",\n";
	}
      }
    }
    if(!generator){
	  out <<"\t\t\t \"p\": "<<0<<",\n";
    }
    out <<"\t\t\t   \"ni\": "<<ni<<" ";
    if(i!=nB-1) out<<"},"<<endl;
    else  out<<"} \n"
	      <<"\t\t],"<<endl;

  }

  //Branches
  out<<"\t\t \"edges\": [ \n";
  for(int i =0; i<nR; i++){
    branch b = gr->getBranch(i);
    int from= gr->getBusNum(b.getFrom());
    int to= gr->getBusNum(b.getTo());
    double X=gr->getBranch(i).getX();
    double cap=gr->getBranch(i).getRateA();
    double flow=double(rg->getF()[i]);
    out<<"\t\t\t { \"source\": "<<from<<", \"target\": "<<to<<","<<endl;
    out<<"\t\t\t   \"index\": "<<i<<",\n";
    out<<"\t\t\t   \"X\": "<<X<<",\n";
    out<<"\t\t\t   \"flow\": "<<flow<<",\n";
    out<<"\t\t\t   \"sd\": "<<sqrt(sd[i])<<",\n";
    out<<"\t\t\t   \"z\": "<<z[i]<<",\n";
    out<<"\t\t\t   \"prob\": "<<p[i]<<",\n";
    out<<"\t\t\t   \"cap\": "<<cap;
    if (i!=nR-1) out<<" },"<<endl;
    else out<<" }"<<endl;
  }
  out<<"\t\t ],\n"
      <<"\t\t \"gens\": [ \n";
  for(int i =0; i<nG; i++){
    gen g =gr->getGen(i); 
    int bus = g.getBus();
    double pow = rg->getG()[i];
    double pmax = g.getPmax();
    double pmin = g.getPmin();
    double c2 = g.getC2();
    double c1 = g.getC1();
    double c0 = g.getC0();
    double cost = pow*pow*c2 + pow*c1 + c0;
    if( pow > 0.001 ){
      out<<"\t\t\t { \"bus\": "<<bus<<",\n"
	  <<"\t\t\t   \"busindex\": "<<gr->getGenBus(i)<<",\n"
	  <<"\t\t\t   \"p\": "<<pow<<",\n"
	  <<"\t\t\t   \"pmax\": "<<pmax<<",\n"
	  <<"\t\t\t   \"pmin\": "<<pmin<<",\n"
	  <<"\t\t\t   \"cost\": "<<cost<<",\n"
	  <<"\t\t\t   \"beta\": "<<beta[i]<<",\n"
	  <<"\t\t\t   \"c2\": "<<c2<<",\n"
	  <<"\t\t\t   \"c1\": "<<c1<<",\n"
	  <<"\t\t\t   \"c0\": "<<c0;
      if (i!=nG-1) out<<" },"<<endl;
      else out<<" }"<<endl;
    }
  }
  out<<"\t\t ]\n"
      <<"\t} \n"
      <<"}"<<endl;
      
}

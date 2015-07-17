#include <stdlib.h>

//#include <ctime>
#include "sqlinter.h"
//#include "jsoninter.h"

#include "gridcalc.h"
//#include "icc.h"
//#include "isj.h"
//#include "isjn.h"
//#include "ioj.h"
//#include "ijn1.h"
//#include "in1.h"
//#include "iopa.h"

using namespace std;

int main(int argc, char* argv[]){

  if(argc<=2){
    cout<<"cmd: design case/30.db <name input>\n"
	<<"\trun main for case30\n"
	<<"\t<name input> \n"<<endl;
    return 1;
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


  gridcalc gc(gr);  

  mat A = gc.getH();
  mat Lo = gc.getL(A);  

  mat C = gc.getC();


  vec check(Nl,fill::ones);
  mat Hb=A*trans(C);

  for(int n=0;n<Nl;n++){  //Line i outage
    for(int j=0;j<Nl;j++){
      if(Hb(n,n)<=1-.000001 || Hb(n,n)>=1+.0000001){
      }
      else{
	check(n)=0;
      }
    }
  }
  cout<<sum(check)<<" / "<<gr->numBranches()<<endl;






  // READ IN INPUT FILE and calculating weighting factor
  int ctr=0;
  vec indexmap(Nl,fill::zeros);
  vec indexmapback(Nl,fill::zeros);
  for(int i=0;i<Nl;i++){
    if(check(i)){
      indexmap(i)=ctr;
      indexmapback(ctr)=i;
      ctr++;
    }
  }
  cout<<ctr<<endl;
  //  return 0;

  int ctr2=0;
  vec indexctr(ctr,fill::zeros);
  vec riskctr(ctr,fill::zeros);
  vec risk2ctr(ctr,fill::zeros);
  vec risk3ctr(ctr,fill::zeros);
  vec trialctr(ctr,fill::zeros);
  vec meanctr(ctr,fill::zeros);
  vec stdvctr(ctr,fill::zeros);

  mat R(ctr,Nl,fill::zeros);

  string line;

  string sstr("data/s-");
  string dstr("data/d-");
  string desstr("data/des-");
  sstr += argv[2];
  dstr += argv[2];
  desstr += argv[2];
  sstr += ".out";
  dstr += ".out";
  desstr += ".out";

  ifstream sfile (sstr.c_str());
  ifstream dfile (dstr.c_str());
  istringstream infilestream;
  istringstream instream2;
  if (sfile.is_open() && dfile.is_open())
    {
      while ( getline (sfile,line) )
	{
	  cout << line << '\n';
	}
      sfile.close();
      while ( getline (dfile,line) )
	{
	  infilestream.clear();
	  infilestream.str(line);
	  	  cout << line << '\n';
	  infilestream >> indexctr(ctr2) 
		       >> riskctr(ctr2)   // Posible change to file structure, two riskctrs
		   >> trialctr(ctr2) 
		   >> meanctr(ctr2) 
		   >> stdvctr(ctr2);
	  	  cout << indexctr(ctr2) <<"\t"
	  	       << riskctr(ctr2) <<"\t"
	  	       << trialctr(ctr2) <<"\t"
	  	       << meanctr(ctr2) <<"\t"
	  	       << stdvctr(ctr2)<<endl;
	  while(!infilestream.eof()){
	    string a;
	    infilestream >> a;
	    instream2.clear();
	    instream2.str(a);
	    string pos;
	    string val;
	    double ctrpos;
	    double ctrmap;
	    double ctrvalue;
	    getline( instream2, pos, ',');
	    instream2 >> val;
	    //	    	    cout<<pos<<" "<<val<<"\t";
	    ctrpos=atoi(pos.c_str());
	    ctrmap=indexmap(ctrpos);
	    ctrvalue=atof(val.c_str());
	    	    cout<<ctr2<<","<<ctrpos<<" "<<ctrvalue<<"\t";
	    R(ctr2,ctrpos)=ctrvalue;
	  }
	  	  cout<<endl;
		  //		  	  int jd;
		  //	  cin>>jd;

	  ctr2++;
	}
      dfile.close();
    }

  else cout << "Unable to open file";





  //  int test;
  //  cin>>test;

    sp_mat sR(R);
  
    //    sR.print("R mat: ");
    //    meanctr.t().print("mean: ");
    
    sR.row(0).print("row0: ");
    sR.row(1).print("row1: ");
    sR.row(2).print("row2: ");
    sR.row(3).print("row3: ");

    vec xdesign = solve(R,meanctr);

  //  xdesign.t().print("xdesign: ");
  cout<<"Design weighting (first 10)"<<endl;
  for(int ci=0;ci<10;ci++) {
    cout<<xdesign(ci)<<"\t"<<meanctr(ci)<<endl;
  }
  cout<<endl;
  
  double avgx = mean(xdesign);

  vec xdes(Nl,fill::zeros);
  
  int ctr3=0;
  for(int i=0;i<Nl;i++){
    if(check(i)){
      xdes(i) = xdesign(ctr3);
      ctr3++;
    }
    else
      xdes(i) = avgx;

  }




  vec resid = R*xdesign - meanctr;

  //  cin>>test;
  
  if(accu(resid)>.000001)  resid.t().print("resid: ");
  
  ofstream mydes( desstr.c_str() );
  
  int ctr4=0;
  for(int i=0;i<Nl;i++){
    if(check(i)){
            cout<<i<<":("<<xdesign(i)<<","<<meanctr(ctr4)<<")"<<endl;
      mydes<<i<<" "<<xdesign(i)<<endl;
      ctr4++;
    }
    else{
      cout<<i<<":"<<xdesign(i)<<"\t"<<avgx<<endl;
      //      mydes<<i<<" "<<xdesign(i)<<endl;
      mydes<<i<<" "<<avgx<<endl;
    }
  }


  return 0;

}  

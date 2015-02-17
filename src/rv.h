#ifndef RV_H
#define RV_H


#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <random>
#include <vector>
#include <sstream>
#include <iomanip>
#include <stdexcept>

using namespace std;



class ranvar
{

 public:
 ranvar(int seed=-372):_seed(seed) { 
    _tol=pow(10,-7);
    if(seed==-372) seed=time(NULL); 
  }
  ~ranvar() {}

  void createRV(int N,double mean, double stdv);
  int getNum(){ return random_variable.size(); }
  double getValue(int N){ return random_variable[N]; }
  double phi(double x);
  double PHI(double x);
  double sampleMean();
  double sampleStDv();
  double testPhi();
  
  double g(double fnormalized, double L, double p, double pc);
  double ginv(double eps, double L, double p, double pc);
  double simProb(double L, double p, double pc);
  double anaProb(double L, double p, double pc){ return anaProb(L,p,pc,_mean,_stdv);}
  double anaProb(double L, double p, double pc, double mean, double stdv);
  double deriveMu(double L, double p, double pc, double mean, double stdv);
  double deriveSigma(double L, double p, double pc, double mean, double stdv);
  double d2Mu(double L, double p, double pc, double mean, double stdv);
  double d2Sigma(double L, double p, double pc, double mean, double stdv);
  double dSigmaMu(double L, double p, double pc, double mean, double stdv);

  double PHIInverse(double p);
  double RationalApproximation(double t);
  void demoInverse();

  class analerror: public exception 
  {
    virtual const char* what() const throw()
    {
      return "Analytic approximation failed!";
    }
  }analerr;
  


  friend ostream& operator<<(ostream& os, const ranvar& rv);

 private:

  vector<double> random_variable;
  int _seed;

  int _N;
  double _mean;
  double _stdv;
  
  double _tol;
  
};

#endif

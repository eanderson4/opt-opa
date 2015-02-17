#ifndef RBASE_H
#define RBASE_H

#include <ilcplex/ilocplex.h>

using namespace std;


class rbase {

 public:

  rbase ():env(IloEnv()){
    _f = IloNumArray(env);
    _g = IloNumArray(env);
    _theta = IloNumArray(env);
    _ls = IloNumArray(env);
  }
  ~rbase() {
    _f.end();
    _g.end();
    _theta.end();
    _ls.end();
    env.end();
  }

  IloNumArray getF(){ return _f; }
  IloNumArray getG(){ return _g; }
  IloNumArray getTheta(){ return _theta; }
  IloNumArray getLS(){ return _ls; }

 private:

  IloEnv env;

  IloNumArray _f;
  IloNumArray _g;
  IloNumArray _theta;
  IloNumArray _ls;

};

#endif

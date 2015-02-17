#ifndef GEN_H
#define GEN_H

#include <iostream>
using namespace std;

class gen {

 public:
  gen() {}
 gen(int num, int bus, double pg, double qg, double qmax, double qmin, double vg, double mbase, double status, double pmax, double pmin) : _num( num ), _bus( bus ), _pg( pg ), _qg( qg ), _qmax( qmax ), _qmin( qmin ), _vg( vg ), _mbase( mbase ), _status( status ), _pmax( pmax ), _pmin( pmin ) { }
  void addCost(int model, double startup, double shutdown, int ncost, double c2, double c1, double c0){
    _model=model;
    _startup=startup;
    _shutdown=shutdown;
    _ncost=ncost;
    _c2=c2;
    _c1=c1;
    _c0=c0;
  }

  int getBus(){ return _bus; }
  double getP(){ return _pg; }
  double getPmax(){ return _pmax; }
  double getPmin(){ return _pmin; }
  double getC2(){ return _c2; }
  double getC1(){ return _c1; }
  double getC0(){ return _c0; }
  double getStatus(){ return _status; }
  void setPmax(double pmax){ _pmax=pmax; }
  void setPmin(double pmin){ _pmin=pmin; }

  friend ostream& operator<<(ostream& os, const gen& g);

 private:
  int _num;
  int _bus;
  double _pg;
  double _qg;
  double _qmax;
  double _qmin;
  double _vg;
  double _mbase;
  double _status;
  double _pmax;
  double _pmin;

  int _model;
  double _startup;
  double _shutdown;
  int _ncost;
  double _c2;
  double _c1;
  double _c0;
  
};

#endif

#ifndef BRANCH_H
#define BRANCH_H

#include <iostream>
using namespace std;

class branch {

 public:
  branch(){}
 branch(int num, int fbus, int tbus, double br_r, double br_x, double br_b, double rate_a, double rate_b, double rate_c, double tap, double shift, int br_status, double ang_min, double ang_max) : _num( num ), _fbus( fbus ), _tbus( tbus ), _br_r( br_r ), _br_x( br_x ), _br_b( br_b ), _rate_a( rate_a ), _rate_b( rate_b ), _rate_c( rate_c ), _tap( tap ), _shift( shift ), _br_status( br_status ), _ang_min( ang_min ), _ang_max( ang_max ) {}

  void setRateA(double U){ _rate_a = U; }
  void setStatus(int status){ _br_status=status; }

  int getFrom(){ return _fbus; }
  int getTo(){ return _tbus; }
  double getShift(){ return _shift; }
  double getX(){    return _br_x;  }
  double getR(){ return _br_r; }
  double getTap(){ return _tap; }
  double getRateA(){ return _rate_a; }
  double getAngMax(){ return _ang_max; }
  double getAngMin(){ return _ang_min; }
  double getStatus(){ return _br_status; }

  friend ostream& operator<<(ostream& os, const branch& b);

 private:
  int _num;
  int _fbus;
  int _tbus;
  double _br_r;
  double _br_x;
  double _br_b;
  double _rate_a;
  double _rate_b;
  double _rate_c;
  double _tap;
  double _shift;
  int _br_status;
  double _ang_min;
  double _ang_max;
};


#endif

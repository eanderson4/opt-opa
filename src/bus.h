#ifndef BUS_H
#define BUS_H

#include <iostream>
using namespace std;

class bus {

 public:
   bus() {}
 bus(int num, int type, double pd, double qd, double gs, double bs, int area, double vm, double va, double basekv, int zone, double vmax, double vmin) : _num( num ), _type( type ), _pd( pd ), _qd( qd ), _gs( gs ), _bs( bs ), _area( area ), _vm( vm ), _va( va ), _basekv( basekv ), _zone( zone ), _vmax( vmax ), _vmin( vmin ) {} 

   void addP(double ap){ _pd=_pd+ap; }
   void setP(double pd){ _pd=pd; }
   
  int getNum(){ return _num; }
  int getType(){ return _type; }
  double getP(){    return _pd;  }
  double getGs(){ return _gs; }

  friend ostream& operator<<(ostream& os, const bus& b);

 private:

  int _num;
  int _type;
  double _pd;
  double _qd;
  double _gs;
  double _bs;
  int _area;
  double _vm;
  double _va;
  double _basekv;
  int _zone;
  double _vmax;
  double _vmin;

  double _p2;

};

#endif

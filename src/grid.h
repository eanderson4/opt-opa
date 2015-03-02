#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <vector>
#include <map>
#include "bus.h"
#include "branch.h"
#include "gen.h"
#include "gentype.h"

using namespace std;

class grid
{
  friend class bus;
  friend class branch;
  friend class gen;

 public:
 grid():mapBus( false ), loadShedPenalty(1000) {}
  ~grid() {}

  void setCapacity(int n, double U){    _branches[n].setRateA(U);  }
  void setGenPmax(int n, double U){    _gens[n].setPmax(U);  }
  void setGenPmin(int n, double U){    _gens[n].setPmin(U);  }
  void addPd(int n, double pd){ _buses[n].addP(pd); }
  void addBus(bus b){    _buses.push_back(b);  }
  void addBranch(branch b){    _branches.push_back(b);  }
  void addGen(gen g){    _gens.push_back(g);  }
  void addGenCost(int num, int model, double startup, double shutdown, int ncost, double c2, double c1, double c0);
  //  void modGrid(del_g mod);
  branch getBranch(int num){    return _branches[num];  }
  bus getBus(int num){    return _buses[num];  }
  gen getGen(int num){    return _gens[num];  }
  int numBuses(){    return _buses.size();  }
  int numBranches(){    return _branches.size();  }
  int numGens(){ return _gens.size();  }
  int getBusNum(int num);
  int getGenBus(int n){return getBusNum( _gens[n].getBus() ); }
  int getFromBus(int n){ return getBusNum( _branches[n].getFrom() ); }
  int getToBus(int n){ return getBusNum( _branches[n].getTo() ); }
  void buildMap();
  void printNums(ostream& stream);
  double getLoadShedPenalty(){ return loadShedPenalty; }
  double getTotalDemand();

  void addGenType(gentype gt){ _gt.push_back(gt); }
  double getC2(string gid){ return _gt[getGenType(gid)].getC2(); }
  double getC1(string gid){ return _gt[getGenType(gid)].getC1(); }
  int getGenType(string gid){ 
    for(unsigned int i=0;i<_gt.size();i++){
      if (gid == _gt[i].getGid()){
	return i;
      }
    }
    return 0;
  }


  
  friend ostream& operator<<(ostream& os, const grid& gr);

 private:
  vector<bus> _buses;
  vector<branch> _branches;
  vector<gen> _gens;
  vector<gentype> _gt;

  static map<int, int> map_busNum;

  bool mapBus;

  double loadShedPenalty;

};

#endif

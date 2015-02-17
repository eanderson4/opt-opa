#include "grid.h"

map<int, int> grid::map_busNum;



double grid::getTotalDemand(){
  double d=0;
  for(int i=0;i<numBuses();i++){
    d=d+_buses[i].getP();
  }
  return d;
}

void grid::addGenCost(int num, int model, double startup, double shutdown, int ncost, double c2, double c1, double c0){
    _gens[num-1].addCost(model,startup,shutdown,ncost,c2,c1,c0);
}

int grid::getBusNum(int num){    
    if(!mapBus) buildMap();    
    return map_busNum[num];
}


void grid::buildMap(){
    for(int i=0;i<numBuses();i++){
      map_busNum[_buses[i].getNum()] = i;
    }
    mapBus=true;
}

void grid::printNums(ostream& stream){
    stream<<"Buses: "<<_buses.size()<<endl;
    stream<<"Branches: "<<_branches.size()<<endl;
    stream<<"Generators: "<<_gens.size()<<endl;
}
ostream& operator<<(ostream& os, const grid& gr)
{

  int nB=gr._buses.size();
  int nR=gr._branches.size();
  int nG=gr._gens.size();
  for(int i =0; i<nB; i++){
    os<<gr._buses[i]<<endl;
  }
  for(int i =0; i<nR; i++){
    os<<gr._branches[i]<<endl;
  }
  for(int i =0; i<nG; i++){
    os<<gr._gens[i]<<endl;
  }
  return os;
}

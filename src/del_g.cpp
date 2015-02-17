#include "del_g.h"

void del_g::addDemand(int num, double demand){
  del_demand[num]=demand;
}
    
    

void del_g::setStatus(int num, bool status){ 
  if(!have_topo){
    cerr<<"Set base topology before allowing changes to it"<<endl;
    return;
  }
  del_topo[num]=status; 
}
void del_g::baseTopo(grid * gr){
  int nB = gr->numBranches();
  
  for(int i=0; i<nB;i++){
    del_topo.push_back(true);
  }

  have_topo=true;
}


ostream& operator<<(ostream& os, const del_g& dg)
{
  for(int i=0; i<dg._gr->numBuses(); i++){
      if(dg.del_demand[i]!=0){
	os<<i<<": "<<dg.del_demand[i]<<endl;
      }
  }
  for(int j=0; j<dg._gr->numBranches(); j++){
    if(!dg.del_topo[j]){
      os<<"Line "<<j<<" is out"<<endl;
    }
  }

  return os;
}

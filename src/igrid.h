#ifndef IGRID_H
#define IGRID_H

#include "ibase.h"
#include "rgrid.h"
#include "del_g.h"

class igrid : public ibase {

 public:

 igrid(grid * gr ) : ibase ( ), have_cost(false),have_loadshed(false),have_slack(false){
    _gr = gr;
    if(buildModel( gr  )) cout<<"IloModel Built\n";
    else cout<<"Could not build IloModel\n";
  }
  virtual ~igrid() {  }

  virtual rgrid * solveModel( isolve * is=NULL);
  int addCost();
  int addCost(IloNumVarArray beta,double sigD);
  void addSlack(IloNumArray _g_nom, IloNumArray slack);
  islack * getSlack(){ return &_islk; }
  void slackMismatch(){ _islk.fixMismatch(_gr,getG()); }
  void allowLoadShed();
  void modGrid( del_g mod );
  void removeDemand( del_g mod );
  void unmodGrid( del_g mod );

  double getCost(IloNumArray g){ return _ic.getCost(_gr,g);  }
  grid * getGrid(){ return _gr; }
  icost getIcost(){ return _ic; }
  ished getIshed(){ return _ils; }
  
  bool isLoadShed(){ return have_loadshed; }
  void relayInfo(rgrid * rg);  

 private:
  
  grid * _gr;

  icost _ic;
  islack _islk;
  ished _ils;
  
  del_g _mod;

  bool have_cost;
  bool have_loadshed;
  bool have_slack;

};


#endif

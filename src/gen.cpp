#include "gen.h"

ostream& operator<<(ostream& os, const gen& g)
{
  os << g._num <<"("<<g._bus<<"): "<<g._pg<<" in ["<<g._pmin<<", "<<g._pmax<<"]\t"<<g._status<<"\t"<<g._c2<<"x^2 + "<<g._c1<<"x + "<<g._c0;
  return os;
}

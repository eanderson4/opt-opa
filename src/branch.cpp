#include "branch.h"


ostream& operator<<(ostream& os, const branch& b)
{
  os << b._num <<": "<<b._fbus<<" -> "<<b._tbus<<", "<<b._br_x<<", "<<b._rate_a<<"/"<<b._rate_b<<"/"<<b._rate_c;
  return os;
}

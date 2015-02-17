#include "bus.h"


ostream& operator<<(ostream& os, const bus& b)
{
  os << b._num <<": "<<b._type<<" - "<<b._pd<<", "<<b._gs;
  return os;
}


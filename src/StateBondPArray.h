#ifndef FINITEELEMENT_STATEBONDPARRAY_H
#define FINITEELEMENT_STATEBONDPARRAY_H

//#include <Pointers/NodeP.h>
#include <StateBondP.h>
#include <vector>

namespace FiniteElement {
  
  typedef std::vector<StateBondP> StateBondPArray;
  typedef std::vector<StateBondP>::iterator StateBondPIterator;
  typedef std::vector<StateBondP>::const_iterator constStateBondPIterator;
}

#endif

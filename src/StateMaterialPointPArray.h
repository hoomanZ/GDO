#ifndef FINITEELEMENT_STATEMATERIALPOINTPARRAY_H
#define FINITEELEMENT_SATEMATERIALPOINTPARRAY_H

//#include <Pointers/NodeP.h>
#include <StateMaterialPointP.h>
#include <vector>

namespace FiniteElement {
  
  typedef std::vector<StateMaterialPointP> StateMaterialPointPArray;
  typedef std::vector<StateMaterialPointP>::iterator StateMaterialPointPIterator;
  typedef std::vector<StateMaterialPointP>::const_iterator constStateMaterialPointPIterator;
}

#endif

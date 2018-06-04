#ifndef FINITEELEMENT_PERIBONDPARRAY_H
#define FINITEELEMENT_PERIBONDPARRAY_H

//#include <Pointers/NodeP.h>
#include <PeriBondP.h>
#include <vector>

namespace FiniteElement {
  
  typedef std::vector<PeriBondP> PeriBondPArray;
  typedef std::vector<PeriBondP>::iterator PeriBondPIterator;
  typedef std::vector<PeriBondP>::const_iterator constPeriBondPIterator;
}

#endif

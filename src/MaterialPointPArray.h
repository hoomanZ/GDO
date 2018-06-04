#ifndef FINITEELEMENT_MATERIALPOINTPARRAY_H
#define FINITEELEMENT_MATERIALPOINTPARRAY_H

//#include <Pointers/NodeP.h>
#include <MaterialPointP.h>
#include <vector>

namespace FiniteElement {
  
  typedef std::vector<MaterialPointP> MaterialPointPArray;
  typedef std::vector<MaterialPointP>::iterator MaterialPointPIterator;
  typedef std::vector<MaterialPointP>::const_iterator constMaterialPointPIterator;
}

#endif

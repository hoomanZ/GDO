#ifndef FINITEELEMENT_PERIMATERIALPOINTPARRAY_H
#define FINITEELEMENT_PERIMATERIALPOINTPARRAY_H

//#include <Pointers/NodeP.h>
#include <PeriMaterialPointP.h>
#include <vector>

namespace FiniteElement {
  
  typedef std::vector<PeriMaterialPointP> PeriMaterialPointPArray;
  typedef std::vector<PeriMaterialPointP>::iterator PeriMaterialPointPIterator;
  typedef std::vector<PeriMaterialPointP>::const_iterator constPeriMaterialPointPIterator;
}

#endif

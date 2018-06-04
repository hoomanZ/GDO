#ifndef FINITEELEMENT_PERI_CRACK_H
#define FINITEELEMENT_PERI_CRACK_H

#include <GeometryMath/Vector3D.h>
#include <PeriMaterialPointPArray.h>
#include <PeriBondP.h>

#include <StateMaterialPoint.h> 
#include <StateMaterialPointP.h>  

namespace FiniteElement {

  class PeriCrack {
    public:


      PeriCrack();
      PeriCrack(Vector3D crackMinPoint, Vector3D crackMaxPoint); 

      inline void setMinPoint(const Vector3D& minPoint) {d_min_point = minPoint;}
      inline Vector3D getMinPoint() const {return d_min_point;}

      inline void setMaxPoint(const Vector3D& maxPoint) {d_max_point = maxPoint;}
      inline Vector3D getMaxPoint() const {return d_max_point;}

      








      void createInitialCrack(PeriMaterialPointPArray& pointsArray);
      void createInitialCrack(std::vector<StateMaterialPointP>& statesArray);


      bool intersectWithInitialCrack(const PeriBondP& bond);
      bool intersectWithInitialCrack(const StateBondP& bond);
      bool interval(double Min, double Max, double A, double B);



    private:
      Vector3D d_min_point;
      Vector3D d_max_point;





  };
}


#endif 

#ifndef INTERPOLATEFUNCTION_H_
#define INTERPOLATEFUNCTION_H_

#include <NodeP.h>
#include <GeometryMath/Vector3D>

namespace FiniteElement {
  
  class InterpolateFunction {

    public:
      InterpolateFunction();
      InterpolateFunction(NodeP node);

      inline void setID(const int& id) {d_id = id;}
      inline const int& getID() const {return d_id;}

      inline void xi1(const double& xi_1) {d_xi_one = xi_1;}
      inline const double& xi1() const {return d_xi_one;}

      inline void xi2(const double& xi_2) {d_xi_two = xi_2;}
      inline const double& xi2() const {return d_xi_two;} 

      inline void xi3(const double& xi_3) {d_xi_three = xi_3;}
      inline const double& xi3() const {return d_xi_three;} 

      inline void localNode(const Vector3D& vec) {d_local_node = vec;}
      inline const double& localNode() const {return d_local_node;} 

      inline void setNode(const NodeP& node) {d_node = node;}
      inline const NodeP& getNode() const {return d_node;}  
       
       
       
            
    private:
      NodeP d_node;
      int d_id;
      double d_xi_one;
      double d_xi_two;
      double d_xi_three;
      Vector3D d_local_node;






  } //end class
} //end namespace
#endif

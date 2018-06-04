#ifndef FINITEELEMENT_TRILINEARVOLUMEELEMENT_H_
#define FINITEELEMENT_TRILINEARVOLUMEELEMENT_H_

#include <iostream>
#include <Tensor.h>
#include <Node.h>
#include <NodeP.h>
#include <GeometryMath/Matrix3D.h>
#include <GeometryMath/Vector3D.h>
#include <GeometryMath/Types.h>
#include <NodePArray.h>
#include <array>

#include <TrilinearVolumeElementSP.h>

namespace FiniteElement {

  class TrilinearVolumeElement {

   public:
     TrilinearVolumeElement();
     TrilinearVolumeElement(NodeP node1, NodeP node2,
                           NodeP node3, NodeP node4,
                           NodeP node5, NodeP node6,
                           NodeP node7, NodeP node8);     
//     TrilinearVolumeElement(NodePArray nodeList);
     TrilinearVolumeElement(std::array <NodeP, 8> nodes);
     TrilinearVolumeElement(std::vector <NodeP> nodes,
                            int n1, int n2, int n3, int n4,
                            int n5, int n6, int n7, int n8);     
     ~TrilinearVolumeElement();

     void setArray(NodeP node1, NodeP node2,
                   NodeP node3, NodeP node4,
                   NodeP node5, NodeP node6,
                   NodeP node7, NodeP node8);

     void setNodeP1 (const NodeP& node1) {d_node_1 = node1;}
     NodeP getNodeP1() const {return d_node_1;}

     void setNodeP2 (const NodeP& node2) {d_node_2 = node2;}
     NodeP getNodeP2() const {return d_node_2;}

     void setNodeP3 (const NodeP& node3) {d_node_3 = node3;}
     NodeP getNodeP3() const {return d_node_3;}

     void setNodeP4 (const NodeP& node4) {d_node_4 = node4;}
     NodeP getNodeP4() const {return d_node_4;}

     void setNodeP5 (const NodeP& node5) {d_node_5 = node5;}
     NodeP getNodeP5() const {return d_node_5;}

     void setNodeP6 (const NodeP& node6) {d_node_6 = node6;}
     NodeP getNodeP6() const {return d_node_6;}

     void setNodeP7 (const NodeP& node7) {d_node_7 = node7;}
     NodeP getNodeP7() const {return d_node_7;}

     void setNodeP8 (const NodeP& node8) {d_node_8 = node8;}
     NodeP getNodeP8() const {return d_node_8;}

     
     void setArray(const std::array <NodeP, 8>& nodesArray);
     std::array <NodeP, 8>  getArray() const
     { return d_nodes_array;}

     void operator = (const TrilinearVolumeElement& triElement);
     void operator = (const TrilinearVolumeElementSP& triElement);


     NodeP nodeFunction(std::vector <NodeP> nodes, int i);

     double testFunction(int a, double xi1, double xi2, double xi3); // a=1,...,8     interval= [0, 1]
                                

     double derivativeTestFunctionToXi(int a, int j,
                                       double xi1, double xi2, double xi3); // a=1,...,8 and j=1,2,3

     Matrix3D jacobianMatrix(double xi1, double xi2, double xi3); //derivativePositionToXi

     double derivativeTestFunctionToPosition(int a, int i,
                                             double xi1, double xi2, double xi3);  
 


                                       
     double symTestFunction(int a, double xi1, double xi2, double xi3); // a=1,...,8     interval= [-1, 1]
                                

     double symDerivativeTestFunctionToXi(int a, int j,
                                       double xi1, double xi2, double xi3); // a=1,...,8 and j=1,2,3     

     Matrix3D symJacobianMatrix(double xi1, double xi2, double xi3); //derivativePositionToXi

     double symDerivativeTestFunctionToPosition(int a, int i,
                                             double xi1, double xi2, double xi3); 


     double symXtoXi1(const double& x);
     double symYtoXi2(const double& y);
     double symZtoXi3(const double& z);

     double symXi1toX(const double& xi1);
     double symXi2toY(const double& xi2);
     double symXi3toZ(const double& xi3); 
 
     
     double symInterpolatedValue(const std::array<double, 8>& nodesValues, const double& xi1,
                                                       const double& xi2, const double& xi3);

     Vector3D symInterpolatedVector(const std::array<Vector3D, 8>& nodesVectors, const double& xi1,
                                                              const double& xi2, const double& xi3);

 
     Vector3D symNormVecXi3Max(const double& xi1, const double& xi2);
     Vector3D symNormVecXi3Min(const double& xi1, const double& xi2);
     Vector3D symNormVecXi2Max(const double& xi1, const double& xi3);
     Vector3D symNormVecXi2Min(const double& xi1, const double& xi3);
     Vector3D symNormVecXi1Max(const double& xi2, const double& xi3);
     Vector3D symNormVecXi1Min(const double& xi2, const double& xi3);


     std::array<double, 8> densityArray();



     bool face(const int& i, const int& j,
               const int& k, const int& l);
/************************************************************ 
 * xi3Max, plane xi3 = 1   means    face(5, 6, 7, 8) = true *
 * xi3Min, plane xi3 = -1  means    face(1, 2, 3, 4) = true *
 * xi2Max, plane xi2 = 1   means    face(3, 4, 5, 8) = true *
 * xi2Min, plane xi2 = -1  means    face(1, 2, 6, 7) = true *
 * xi1Max, plane xi1 = 1   means    face(2, 3, 7, 8) = true *
 * xi1Min, plane xi1 = -1  means    face(1, 4, 5, 6) = true *
 ************************************************************/

    std::array<Vector3D, 8> surfaceForceArray(const int& i, const int& j,
                                              const int& k, const int& l);

    std::array<Vector3D, 8> bodyForceArray();
    std::array<Vector3D, 8> oldDisplacementArray();
    std::array<Vector3D, 8> displacementArray();
    std::array<Vector3D, 8> newDisplacementArray();


    double twoDJacobian(const int& i, const int& j,
                        const int& k, const int& l,
                        const double& xi1, const double& xi2);
                                                              




   private:
     Vector3D noTractionForce;
     std::array <NodeP, 8>  d_nodes_array;
     NodeP d_node_1;     //  xi1 = -1,  xi2 = -1,  xi3 = -1             (-1, -1, -1)
     NodeP d_node_2;     //  xi1 = 1,   xi2 = -1,  xi3 = -1             (1, -1, -1)
     NodeP d_node_3;     //  xi1 = 1,   xi2 = 1,   xi3 = -1             (1, 1, -1)
     NodeP d_node_4;     //  xi1 = -1,  xi2 = 1,   xi3 = -1             (-1, 1, -1)
     NodeP d_node_5;     //  xi1 = -1,  xi2 = 1,   xi3 = 1              (-1, 1, 1)
     NodeP d_node_6;     //  xi1 = -1,  xi2 = -1,  xi3 = 1              (-1, -1, 1)
     NodeP d_node_7;     //  xi1 = 1,   xi2 = -1,  xi3 = 1              (1, -1, 1)
     NodeP d_node_8;     //  xi1 = 1,   xi2 = 1,   xi3 = 1              (1, 1, 1)
 

  }; // end of class
} // end of namespace
#endif

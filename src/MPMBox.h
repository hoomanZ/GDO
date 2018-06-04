#ifndef FINITEELEMENT_MPMBOX_H
#define FINITEELEMENT_MPMBOX_H


#include <Body.h>
#include <GeometryMath/Vector3D.h>
#include <GeometryMath/Matrix3D.h>
#include <Box.h>
#include <BoxSP.h>
#include <ComplicatedGeometry.h>
#include <ComplicatedGeometrySP.h>
#include <string>

namespace FiniteElement {

  class MPMBox {
    public:

      enum basisFunction {LinearHexahedral, PiecewiseLinear, GIMPfunction, BSpline};

      MPMBox();
      MPMBox(BoxSP box);
      MPMBox(BoxSP box, basisFunction basis);
      MPMBox(ComplicatedGeometrySP geometry, basisFunction basis);

      double interpolate(NodeP node, MaterialPointP material);
      Vector3D gradientInterpolate(NodeP node, MaterialPointP material);


//***************************************************** Linear Hexahedral basis function *********************************
      bool isPointInElement(MaterialPointP material, SymmetricMaterialTrilinearElementSP element);
      double linearHexahedral(NodeP node, MaterialPointP material);
      Vector3D gradientLinearHexahedral(NodeP node, MaterialPointP material);
//************************************************************************************************************************

//****************************************************** Piecewise-linear basis function *********************************
      double tentFunction(NodeP node, MaterialPointP material);
      double tent(double diff, double length);
      Vector3D grtadientTentFunction(NodeP node, MaterialPointP material);
      double gradientTent(double diff, double length); 
//************************************************************************************************************************

//******************************************************* GIMP ***********************************************************
//      double tentContiguousParticlesGIMPFunction(NodeP node, MaterialPointP material, double particleHalfSize);
      double tentContiguousParticlesGIMPFunction(NodeP node, MaterialPointP material);
      double GIMP(double diff, double length, double lp);
//      Vector3D gradientTentContiguousParticlesGIMPFunction(NodeP node, MaterialPointP material, double particleHalfSize);
      Vector3D gradientTentContiguousParticlesGIMPFunction(NodeP node, MaterialPointP material);
      double gradientGIMP(double diff, double length, double lp);
//************************************************************************************************************************


//******************************************************* B-spline *******************************************************
      double bSplineFunction(NodeP node, MaterialPointP material);
      double bSpline(double diff, double length);
      Vector3D gradientBSplineFunction(NodeP node, MaterialPointP material);
      double gradientBSpline(double diff, double length);
//************************************************************************************************************************

      void updatePointsVolume();
      Vector3D internalForce(NodeP node);
      
      
      double consistentMassComponent(NodeP node1, NodeP node2);
      double lumpedMassComponent(NodeP node);
      Vector3D momentumComponent(NodeP node);
      Vector3D newMomentumComponent(NodeP node);
      Vector3D bodyForceComponent(NodeP node);
      Vector3D pointForceComponent(NodeP node);

      double matrixElement(NodeP node, MaterialPointP material, int t, int s);
      Matrix3D forceIncrementComponent(NodeP node);      

//      Vector3D accelerationComponent(NodeP node);
//      Vector3D velocityComponent(NodeP node);

      inline void setBox(const BoxSP& box) {d_box = box;}
      inline BoxSP getBox() const {return d_box;}

      inline void setGeometry(const ComplicatedGeometrySP& geometry) {d_geometry = geometry;}
      inline ComplicatedGeometrySP getGeometry() const {return d_geometry;}

      inline void setBoxFlag(const bool& boxFlag) {d_box_flag = boxFlag;}
      inline bool getBoxFlag() const {return d_box_flag;}




      BoxSP Boxx();
      std::vector<MaterialPointP> getPoints();


 

    private:
      BoxSP d_box;
      basisFunction d_basis_function;

      ComplicatedGeometrySP d_geometry;

      bool d_box_flag;
      

  };
}
#endif

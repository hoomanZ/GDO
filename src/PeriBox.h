#ifndef FINITEELEMENT_PERI_BOX_H
#define FINITEELEMENT_PERI_BOX_H

#include <Box.h>
//#include <PeriMaterialPointPArray.h>
#include <PeriPMBMaterialP.h>
#include <StateMaterialPointP.h>

//#include <>

namespace FiniteElement {

  class PeriBox : public Box {
    public:

      enum PeriType {Bond, State};
      
      PeriBox();
      PeriBox(Vector3D minPoint, Vector3D maxPoint, int xNumElements, int yNumElements, int zNumElements,
              double horizon);
      PeriBox(Vector3D minPoint, Vector3D maxPoint, int xNumElements, int yNumElements, int zNumElements,
                                           double density, double poission, double young, double horizon);
/*      PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint, 
              const int& xNumElements, const int& yNumElements, const int& zNumElements,
              const double& density, const double& poission, const double& young, const double& bodyForce, 
              const double& surForce, const ForceSurface& surface,
              const FixedDirichletBoundary& boundary, const Vector3D& boundary_position, double horizon,
              PeriPMBMaterialP pmb_material);*/

      PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint,
              const int& xNumElements, const int& yNumElements, const int& zNumElements,
              const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
              const double& density, const double& poission, const double& young,
              const double& bodyForce, const Vector3D& pointForce,
              const double& surForce, const ForceSurface& tracSurface, const ForceSurface& pointSurface,
              const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition, double horizon,
              PeriPMBMaterialP pmb_material);


//********************************************** Constructor for Peridynamics **************************************** 

      PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint,
          const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
          const double& density, const double& poission, const double& young,
          const double& bodyForce, 
          const double& surForce, const ForceSurface& tracSurface,
          const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition, double horizon,
          PeriPMBMaterialP pmb_material);

//*********************************************************************************************************************

      PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint,
              const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
              const double& density, const double& poission, const double& young,
              const double& bodyForce, 
              const double& surForce, const ForceSurface& tracSurface,
              const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition, double horizon,
              PeriPMBMaterialP pmb_material, const ForceSurface& VelocitySurface, 
              const Vector3D& BoundaryVelocity);


      PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint,
              const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
              const double& density, const double& poission, const double& young,
              const double& bodyForce, 
              const double& surForce, const ForceSurface& tracSurface,
              const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition, double horizon,
              PeriPMBMaterialP pmb_material, const PeriType& type, const ForceSurface& VelocitySurface, 
              const Vector3D& BoundaryVelocity);


      PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint,
              const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
              const double& density, const double& poission, const double& young,
              const double& yield, const double& bodyForce, 
              const double& surForce, const ForceSurface& tracSurface,
              const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition, double horizon,
              const PeriType& type, const ForceSurface& VelocitySurface, 
              const Vector3D& BoundaryVelocity);

//********************************************* State-based Peridynamics ***************************************

      PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint,
              const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
              const double& density, const double& poission, const double& young,
              const double& yield, const double& bodyForce, 
              const double& surForce, const ForceSurface& tracSurface,
              const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition, double horizon,
              const PeriType& type, const ForceSurface& VelocitySurface, 
              const Vector3D& BoundaryVelocity, const Matrix3D& InitialVelocityGradient);

//***************************************************************************************************************

      ~PeriBox();
  
      




//      inline void setPeriPointsArray(const std::vector<PeriMaterialPointP>& pointsArray)
//                                                                      {d_peri_points_array = pointsArray;}
//      inline std::vector<PeriMaterialPointP> getPeriPointsArray() const {return d_peri_points_array;}


      inline void setHorizon(const double& horizon) {d_horizon = horizon;}
      inline double getHorizon() const {return d_horizon;}

      inline void setAreaX(const double& area) {d_areaX = area;}
      inline double getAreaX() const {return d_areaX;}

      inline void setAreaY(const double& area) {d_areaY = area;}
      inline double getAreaY() const {return d_areaY;}

      inline void setAreaZ(const double& area) {d_areaZ = area;}
      inline double getAreaZ() const {return d_areaZ;}



      inline void setPMB(const PeriPMBMaterialP& pmb) {d_pmb_material = pmb;}
      inline PeriPMBMaterialP getPMB() const {return d_pmb_material;}


      inline void setBoundaryVelocity(const Vector3D& boundaryVelocity) {d_boundary_velocity = boundaryVelocity;}
      inline Vector3D getBoundaryVelocity() const {return d_boundary_velocity;}


      inline void setVelocityBoundarySurface(const ForceSurface& surface) {d_velocity_surface = surface;}
      inline ForceSurface getVelocityBoundarySurface() const {return d_velocity_surface;}

      inline void setPeriType(const PeriType& type) {d_peri_type = type;}
      inline PeriType getPeriType() const {return d_peri_type;}


      inline void setZeroEnergyConstantTwo(const double& constant) {d_zero_energy_constant_two = constant;}
      inline double getZeroEnergyConstantTwo() const {return d_zero_energy_constant_two;}

      inline void setPmbConstantSpring(const double& constant) {d_pmb_constant_spring = constant;}
      inline double getPmbConstantSpring() const {return d_pmb_constant_spring;}


      inline void setPointsPositiveField(const std::vector <StateMaterialPointP>& field) {d_points_positive_field = field;}
      inline std::vector <StateMaterialPointP> getPointsPositiveField() const {return d_points_positive_field;}


      inline void setPointsNegativeField(const std::vector <StateMaterialPointP>& field) {d_points_negative_field = field;}
      inline std::vector <StateMaterialPointP> getPointsNegativeField() const {return d_points_negative_field;}




//      createPeriMaterialPointArray();


      void findPointsMassAndVolume(); 
      void findPointsAreaBond();
      void findPointsAreaState();
      void calculatePointsBodyForceBond();
      void calculatePointsBodyForceState();
      void calculateBondsSurfaceCorrectionFactor();

      
      void findZeroEnergyConstantTwo();
      

      void calculateMicromodulusFunction();  //Discretized bond-based peridynamics for solid mechanics, W. Liu
      void calculateBondsVolumeBond();
      void calculateBondsVolumeState();

      void applyBoundaryVelocity();
      void applyBoundaryVelocityState();

      void fillBoundaryPointsFields();
      void velocityBoundaryCondition();
      void showBoundaryVelocityPoints();

      void calcualtePmbConstantSpring();


    private:

//       PeriMaterialPointPArray d_peri_points_array;
      double d_horizon;

      double d_areaX;
      double d_areaY;
      double d_areaZ;     

      PeriPMBMaterialP d_pmb_material;

      Vector3D d_boundary_velocity;
      ForceSurface d_velocity_surface;

      PeriType d_peri_type;

      double d_zero_energy_constant_two;  // c2 = (5*10^(13)*E)/(delX*delX)  Quasi-Static Non-Ordinary State-Based Peridynamics for the modeling of 3D fracture. M.S.Britenfeld

      double d_pmb_constant_spring;


      std::vector <StateMaterialPointP>  d_points_positive_field;
      std::vector <StateMaterialPointP>  d_points_negative_field;



  };
}

#endif

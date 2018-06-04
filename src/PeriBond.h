#ifndef FINITEELEMENT_PERIBOND_H__
#define FINITEELEMENT_PERIBOND_H__

#include <PeriMaterialPointP.h>
#include <GeometryMath/Vector3D.h>

namespace FiniteElement {
  class PeriBond {
    public:
      PeriBond();
      PeriBond(PeriMaterialPointP center, PeriMaterialPointP second, bool broken);
      
      inline void setCenterPoint(const PeriMaterialPointP& center) {d_center_point = center;}
      inline PeriMaterialPointP getCenterPoint() const {return d_center_point;}


      inline void setSecondPoint(const PeriMaterialPointP& second) {d_second_point = second;}
      inline PeriMaterialPointP getSecondPoint() const {return d_second_point;}

      inline void setPairwiseForce(const Vector3D& force) {d_pairwise_force_function = force;}
      inline Vector3D getPairwiseForce() const {return d_pairwise_force_function;}

      inline void setBrokenBond(const bool& broken) {d_broken_bond = broken;}
      inline bool getBrokenBond() const {return d_broken_bond;}

      inline void setSurfaceCorrectionFactor(const double& factor) {d_surface_correction_factor = factor;}
      inline double getSurfaceCorrectionFactor() const {return d_surface_correction_factor;}


      inline void setVolume(const double& volume) {d_volume = volume;}
      inline double getVolume() const {return d_volume;}

      inline void setInfluenceFunction(const double& influence) {d_influence_function = influence;}
      inline double getInfluenceFunction() const {return d_influence_function;}

      inline void setForceVectorState(const Vector3D& state) {d_force_vector_state = state;}
      inline Vector3D getForceVectorStat() const {return d_force_vector_state;}




      void clearBond();

      double lengthOfSai();// const {return (d_second_point->getPosOld() - d_center_point->getPosOld()).length();}

      
      void calculateVolume(double horizon, double rj, double delX); 
                                          //Discretized bond-based peridynamics for solid mechanics, W. Liu

    private:
      PeriMaterialPointP d_center_point;
      PeriMaterialPointP d_second_point;
      bool d_broken_bond;
      

      Vector3D d_pairwise_force_function;

      double d_surface_correction_factor;

      double d_volume;

      double d_influence_function;    // State_based Peridynamics
      Vector3D d_force_vector_state;    // State_based Peridynamics

  };
}
#endif

#ifndef FINITEELEMENT_STATEBOND_H__
#define FINITEELEMENT_STATEBOND_H__

#include <StateMaterialPointP.h>
#include <GeometryMath/Vector3D.h>
#include <GeometryMath/Matrix3D.h>
#include <GeometryMath/Types.h>


#include <string>


namespace FiniteElement {
  class StateBond { //: public PeriBond {
    public:
      StateBond();
      StateBond(StateMaterialPointP center, StateMaterialPointP second, bool broken);
      StateBond(StateMaterialPointP center, StateMaterialPointP second, const double& horizon, bool broken);
      
      inline void setCenterPoint(const StateMaterialPointP& center) {d_center_point = center;}
      inline StateMaterialPointP getCenterPoint() const {return d_center_point;}


      inline void setSecondPoint(const StateMaterialPointP& second) {d_second_point = second;}
      inline StateMaterialPointP getSecondPoint() const {return d_second_point;}

//      inline void setPairwiseForce(const Vector3D& force) {d_pairwise_force_function = force;}
//      inline Vector3D getPairwiseForce() const {return d_pairwise_force_function;}

      inline void setBrokenBond(const bool& broken) {d_broken_bond = broken;}
      inline bool getBrokenBond() const {return d_broken_bond;}

      inline void setSurfaceCorrectionFactor(const double& factor) {d_surface_correction_factor = factor;}
      inline double getSurfaceCorrectionFactor() const {return d_surface_correction_factor;}


      inline void setVolume(const double& volume) {d_volume = volume;}
      inline double getVolume() const {return d_volume;}

      inline void setInfluenceFunction(const double& influence) {d_influence_function = influence;}
      inline double getInfluenceFunction() const {return d_influence_function;}


      inline void setForceVectorState(const Vector3D& state) {d_force_vector_state = state;}
      inline Vector3D getForceVectorState() const {return d_force_vector_state;}


      inline void setFirstPiolaDivergenceIntegrand(const Vector3D& state) {d_first_piola_divergence_integrand = state;}
      inline Vector3D getFirstPiolaDivergenceIntegrand() const {return d_first_piola_divergence_integrand;}



      inline void setAverageStrainTensor(const Matrix3D& strain) {d_average_strain_tensor = strain;}
      inline Matrix3D getAverageStrainTensor() const {return d_average_strain_tensor;}


      inline void setEquivalentStrain(const double& equivalent) {d_equivalent_strain = equivalent;}
      inline double getEquivalentStrain() const {return d_equivalent_strain;}


      inline void setVolumetricStrain(const double& volume) {d_volumetric_strain = volume;}
      inline double getVolumetricStrain() const {return d_volumetric_strain;}


      inline void setSai(const Vector3D& sai) {d_sai = sai;}
      inline Vector3D getSai() const {return d_sai;}

      inline void setDeformationVector(const Vector3D& vec) {d_deformation_vector = vec;}
      inline Vector3D getDeformationVector() const {return d_deformation_vector;}


      inline void setVelocityVector(const Vector3D& vec) {d_velocity_vector = vec;}
      inline Vector3D getVelocityVector() const {return d_velocity_vector;}


      inline void setEta(const Vector3D& eta) {d_eta = eta;}
      inline Vector3D getEta() const {return d_eta;}


      inline void setPolynomialBasisOrder1(const Vector3D& polynomial) {d_polynomial_basis_order_1 = polynomial;}
      inline Vector3D getPolynomialBasisOrder1() const {return d_polynomial_basis_order_1;}


      inline void setPolynomialBasisOrder2(const std::array <double, 9>& polynomial) {d_polynomial_basis_order_2 = polynomial;}
      inline std::array <double, 9> getPolynomialBasisOrder2() const {return d_polynomial_basis_order_2;}




      void clearBond();

      double lengthOfSai();// const {return (d_second_point->getPosOld() - d_center_point->getPosOld()).length();}

      
      void calculateVolume(double horizon, double rj, double delX); 
                                          //Discretized bond-based peridynamics for solid mechanics, W. Liu

      void averageStrainTensor();
      void equivalentStrain();
      void volumetricStrain();

      void findInfluenceFunction(const double& horizon);
      void normalExponentialInfluenceFunction(const double& horizon);
      void calculateEta();


      void printScreenMatrix3D(const Matrix3D& matrix, const int& id, const std::string& name);


      Matrix3D tensorProduct(const Vector3D& vec1, const Vector3D& vec2);
      Matrix3D saiTensorSai();  // sai tensor product to itself
      Matrix3D deformVecTensorSai();  // Y<sai> tensor product to sai
      Matrix3D velocityVecTensorSai();  // Ydot<sai> tensor product to sai

      Matrix9  polynomialsOrder2TensorPolynomialsOrder2();

      Matrix3_9 deformVecTensorPolynomialsOrder2();

      Vector3D hourglassForceDensity(const double& springConstant, const Matrix3D& deformation);

     void fillPolynomialBasisOrder1();  
     void fillPolynomialBasisOrder2();


    private:
      StateMaterialPointP d_center_point;
      StateMaterialPointP d_second_point;
      bool d_broken_bond;
      

      Vector3D d_pairwise_force_function;

      double d_surface_correction_factor;

      double d_volume;

      double d_influence_function;    
      Vector3D d_force_vector_state; 


      Vector3D d_first_piola_divergence_integrand; 


      Matrix3D d_average_strain_tensor; // (LagrangianStrain(centre) + LagrangianStrain(second))/2 

      double d_equivalent_strain;
      double d_volumetric_strain; 


      Vector3D d_sai;
      Vector3D d_deformation_vector; // Y<sai>
      Vector3D d_velocity_vector;

      Vector3D d_eta;

      Vector3D  d_polynomial_basis_order_1; // [sai_x, sai_y, sai_z]

      std::array <double, 9>  d_polynomial_basis_order_2; // [sai_x, sai_y, sai_z, (sai_x)^2, (sai_y)^2, (sai_z)^2, sai_x*sai_y, sai_x*sai_z, sai_y*sai_z]
         
   

  };
}
#endif

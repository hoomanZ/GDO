#ifndef FINITEELEMENT_STATEMATERIALPOINT_H__
#define FINITEELEMENT_STATEMATERIALPOINT_H__

#include <PeriMaterialPointP.h>
#include <PeriMaterialPoint.h>
#include <GeometryMath/Vector3D.h>
#include <GeometryMath/Point3D.h>
#include <GeometryMath/Matrix3D.h>
#include <GeometryMath/Types.h>

#include <StateMaterialPointP.h>

#include <StateBondPArray.h>

#include <eigen/Eigen/Eigenvalues>
#include <string>

namespace FiniteElement {

  class StateMaterialPoint : public PeriMaterialPoint {

    typedef Vector3D (StateMaterialPoint::*field)();
    typedef Matrix3D (StateMaterialPoint::*tensor)()const;
//    typedef Matrix3D (StateMaterialPoint::*tensor)();


    public:
      StateMaterialPoint();

      StateMaterialPoint(PeriMaterialPointP point);
      StateMaterialPoint(PeriMaterialPointP point, const Matrix3D& InitialVelocityGradient);

      inline void setStateBonds(const StateBondPArray& bonds) {d_state_bonds_array = bonds;}
      inline StateBondPArray getStateBonds() const {return d_state_bonds_array;}


      inline void setShapeTensor(const Matrix3D& shape) {d_shape_tensor = shape;}
      inline Matrix3D getShapeTensor() const {return d_shape_tensor;}

      inline void setShapeTensorPolynomialOrder2(const Matrix9& shape) {d_shape_tensor_polynomial_order_2 = shape;}
      inline Matrix9 getShapeTensorPolynomialOrder2() const {return d_shape_tensor_polynomial_order_2;}


      inline void setDeformationGradient(const Matrix3D& deform) {d_deformation_gradient = deform;}
      inline Matrix3D getDeformationGradient() const {return d_deformation_gradient;}

      inline void setVelocityGradient(const Matrix3D& velocity) {d_velocity_gradient = velocity;}
      inline Matrix3D getVelocityGradient() const {return d_velocity_gradient;}

      inline void setDeformationRateTensor(const Matrix3D& rate) {d_deformation_rate_tensor = rate;}
      inline Matrix3D getDeformationRateTensor() const {return d_deformation_rate_tensor;}
      
      inline void setSpinTensor(const Matrix3D& spin) {d_spin_tensor = spin;}
      inline Matrix3D getSpinTensor() const {return d_spin_tensor;}

      inline void setRotationTensor(const Matrix3D& rotation) {d_rotation_tensor = rotation;}
      inline Matrix3D getRotationTensor() const {return d_rotation_tensor;}
      
      inline void setRightStretch(const Matrix3D& right) {d_right_stretch = right;}
      inline Matrix3D getRightStretch() const {return d_right_stretch;}

      inline void setLeftStretch(const Matrix3D& left) {d_left_stretch = left;}
      inline Matrix3D getLeftStretch() const {return d_left_stretch;}     

      inline void setRateRigidBodyRotation(const Matrix3D& rate) {d_rate_rigid_body_rotation = rate;}
      inline Matrix3D getRateRigidBodyRotation() const {return d_rate_rigid_body_rotation;}          

      inline void setAngularVelocity(const Vector3D& angular) {d_angular_velocity = angular;}
      inline Vector3D getAngularVelocity() const {return d_angular_velocity;}          

      inline void setLeftStretchRate(const Matrix3D& left) {d_left_stretch_rate = left;}
      inline Matrix3D getLeftStretchRate() const {return d_left_stretch_rate;}     

      inline void setUnrotatedDeformationGradient(const Matrix3D& deform) 
                                                 {d_unrotated_defromation_rate = deform;}
      inline Matrix3D getUnrotatedDeformationGradient()
                                                 const {return d_unrotated_defromation_rate;}

      inline void setTotalElasticStrianIncrement(const Matrix3D& strain) 
                                                {d_total_elastic_strian_increment = strain;}
      inline Matrix3D getTotalElasticStrianIncrement() 
                                                const {return d_total_elastic_strian_increment;}     


      inline void setDeviatoricStrianIncrement(const Matrix3D& strain) 
                                              {d_deviatoric_strian_increment = strain;}
      inline Matrix3D getDeviatoricStrianIncrement() 
                                              const {return d_deviatoric_strian_increment;}     


      inline void setUnrotatedCauchyStress(const Matrix3D& stress) 
                                          {d_unrotated_cauchy_stress = stress;}
      inline Matrix3D getUnrotatedCauchyStress()
                                          const {return d_unrotated_cauchy_stress;}


      inline void setRotatedCauchyStress(const Matrix3D& stress) 
                                        {d_rotated_cauchy_stress = stress;}
      inline Matrix3D getRotatedCauchyStress()
                                        const {return d_rotated_cauchy_stress;}

      inline void setCauchyStress(const Matrix3D& stress) 
                                        {d_cauchy_stress = stress;}
      inline Matrix3D getCauchyStress()
                                        const {return d_cauchy_stress;}




      inline void setFirstPiolaKirchhoffStress(const Matrix3D& stress) 
                                              {d_first_Piola_Kirchhoff_stress = stress;}
      inline Matrix3D getFirstPiolaKirchhoffStress()
                                               const {return d_first_Piola_Kirchhoff_stress;}




      inline void setRightCauchyGreenDeformation(const Matrix3D& deformation) 
                                               {d_right_Cauchy_Green_deformation = deformation;}
      inline Matrix3D getRightCauchyGreenDeformation()
                                        const {return d_right_Cauchy_Green_deformation;}



      inline void setLeftCauchyGreenDeformation(const Matrix3D& deformation) 
                                               {d_left_Cauchy_Green_deformation = deformation;}
      inline Matrix3D getLeftCauchyGreenDeformation()
                                        const {return d_left_Cauchy_Green_deformation;}





      inline void setLagrangianStrainTensor(const Matrix3D& strain) 
                                        {d_Lagrangian_strain = strain;}
      inline Matrix3D getLagrangianStrainTensor()
                                        const {return d_Lagrangian_strain;}

      inline void setInfinitesimalStrainTensor(const Matrix3D& strain) 
                                        {d_infinitesimal_strain = strain;}
      inline Matrix3D getInfinitesimalStrainTensor()
                                        const {return d_infinitesimal_strain;}



      inline void setOrderFlag(const int& order) {d_order_flag = order;}
      inline int getOrderFlag() const {return d_order_flag;}


      inline void setInternalForce(const Vector3D& force) {d_internal_force = force;}
      inline Vector3D getInternalForce() const {return d_internal_force;} 



//************************************* Implicit Time Integration ************************

       inline void setAcceleNew(const Vector3D& accele) {d_accele_new = accele;}
       inline Vector3D getAcceleNew() const {return d_accele_new;}
         

       inline void setBroydenMatrix(const Matrix3D& matrix) {d_Broyden_matrix = matrix;}
       inline Matrix3D getBroydenMatrix() const {return d_Broyden_matrix;} 


      inline void setResidualVector(const Vector3D& residual) {d_residual_vector = residual;}
      inline Vector3D getResidualVector() const {return d_residual_vector;} 


      inline void setDispNew(const Vector3D& disp) {d_disp_new = disp;}
      inline Vector3D getDispNew() const {return d_disp_new;} 

      inline void setNumIter(const int& num) {d_num_iter = num;}
      inline int getNumIter() const {return d_num_iter;}

    

//****************************************************************************************



      MatrixXd matrix3dToEigenMatrix(const Matrix3D& matrix);
      Matrix3D eigenMatrixToMatrix3D(const MatrixXd& matrix);

      Matrix3D invMatrix3D(const Matrix3D& matrix);
      

      int permutationTensor(const int& i, const int& j, const int& k);
      Vector3D dualVector(Matrix3D skew);
      Matrix3D dualSkewMatrix(Vector3D vector);


      void shapeTensor(); // K = X*X
      void deformationGradient();  // F = Y*X^(-1)

      Matrix3D derivativeOperatorApproximationUsingFirstOrderPolynomialBasis(field f);

      void shapeTensorPolynomialOrder2();
      Matrix3D derivativeOperatorApproximationUsingSecondOrderPolynomialBasis(field f);

      void deformationGradientUsingPolynomials();

         

      Matrix3D materialTimeDerivative(); // Fdot = Ydot*X^(-1)
      Matrix3D materialTimeDerivativeUsingPolynomials(); // Fdot = Ydot*X^(-1)

      void spatialVelocityGradient(); // L = (Fdot)(Finverce)
      void spatialVelocityGradient(const int& numIter,
                                   const Matrix3D& initialVelocityGradient); // L = (Fdot)(Finverce)
      void deformationRateTensor(); // D = sym(L) = (L + Ltranspose)/2 
      void spinTensor(); // W = skew(L) = (L - Ltranspose)/2          "Antisymmetric tensor"
      void polarDecomposition(const Matrix3D& deformation,
                                    Matrix3D& rotation,
                                    Matrix3D& right,
                                    Matrix3D& left);
      void polarDecompositionOfDeformationGradient(); 
      void angularVelocity();
      void rateRigidBodyRotation();

      void updateRotationTensorMoreEfficient(const double& delT);  // Hughes and Winget (1980) 
      void updateRotationTensorMoreStable(const double& delT);  // Rubinstein and Atluri (1983)
      void leftStretchRateTensor(); // Vdot = (L)(V) - (V)(omega) 
      void updateLeftStretch(const double& delT);  // V(t+delT) = V(t) + delT*Vdot
      void unrotatedDeformationRate(); // d = (Rtranspose)(D)(R)
      void strainIncrement(const double& delT); // delE = d(delT), delEdev = delE - (1/3)(trace(delE))(I)
      void updateUnrotatedCauchyStress(const double& firstLamme,
                                       const double& secondLamme); 
                           // taw(t+delT) = taw(t) + (firstLamme)(trace(delE))I + 2(secondLamme)(delEdev)
      void rotatedCauchyStress(); // sigma = (R)(taw)(Rtranspose)
      void firstPiolaKirchhoffStress(); // P = (J)(sigma)(FminusTtranspose)
      void forceVectorStateOfBonds();  // T<sai> = (influenceFunction<sai>)(P)(Kinverse)(sai)
      void forceVectorStateOfBonds(const double& constant);  // T<sai> = (influenceFunction<sai>)(P)(Kinverse)(sai)


//************************************* Implicit Time Integration ************************

      void updateAccelerationDisplacementImplicitly(const double& firstLamme,
                                                    const double& secondLamme,
                                                    const double& alpha,
                                                    const double& beta,
                                                    const double& density,
                                                    const double& delT,
                                                    const double& tolerance,
                                                    const int& governingType);


      void updateAccelerationDisplacementImplicitly(const double& firstLamme,
                                                    const double& secondLamme,
                                                    const double& beta,
                                                    const double& alpha,
                                                    const double& density,
                                                    const double& delT,
                                                    const double& tolerance,
                                                    const int& governingType,
                                                    const int& num_iter,
                                                    const std::string& currentFolder);




      Vector3D internalForceLinearElasticImplicit(const double& firstLamme,
                                                  const double& secondLamme,
                                                  const int& governingType);

      void updateBroydenMatrix(const double& alpha);

      void getInformation(Matrix3D& BroydenCurrent,
                          Matrix3D& BroydenNew,
                          Vector3D& residualVectorCurrent,
                          Vector3D& residualVectorNew);

      void updateNumIter(const int& numIter) {setNumIter(numIter);}


//****************************************************************************************



      void calculateForceVectorStateOfBonds(const double& delT,
                                            const double& firstLamme,
                                            const double& secondLamme,
                                            const double& yield);


      void calculateForceVectorStateOfBonds(const double& delT,
                                            const double& firstLamme,
                                            const double& secondLamme,
                                            const double& yield,
                                            const int& numIter, 
                                            const Matrix3D& initialVelocityGradient);

      void calculateForceVectorStateOfBonds(const double& delT,
                                            const double& firstLamme,
                                            const double& secondLamme,
                                            const double& yield,
                                            const int& numIter, 
                                            const Matrix3D& initialVelocityGradient,
                                            const int& orderFlag);


      void calculateFirstPiolaKirchhoffStress(const double& delT,
                                              const double& firstLamme,
                                              const double& secondLamme,
                                              const double& yield,
                                              const int& numIter, 
                                              const Matrix3D& initialVelocityGradient,
                                              const int& orderFlag);

     



      void calculateForceVectorStateOfBonds(const double& delT,
                                            const double& firstLamme,
                                            const double& secondLamme,
                                            const double& yield,
                                            const int& numIter, 
                                            const Matrix3D& initialVelocityGradient,
                                            const double& constant);

      



      void setForceBondsZero();


      void rightCauchyGreenDeformation();
      void leftCauchyGreenDeformation();


      void LagrangianStrainTensor();
      void infinitesimalStrainTensor();
      void unrotatedDeformationRateTensor(const double& delT);

      void isotropicLinearElasticityCauchyStressTensor(const double& firstLamme, const double& secondLamme);
      void hypoelasticRotatedCauchyStressTensor(const double& firstLamme, const double& secondLamme);
      void firstPiolaKirchhoffStress(const Matrix3D& CauchyStressTensor);


      Vector3D divergenceMatrix3D(tensor f);
      void internalForce();
      void internalForceUsingStateBasedPeridynamics();
      void setBondsInfluenceFunctionsZero();

      Vector3D nonLocalDivergenceMatrix3D(tensor f);

      void currentPositionVelocityAcceleration(const double& density, const double& delT);

      void cauchyStress(const double& firstLamme, const double& secondLamme);
      void updatePosition(const double& density, const double& delT);
      void updatePositionUsingStateBasedPeridynamics(const double& density, const double& delT);

      StateBondP twinBond(const StateBondP& bond);
      


      void calculateStateLocalDamage();
      void checkBrokenBonds(const double& yield);
      void checkBrokenBonds(const double& yield, bool& flag);

      void printScreenMatrix3D(const Matrix3D& matrix);
      void printScreenMatrix3D(const Matrix3D& matrix, const std::string& name); 

      double ramp(const int& num_iter, const int& num_initial);
      double triangle(const int& num_iter, const int& num_initial);

      Vector3D returnCurrentPosition();
      Vector3D returnCurrentVelocity();
      Vector3D returnDisplacement();

      Vector3D returnDisplacementNew();

      Vector3D returnVectorExample();

      Matrix3D returnMatrixExample();

      void derivative(const bool& secondOrder);

//      void divergance(const bool& diverganceMethod, const bool& secondOrder);

      Matrix9 pseudoInverseMoorePenrose9by9(const Matrix9& matrix9);
      Matrix2 pseudoInverseMoorePenrose2by2(const Matrix2& matrix2);
      Matrix3 pseudoInverseMoorePenrose3by3(const Matrix3& matrix3);

      double machine_eps (double value);

   

    private:

      StateBondPArray d_state_bonds_array;
      
      Matrix3D d_shape_tensor;
      Matrix3D d_deformation_gradient;
      Matrix3D d_velocity_gradient;            
      Matrix3D d_deformation_rate_tensor;
      Matrix3D d_spin_tensor;
      Matrix3D d_rotation_tensor;   // F = RU,  R is the rotation tensor
      Matrix3D d_right_stretch;     // U is the right stretch tensor 
      Matrix3D d_left_stretch;      // F = VR,  V is the left stretch tensor           
      Matrix3D d_rate_rigid_body_rotation;  // omega = (Rdot)(Rinverse)
      Matrix3D d_left_stretch_rate;  // Vdot = (L)(V) - (V)(omega)

      Vector3D d_angular_velocity;  // angular veloctiy vector corresponding to omega

      Matrix3D d_unrotated_defromation_rate; // d = (Rtranspose)(D)(R)
      Matrix3D d_total_elastic_strian_increment; // delE = (d)(delT)
      Matrix3D d_deviatoric_strian_increment; // delEdev = delE - (1/3)(trace(delE))(I)
      Matrix3D d_unrotated_cauchy_stress; // taw
      Matrix3D d_rotated_cauchy_stress; // sigma = (R)(taw)(Rtranspose)
      Matrix3D d_first_Piola_Kirchhoff_stress; // P = (J)(sigma)(FminusTranspose)


      Matrix3D d_right_Cauchy_Green_deformation; // C = Ftranspose*F
      Matrix3D d_left_Cauchy_Green_deformation; // B = F*Ftranspose

      Matrix3D d_Lagrangian_strain; // E = 0.5*((Ftranspose)(F) - I)
      Matrix3D d_infinitesimal_strain; // epsilon = 0.5*(F + Ftranspose) - I



      Matrix9 d_shape_tensor_polynomial_order_2 = Matrix9::Zero();

      int d_order_flag;     // the maximum order of the polynomial basis functions


      Matrix3D d_cauchy_stress; //  sigma[ij] = firstLamme*delta[ij]*epsilon[kk] + 2*secondLamme*epsilon[ij]

      Vector3D d_internal_force;


//********************************** Implicit Time Integration *****************************

      Vector3D d_accele_new;
      Matrix3D d_Broyden_matrix;   // Broyden's method  
      Vector3D d_residual_vector;
      Vector3D d_disp_new;

      int d_num_iter;
      

//*******************************************************************************************

  };
}

#endif

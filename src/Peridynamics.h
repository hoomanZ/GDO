#ifndef FINITEELEMENT_PERIDYNAMICS_H__
#define FINITEELEMENT_PERIDYNAMICS_H__

#include <PeriPMBMaterialP.h>
#include <BodySP.h>
#include <PeriBoxP.h>

#include <PeriCrackP.h>
#include <PeriCrack.h>




namespace FiniteElement {

  class Peridynamics {
    public:

      enum BasisPolynomialOrders {NoPolynomial, One, Two};
      enum GoverningEquationType {BondBasedPeridynamics, StateBasedPeridynamics, Divergance};
      enum TimeIntegration {Explicit, Implicit};


//      enum PreCrackPlane
//      {
//        
//      };


      Peridynamics();
//      Peridynamics(PeriPMBMaterialP pmbMaterial);
 //     Peridynamics(PeriPMBMaterialP pmbMaterial, BodySP body);
      Peridynamics(double totalT, double delT);
      Peridynamics(double totalT, double delT, PeriBoxP periBox);
      Peridynamics(double totalT, double delT, PeriBoxP periBox, PeriCrackP periCrack);
      Peridynamics(double totalT, double delT, PeriBoxP periBox, const Matrix3D& initialVelocityGradient);
      Peridynamics(double totalT, double delT, PeriBoxP periBox,
                   const Matrix3D& initialVelocityGradient, PeriCrackP periCrack);

      Peridynamics(double totalT, double delT, PeriBoxP periBox,
                   const Matrix3D& initialVelocityGradient, PeriCrackP periCrack, const BasisPolynomialOrders& order);

      Peridynamics(double totalT, double delT, PeriBoxP periBox,
                   const Matrix3D& initialVelocityGradient, PeriCrackP periCrack,
                   const BasisPolynomialOrders& order, const GoverningEquationType& type);


      Peridynamics(double totalT, double delT, PeriBoxP periBox,
                   PeriCrackP periCrack,
                   const BasisPolynomialOrders& order, const GoverningEquationType& type);


      Peridynamics(double totalT, double delT, PeriBoxP periBox,
                   PeriCrackP periCrack, const BasisPolynomialOrders& order, 
                   const GoverningEquationType& type, const TimeIntegration& timeIntegration,
                   double alpha, double beta, double tolerance);




      ~Peridynamics();


//     inline void setPMBMaterial(const PeriPMBMaterialP& pmbMaterial) {d_pmb_material = pmbMaterial;}
//     inline PeriPMBMaterialP getPMBMaterial() const {return d_pmb_material;}

     inline void setPeriBox(const PeriBoxP& periBox) {d_peri_box = periBox;}
     inline PeriBoxP getPeriBox() const {return d_peri_box;}


     inline void setTimeStep(const double& delT) {d_del_t = delT;}
     inline double getTimeStep() const {return d_del_t;}

     inline void setTotalTime(const double& total) {d_total_t = total;}
     inline double getTotalTime() const {return d_total_t;}

     inline void setPeriCrack(const PeriCrackP& periCrack) {d_peri_crack = periCrack;}
     inline PeriCrackP getPeriCrack() const {return d_peri_crack;}

     inline void setBasisPolynomialOrders(const BasisPolynomialOrders& polynomial) {d_basis_polynomial_orders = polynomial;}
     inline BasisPolynomialOrders getBasisPolynomialOrders() const {return d_basis_polynomial_orders;}

     inline void setGoverningEquationType(const GoverningEquationType& type) {d_governing_equation_type = type;}
     inline GoverningEquationType getGoverningEquationType() const {return d_governing_equation_type;}


//***************************** Implicit Time Integration ***************************

     inline void setTimeIntegration(const TimeIntegration& type) {d_time_integration = type;}
     inline TimeIntegration getTimeIntegration() const {return d_time_integration;}


     inline void setAlpha(const double& alpha) {d_alpha = alpha;}
     inline double getAlpha() const {return d_alpha;}
    
     inline void setBeta(const double& beta) {d_beta = beta;}
     inline double getBeta() const {return d_beta;}

     inline void setTolerance(const double& tolerance) {d_tolerance = tolerance;}
     inline double getTolerance() const {return d_tolerance;}

     void printDetails(int& iter_num, const std::string& currentFolder);


//************************************************************************************

     void createGeometry();
     void createBonds();
     void findStressOfPoints();
     void periPMBBoxVelocityVerlet();
     void run();

     void stateLargeRotation();
     void stateLargeRotation(const int& numIter, const Matrix3D& initialVelocityGradient);
     void stateRun();    
     void stateRun(const Matrix3D& initialVelocityGradient);


     void polynomialOrderPoint();
     void linearElasticRun();
     void damageModelingUsingCriticalStrain(const double& yield);
     void damageModelingUsingCriticalStrain(const double& yield, const int& num_iter);




     void pointsFile(int& iter_num, const std::string& currentFolder);
     void pointsFileState(int& iter_num, const std::string& currentFolder);


     bool intersectBondWithPlane(const PeriBondP& bond, const Vector3D& minPoint, const Vector3D& maxPoint);
     bool interval(double Min, double Max, double A, double B);


     double ramp(const int& num_iter, const int& num_initial);

     double triangle(const int& num_iter, const int& num_initial);


     void testDerivative();

     void testDivergance();


   private:
//       PeriPMBMaterialP d_pmb_material; 
//       BodySP d_body;

       PeriBoxP d_peri_box;

       double d_horizon;
       double d_del_t;
       double d_total_t;

       PeriCrackP d_peri_crack;

       BasisPolynomialOrders d_basis_polynomial_orders;
       GoverningEquationType d_governing_equation_type;
       TimeIntegration d_time_integration;


//************* Implicit Time Integration ****************

       double d_alpha;
       double d_beta;
       double d_tolerance;

//********************************************************

  };
}
#endif

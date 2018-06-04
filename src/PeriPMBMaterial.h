#ifndef FINITEELEMENT_PERIPMBMATERIAL_H
#define FINITEELEMENT_PERIPMBMATERIAL_H


#include <BodySP.h>
#include <PeriMaterialPointP.h>
#include <PeriBondP.h>
#include <GeometryMath/Vector3D.h>
#include <GeometryMath/Point3D.h>


#include <vector>


namespace FiniteElement {
 class PeriPMBMaterial {
   public:

     PeriPMBMaterial();
     PeriPMBMaterial(double fractureEnergy, double horizon);
     PeriPMBMaterial(double density, double horizon, double fractureEnergy, double youngModulus);
     PeriPMBMaterial(double density, double horizon, double youngModulus, bool bondBased);

//     PeriPMBMaterial(double density, double horizon, double fractureEnergy, double youngModulus, const BodySP& body);

     inline void setDensity(const double& density) {d_density = density;}
     inline double getDensity() const {return d_density;}
     
     inline void setHorizon(const double& horizon) {d_horizon = horizon;}
     inline double getHorizon() const {return d_horizon;}

     inline void setFractureEnergy(const double& fractureEnergy) {d_fracture_energy = fractureEnergy;}
     inline double getFractureEnergy() const {return d_fracture_energy;}

     inline void setYoungsModulus(const double& youngsModulus) {d_youngs_modulus = youngsModulus;}
     inline double getYoungsModulus() const {return d_youngs_modulus;}

     inline void setCriticalStretch(const double& criticalStretch) {d_critical_stretch = criticalStretch;}
     inline double getCriticalStretch() const {return d_critical_stretch;}

     inline void setMicromodulusFunction(const double& MicromodulusFunction) 
                                        {d_micromodulus_function = MicromodulusFunction;}
     inline double getMicromodulusFunction() const {return d_micromodulus_function;}

//     inline void setBody(const BodySP& body) {d_body = body;}
//     inline BodySP getBody() const {return d_body;}


     void calculateCriticalStretch();
     void calculateMicromodulusFunction();

     void calculatePairwiseForceOfEachBond(const std::vector<PeriMaterialPointP>& pointsArray);


     Vector3D calculatePairwiseForce(Point3D pos_second, Point3D pos_center,
                                     Vector3D disp_second, Vector3D disp_center, bool broken);
     Vector3D coefficientOfPoint(PeriMaterialPointP point); //g_coeff(x) = [g_inf/g_x, g_inf/g_y, g_inf/g_z]
     void surfaceCorrectionFactor(PeriBondP bond); 


   private:
     double d_density;
     double d_horizon;
     double d_fracture_energy;
     double d_youngs_modulus;

     double d_critical_stretch;
     double d_micromodulus_function;
     
//     BodySP d_body;
  

 };
}
#endif

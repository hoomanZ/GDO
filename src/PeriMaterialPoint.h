#ifndef FINITEELEMENT_PERIMATERIALPOINT_H__
#define FINITEELEMENT_PERIMATERIALPOINT_H__


#include <MaterialPoint.h>
#include <MaterialPointP.h>
#include <PeriBondPArray.h>


namespace FiniteElement {
  
  class PeriMaterialPoint : public MaterialPoint {
     public:
       PeriMaterialPoint();
       PeriMaterialPoint(MaterialPointP point);  
       
       inline void setVelMid(const Vector3D& vel) {d_vel_mid = vel;}
       inline Vector3D getVelMid() const {return d_vel_mid;}

       inline void setDispOld(const Vector3D& disp) {d_disp_old = disp;}
       inline Vector3D getDispOld() const {return d_disp_old;}

       inline void setDispNew(const Vector3D& disp) {d_disp_new = disp;}
       inline Vector3D getDispNew() const {return d_disp_new;}

       inline void setAccele(const Vector3D& accele) {d_accele = accele;}
       inline Vector3D getAccele() const {return d_accele;}
    
       inline void setHorizonBonds(const PeriBondPArray& bonds) {d_horizon_bonds_array = bonds;}
       inline PeriBondPArray getHorizonBonds() const {return d_horizon_bonds_array;}

       inline void setLocalDamage(const double& damage) {d_local_damage = damage;}
       inline double getLocalDamage() const {return d_local_damage;}

       inline void setUniaxialTensionX(const Vector3D& tension) {d_uniaxial_tension_x = tension;}
       inline Vector3D getUniaxialTensionX() const {return d_uniaxial_tension_x;}

       inline void setUniaxialTensionY(const Vector3D& tension) {d_uniaxial_tension_y = tension;}
       inline Vector3D getUniaxialTensionY() const {return d_uniaxial_tension_y;}

       inline void setUniaxialTensionZ(const Vector3D& tension) {d_uniaxial_tension_z = tension;}
       inline Vector3D getUniaxialTensionZ() const {return d_uniaxial_tension_z;}

       inline void setVelocityBoundaryFlag(const bool& velBoundary) {d_velocity_boundary = velBoundary;}
       inline bool getVelocityBoundaryFlag() const {return d_velocity_boundary;}

       inline void setPeriStress(const Vector3D& stress) {d_peri_stress = stress;}
       inline Vector3D getPeriStress() const {return d_peri_stress;}

       inline void setPeriStrain(const Vector3D& strain) {d_peri_strain = strain;}
       inline Vector3D getPeriStrain() const {return d_peri_strain;}

       
       void makePointHorizon(double horizon);
       void calculateAcceleration(double density);
       void calculateLocalDamage();

       void setUniaxialTensions(); // For surface correction factors


    private:
        Vector3D d_vel_mid;

        Vector3D d_disp_old;
        Vector3D d_disp_new;

        Vector3D d_accele;

        PeriBondPArray d_horizon_bonds_array;

        double d_local_damage;

        Vector3D d_uniaxial_tension_x;  // For surface correction factors
        Vector3D d_uniaxial_tension_y;
        Vector3D d_uniaxial_tension_z;

        bool d_velocity_boundary;

        Vector3D d_peri_stress;
        Vector3D d_peri_strain;
  
  };
}
#endif

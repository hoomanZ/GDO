//#ifndef FINITEELEMENT_MATERIALPOINT_H__
//#define FINITEELEMENT_MATERIALPOINT_H__


//#include <GeometryMath/Types.h>
//#include <GeometryMath/Point3D.h>
//#include <GeometryMath/Vector3D.h>
//#include <GeometryMath/Matrix3D.h>


//#include <iostream>
                                                                                                                                                                                                                                                                                                                                         

namespace FiniteElement {  

  class MaterialPoint {
    public:
      MaterialPoint();
      MaterialPoint(const int id, const double xx, const double yy, const double zz);
      ~MaterialPoint();



      inline void setPosOld(const Point3D& pos) {d_pos_old = pos;}
      inline Point3D getPosOld() const {return d_pos_old;}

//      inline void setPosNew(const Point3D& pos) {d_pos_new = pos;}
//      inline Point3D getPosNew() const {return d_pos_new;}

      inline void setDisp(const Vector3D& disp) {d_disp = disp;}
      inline Vector3D getDisp() const {return d_disp;}

      inline void setVelOld(const Vector3D& vel) {d_vel_old = vel;}
      inline Vector3D getVelOld() const {return d_vel_old;}

//      inline void setVelNew(const Vector3D& vel) {d_vel_new = vel;}
//      inline Vector3D getVelNew() const {return d_vel_new;}

      inline void setBodyForce(const Vector3D& force) {d_body_force = force;}
      inline Vector3D getBodyForce() const {return d_body_force;}

//      inline void setPointForce(const Vector3D& force) {d_point_force = force;}
//      inline Vector3D getPointForce() const {return d_point_force;}

//      inline void setDeformationGradientOld(const Matrix3D& deform) {d_deformation_gradient_old = deform;}
//      inline Matrix3D getDeformationGradientOld() const {return d_deformation_gradient_old;}

//      inline void setDeformationGradientNew(const Matrix3D& deform) {d_deformation_gradient_new = deform;}
//      inline Matrix3D getDeformationGradientNew() const {return d_deformation_gradient_new;}

//      inline void setStrainOld(const Matrix3D& strain) {d_strain_old = strain;}
//      inline Matrix3D getStrainOld() const {return d_strain_old;}

//      inline void setStrainNew(const Matrix3D& strain) {d_strain_new = strain;}
//      inline Matrix3D getStrainNew() const {return d_strain_new;}

//      inline void setStressOld(const Matrix3D& stress) {d_stress_old = stress;}
//      inline Matrix3D getStressOld() const {return d_stress_old;}

//      inline void setStressNew(const Matrix3D& stress) {d_stress_new = stress;}
//      inline Matrix3D getStressNew() const {return d_stress_new;}

//      inline void setStrainIncrement(const Matrix3D& strain) {d_strain_increment = strain;}
//      inline Matrix3D getStrainIncrement() const {return d_strain_increment;}

//      inline void setStressIncrement(const Matrix3D& stress) {d_stress_increment = stress;}
//      inline Matrix3D getStressIncrement() const {return d_stress_increment;}

      inline void setID(const int& id) {d_id = id;}
      inline int getID() const {return d_id;}

      inline void setMass(const double& mass) {d_mass = mass;}
      inline double getMass() const {return d_mass;}

      inline void setVolume(const double& volume) {d_volume = volume;}
      inline double getVolume() const {return d_volume;}

      inline void setInitialVolume(const double& volume) {d_initial_volume = volume;}
      inline double getInitialVolume() const {return d_initial_volume;}

//      inline void setNoneZeroExtForce(const bool& noneZero) {d_none_zero_ext_force = noneZero;}
//      inline bool getNoneZeroExtForce() const {return d_none_zero_ext_force;}

//      inline void setHasBoundaryCondition(const bool& has) {d_has_boundary_condition = has;}
//      inline bool getHasBoundaryCondition() const {return d_has_boundary_condition;} 



    private:
      Point3D d_pos_old;
      Point3D d_pos_new;

      Vector3D d_disp;  
 
      Vector3D d_vel_old;
      Vector3D d_vel_new;

      Vector3D d_body_force;
      Vector3D d_point_force;

      Matrix3D d_deformation_gradient_old;
      Matrix3D d_deformation_gradient_new;

      Matrix3D d_strain_old;
      Matrix3D d_strain_new;

      Matrix3D d_stress_old;
      Matrix3D d_stress_new;

      Matrix3D d_strain_increment;
      Matrix3D d_stress_increment;


      int d_id;
      double d_mass;
      double d_initial_volume;
      double d_volume;
     
      bool d_none_zero_ext_force;
      bool d_has_boundary_condition;


  };
}

//#endif

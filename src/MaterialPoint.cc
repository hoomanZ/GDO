#include <MaterialPoint.h>


using namespace FiniteElement;


MaterialPoint::MaterialPoint()
  : d_id(0), d_mass(10.0), d_volume(10.0), d_initial_volume(10.0),
    d_pos_old(0.0, 0.0, 0.0), d_pos_new(0.0, 0.0, 0.0), d_disp(0.0, 0.0, 0.0),
    d_vel_old(0.0, 0.0, 0.0), d_vel_new(0.0, 0.0, 0.0), 
    d_body_force(0.0, 0.0, 0.0), d_point_force(0.0, 0.0, 0.0), d_surface_force(0.0, 0.0, 0.0),
    d_deformation_gradient_old(1.0, 0.0, 0.0,
                               0.0, 1.0, 0.0,
                               0.0, 0.0, 1.0), 
    d_deformation_gradient_new(1.0, 0.0, 0.0,
                               0.0, 1.0, 0.0,
                               0.0, 0.0, 1.0), 
    d_strain_old(0.0), d_strain_new(0.0), d_stress_old(0.0), d_stress_new(0.0), 
    d_strain_increment(0.0), d_stress_increment(0.0),
    d_none_zero_ext_force(false), d_has_boundary_condition(false)
{
  Vector3D vec(0.0, 0.0, 0.0);
  vec.x(getPosNew().x());
  vec.y(getPosNew().y());
  vec.z(getPosNew().z());
  setPositionVector(vec);

}

MaterialPoint::MaterialPoint(const int id, const double xx, const double yy, const double zz)
  : d_id(id), d_mass(10.0), d_volume(10.0), d_initial_volume(10.0),
    d_pos_old(xx, yy, zz), d_pos_new(xx, yy, zz), d_disp(0.0, 0.0, 0.0),
    d_vel_old(0.0, 0.0, 0.0), d_vel_new(0.0, 0.0, 0.0), 
    d_body_force(0.0, 0.0, 0.0), d_point_force(0.0, 0.0, 0.0), d_surface_force(0.0, 0.0, 0.0),
    d_deformation_gradient_old(1.0, 0.0, 0.0,
                               0.0, 1.0, 0.0,
                               0.0, 0.0, 1.0), 
    d_deformation_gradient_new(1.0, 0.0, 0.0,
                               0.0, 1.0, 0.0,
                               0.0, 0.0, 1.0), 
    d_strain_old(0.0), d_strain_new(0.0), d_stress_old(0.0), d_stress_new(0.0), 
    d_strain_increment(0.0), d_stress_increment(0.0),
    d_none_zero_ext_force(false), d_has_boundary_condition(false)
{
  Vector3D vec(0.0, 0.0, 0.0);
  vec.x(getPosNew().x());
  vec.y(getPosNew().y());
  vec.z(getPosNew().z());
  setPositionVector(vec);

}

MaterialPoint::~MaterialPoint()
{
}



Vector3D
MaterialPoint::newPositionVector()
{
  Vector3D vec(0.0, 0.0, 0.0);
  vec.x(getPosOld().x()+getDisp().x());
  vec.y(getPosOld().y()+getDisp().y());
  vec.z(getPosOld().z()+getDisp().z());
  return vec;
}



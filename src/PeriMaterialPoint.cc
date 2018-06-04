#include <PeriMaterialPoint.h>
#include <PeriMaterialPointP.h>
#include <PeriBond.h>
#include <PeriBondP.h>
#include <GeometryMath/Vector3D.h>


using namespace FiniteElement;


PeriMaterialPoint::PeriMaterialPoint()
 : MaterialPoint::MaterialPoint(), d_vel_mid(0.0), d_disp_old(0.0), d_disp_new(0.0),
   d_accele(0.0), d_peri_stress(0.0), d_peri_strain(0.0), d_velocity_boundary(false)
{
  setUniaxialTensions();

}

PeriMaterialPoint::PeriMaterialPoint(MaterialPointP point)
 : d_vel_mid(0.0), d_disp_old(0.0), d_disp_new(0.0),
   d_accele(0.0), d_peri_stress(0.0), d_peri_strain(0.0),
   d_velocity_boundary(false)
{
   d_horizon_bonds_array.reserve(6000);
   MaterialPoint::setPosOld(point->getPosOld()); 
   MaterialPoint::setPosNew(point->getPosNew());
   MaterialPoint::setDisp(point->getDisp());
   MaterialPoint::setVelOld(point->getVelOld());
   MaterialPoint::setVelNew(point->getVelNew());
   MaterialPoint::setBodyForce(point->getBodyForce());
   MaterialPoint::setPointForce(point->getPointForce());
   MaterialPoint::setDeformationGradientOld(point->getDeformationGradientOld());
   MaterialPoint::setDeformationGradientNew(point->getDeformationGradientNew());
   MaterialPoint::setStrainOld(point->getStrainOld());
   MaterialPoint::setStrainNew(point->getStrainNew());
   MaterialPoint::setStressOld(point->getStressOld());
   MaterialPoint::setStressNew(point->getStressNew());
   MaterialPoint::setStrainIncrement(point->getStrainIncrement());
   MaterialPoint::setStressIncrement(point->getStressIncrement());
   MaterialPoint::setID(point->getID());
   MaterialPoint::setMass(point->getMass());
   MaterialPoint::setVolume(point->getVolume());
   MaterialPoint::setInitialVolume(point->getInitialVolume());
   MaterialPoint::setNoneZeroExtForce(point->getNoneZeroExtForce());
   MaterialPoint::setHasBoundaryCondition(point->getHasBoundaryCondition());

   setUniaxialTensions();

}


void
PeriMaterialPoint::calculateAcceleration(double density)
{
  Vector3D internal_force(0.0, 0.0, 0.0);
  PeriBondPArray bondsArray = getHorizonBonds();
  int size = bondsArray.size();
  for (int i = 0; i < size; i++)
  { 
    PeriBondP bond = bondsArray[i];
//    PeriMaterialPointP second_point = bond->getSecondPoint();
    










  }

}



void 
PeriMaterialPoint::calculateLocalDamage()
{
  PeriBondPArray bondsArray = getHorizonBonds();
  int size = bondsArray.size();
  double numerator = 0.0;
  double denominator = 0.0;
  for (int i = 0; i < size; i++)
  {
    PeriBondP bond = bondsArray[i];
    PeriMaterialPointP secondPoint = bond->getSecondPoint();
    double volume = secondPoint->getInitialVolume();
    if (!bond->getBrokenBond())
      numerator += volume;
    denominator += volume;
    
  }
  double damage = numerator/denominator;
  setLocalDamage(damage);
  
}



void 
PeriMaterialPoint::setUniaxialTensions()
{
  Vector3D uniaxial_x(this->x()*(0.001), y()*(-0.00025), z()*(-0.00025));
  setUniaxialTensionX(uniaxial_x);

  Vector3D uniaxial_y(x()*(-0.00025), y()*(0.001), z()*(-0.00025));
  setUniaxialTensionY(uniaxial_y);

  Vector3D uniaxial_z(x()*(-0.00025), y()*(-0.00025), z()*(0.001));
  setUniaxialTensionZ(uniaxial_z);


}

#include <PeriBox.h>
#include <PeriBond.h>
#include <StateBond.h>
#include <PeriPMBMaterial.h>
#include <PeriMaterialPoint.h>
//#include <StateMaterialPointP.h>
#include <MaterialPoint.h>

#include <math.h>

using namespace FiniteElement;


PeriBox::PeriBox()
  : Box(), d_horizon(0.1),  // d_pmb_material(new PeriPMBMaterial()),
    d_areaX(0.0), d_areaY(0.0), d_areaZ(0.0),
    d_boundary_velocity(0.0),
    d_velocity_surface(none), d_peri_type(Bond) 
{
  setPeriPoints();
  createPeriBondsOfEachPeriPoints(getHorizon());
    
  findPointsAreaBond();

  findZeroEnergyConstantTwo();
}


PeriBox::PeriBox(Vector3D minPoint, Vector3D maxPoint, int xNumElements, int yNumElements, int zNumElements,
                 double horizon)
  : Box(minPoint, maxPoint, xNumElements, yNumElements, zNumElements), d_horizon(horizon),
    d_pmb_material(new PeriPMBMaterial()),
    d_areaX(0.0), d_areaY(0.0), d_areaZ(0.0),
    d_boundary_velocity(0.0),
    d_velocity_surface(none), d_peri_type(Bond) 
{
  setPeriPoints();
  createPeriBondsOfEachPeriPoints(getHorizon());  
  findPointsAreaBond();

  findZeroEnergyConstantTwo();

}


PeriBox::PeriBox(Vector3D minPoint, Vector3D maxPoint, int xNumElements, int yNumElements, int zNumElements,
                                              double density, double poission, double young, double horizon)
  : Box(minPoint, maxPoint, xNumElements, yNumElements, zNumElements, density, poission, young),
    d_horizon(horizon),
    d_areaX(0.0), d_areaY(0.0), d_areaZ(0.0),
    d_boundary_velocity(0.0),
    d_velocity_surface(none), d_peri_type(Bond) 
{
  setPeriPoints();
  createPeriBondsOfEachPeriPoints(getHorizon());
  getPMB()->setDensity(density);  
  getPMB()->setHorizon(horizon);  
  findPointsAreaBond();

  findZeroEnergyConstantTwo();

}


/*PeriBox::PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint, 
                 const int& xNumElements, const int& yNumElements, const int& zNumElements,
                 const double& density, const double& poission, const double& young, const double& bodyForce, 
                 const double& surForce, const ForceSurface& surface,
                 const FixedDirichletBoundary& boundary, const Vector3D& boundary_position,
                 double horizon, PeriPMBMaterialP pmb_material)
  : Box(minPoint, maxPoint, xNumElements, yNumElements, zNumElements, density, poission, young, bodyForce,
        surForce, surface, boundary, boundary_position), d_horizon(horizon), d_pmb_material(pmb_material),
    d_areaX(0.0), d_areaY(0.0), d_areaZ(0.0)
{
  setPeriPoints();
  createPeriBondsOfEachPeriPoints(getHorizon());  
  findPointsArea();

}*/


PeriBox::PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint,
                 const int& xNumElements, const int& yNumElements, const int& zNumElements,
                 const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
                 const double& density, const double& poission, const double& young,
                 const double& bodyForce, const Vector3D& pointForce,
                 const double& surForce, const ForceSurface& tracSurface, const ForceSurface& pointSurface,
                 const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition,
                 double horizon, PeriPMBMaterialP pmb_material)
   : Box(minPoint, maxPoint, xNumElements, yNumElements, zNumElements,
         xNumPoints, yNumPoints, zNumPoints, density, poission, young,
         bodyForce, pointForce, surForce, tracSurface, pointSurface,
         boundary, boundaryPosition), 
     d_horizon(horizon), d_pmb_material(pmb_material),
     d_areaX(0.0), d_areaY(0.0), d_areaZ(0.0),
     d_boundary_velocity(0.0),
     d_velocity_surface(none), d_peri_type(Bond) 
{
  setPeriPoints();
  createPeriBondsOfEachPeriPoints(getHorizon());  
  findPointsAreaBond();

}


PeriBox::PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint,
                 const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
                 const double& density, const double& poission, const double& young,
                 const double& bodyForce,
                 const double& surForce, const ForceSurface& tracSurface,
                 const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition,
                 double horizon, PeriPMBMaterialP pmb_material)
   : Box(minPoint, maxPoint, xNumPoints, yNumPoints, zNumPoints,
         density, poission, young,
         bodyForce, surForce, tracSurface,
         boundary, boundaryPosition), 
     d_horizon(horizon), d_pmb_material(pmb_material),
     d_areaX(0.0), d_areaY(0.0), d_areaZ(0.0),
     d_boundary_velocity(0.0),
     d_velocity_surface(none), d_peri_type(Bond) 
{
//  std::cout << xNumPoints << " " << yNumPoints << " " << zNumPoints << std::endl;
  setPeriPoints();
//  std::cout << "Hi " << std::endl;
//  std::cout << "Horizon= " << getHorizon() << std::endl; 
  createPeriBondsOfEachPeriPoints(getHorizon());
  calculateMicromodulusFunction();
  findPointsAreaBond();
  calculateBondsVolumeBond();
//  calculateBondsSurfaceCorrectionFactor();
//  findPointsMassAndVolume();
//  findPointsArea();
//  std::cout << "Areas are successfully calculated. " << std::endl; 
  calculatePointsBodyForceBond();

  findZeroEnergyConstantTwo();

}



PeriBox::PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint,
                const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
                const double& density, const double& poission, const double& young,
                const double& bodyForce,
                const double& surForce, const ForceSurface& tracSurface,
                const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition, double horizon,
                PeriPMBMaterialP pmb_material, const ForceSurface& VelocitySurface, const Vector3D& BoundaryVelocity)
   : Box(minPoint, maxPoint, xNumPoints, yNumPoints, zNumPoints,
         density, poission, young,
         bodyForce, surForce, tracSurface,
         boundary, boundaryPosition), 
     d_horizon(horizon), d_pmb_material(pmb_material),
     d_areaX(0.0), d_areaY(0.0), d_areaZ(0.0),
     d_boundary_velocity(BoundaryVelocity),
     d_velocity_surface(VelocitySurface), d_peri_type(Bond) 
{
//  std::cout << xNumPoints << " " << yNumPoints << " " << zNumPoints << std::endl;
  setPeriPoints();
//  std::cout << "Hi " << std::endl;
//  std::cout << "Horizon= " << getHorizon() << std::endl; 
  createPeriBondsOfEachPeriPoints(getHorizon());
  calculateMicromodulusFunction();
  findPointsAreaBond();
  calculateBondsVolumeBond();
//  calculateBondsSurfaceCorrectionFactor();
//  findPointsMassAndVolume();
//  findPointsArea();
//  std::cout << "Areas are successfully calculated. " << std::endl; 
  calculatePointsBodyForceBond();

  findZeroEnergyConstantTwo();

}


PeriBox::PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint,
                const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
                const double& density, const double& poission, const double& young,
                const double& bodyForce,
                const double& surForce, const ForceSurface& tracSurface,
                const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition, double horizon,
                PeriPMBMaterialP pmb_material, const PeriType& type, 
                const ForceSurface& VelocitySurface, const Vector3D& BoundaryVelocity)
   : Box(minPoint, maxPoint, xNumPoints, yNumPoints, zNumPoints,
         density, poission, young,
         bodyForce, surForce, tracSurface,
         boundary, boundaryPosition), 
     d_horizon(horizon), d_pmb_material(pmb_material),
     d_areaX(0.0), d_areaY(0.0), d_areaZ(0.0),
     d_boundary_velocity(BoundaryVelocity),
     d_velocity_surface(VelocitySurface), d_peri_type(type) 
{
//  std::cout << xNumPoints << " " << yNumPoints << " " << zNumPoints << std::endl;
  if (type == Bond)
  {
    setPeriPoints();
//  std::cout << "Hi " << std::endl;
//  std::cout << "Horizon= " << getHorizon() << std::endl; 
    createPeriBondsOfEachPeriPoints(getHorizon());
    calculateMicromodulusFunction();
    findPointsAreaBond();
    calculateBondsVolumeBond();
    calculatePointsBodyForceBond();
  }
  else
  {
    setPeriPoints();
    setStatePoints();
    createStateBondsOfEachStatePoints(getHorizon());
    findPointsAreaState(); 
    calculateBondsVolumeState();
    calculatePointsBodyForceState();

    calcualtePmbConstantSpring();

    findZeroEnergyConstantTwo();

  }   
//  calculateMicromodulusFunction();
//  findPointsArea();
//  calculateBondsVolume();
//  calculateBondsSurfaceCorrectionFactor();
//  findPointsMassAndVolume();
//  findPointsArea();
//  std::cout << "Areas are successfully calculated. " << std::endl; 
//  calculatePointsBodyForce();

}


PeriBox::PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint,
                const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
                const double& density, const double& poission, const double& young,
                const double& yield, const double& bodyForce,
                const double& surForce, const ForceSurface& tracSurface,
                const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition, 
                double horizon, const PeriType& type, 
                const ForceSurface& VelocitySurface, const Vector3D& BoundaryVelocity)
   : Box(minPoint, maxPoint, xNumPoints, yNumPoints, zNumPoints,
         density, poission, young, yield,
         bodyForce, surForce, tracSurface,
         boundary, boundaryPosition), 
     d_horizon(horizon), d_pmb_material(new PeriPMBMaterial(density, horizon, young, false)),
     d_areaX(0.0), d_areaY(0.0), d_areaZ(0.0),
     d_boundary_velocity(BoundaryVelocity),
     d_velocity_surface(VelocitySurface), d_peri_type(type) 
{
//  std::cout << xNumPoints << " " << yNumPoints << " " << zNumPoints << std::endl;
  if (type == Bond)
  {
//    d_pmb_material(new PeriPMBMaterial());
    getPMB()->setFractureEnergy(15000);
    getPMB()->calculateCriticalStretch();
    getPMB()->calculateMicromodulusFunction();
    setPeriPoints();

    
//  std::cout << "Hi " << std::endl;
//  std::cout << "Horizon= " << getHorizon() << std::endl; 
//    createPeriBondsOfEachPeriPoints(getHorizon());
    calculateMicromodulusFunction();
    findPointsAreaBond();
    calculateBondsVolumeBond();
    calculatePointsBodyForceBond();
  }
  else
  {
//    d_pmb_material(new PeriPMBMaterial(density, horizon, young, false));
    setPeriPoints();
    setStatePoints();
//    std::cout << "Hi " << std::endl;
    createStateBondsOfEachStatePoints(getHorizon());
    findPointsAreaState(); 
    calculateBondsVolumeState();
    calculatePointsBodyForceState();

    calcualtePmbConstantSpring();

    findZeroEnergyConstantTwo();

  }   
//  calculateMicromodulusFunction();
//  findPointsArea();
//  calculateBondsVolume();
//  calculateBondsSurfaceCorrectionFactor();
//  findPointsMassAndVolume();
//  findPointsArea();
//  std::cout << "Areas are successfully calculated. " << std::endl; 
//  calculatePointsBodyForce();

}



//*********************************************** State-based Peridynamics ********************************

PeriBox::PeriBox(const Vector3D& minPoint, const Vector3D& maxPoint,
                const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
                const double& density, const double& poission, const double& young,
                const double& yield, const double& bodyForce,
                const double& surForce, const ForceSurface& tracSurface,
                const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition, 
                double horizon, const PeriType& type, 
                const ForceSurface& VelocitySurface, 
                const Vector3D& BoundaryVelocity, const Matrix3D& InitialVelocityGradient)
   : Box(minPoint, maxPoint, xNumPoints, yNumPoints, zNumPoints,
         density, poission, young, yield,
         bodyForce, surForce, tracSurface,
         boundary, boundaryPosition), 
     d_horizon(horizon), d_pmb_material(new PeriPMBMaterial(density, horizon, young, false)),
     d_areaX(0.0), d_areaY(0.0), d_areaZ(0.0),
     d_boundary_velocity(BoundaryVelocity),
     d_velocity_surface(VelocitySurface), d_peri_type(type) 
{
//  std::cout << xNumPoints << " " << yNumPoints << " " << zNumPoints << std::endl;
  if (type == Bond)
  {
//    d_pmb_material(new PeriPMBMaterial());
    getPMB()->setFractureEnergy(15000);
    getPMB()->calculateCriticalStretch();
    getPMB()->calculateMicromodulusFunction();
    setPeriPoints();
//  std::cout << "Hi " << std::endl;
//  std::cout << "Horizon= " << getHorizon() << std::endl; 
    createPeriBondsOfEachPeriPoints(getHorizon());
    calculateMicromodulusFunction();
    findPointsAreaBond();
    calculateBondsVolumeBond();
    calculatePointsBodyForceBond();
  }
  else
  {
//    d_pmb_material(new PeriPMBMaterial(density, horizon, young, false));
    setPeriPoints();
    setStatePoints(InitialVelocityGradient);
//    std::cout << "Hi " << std::endl;
//    std::cout << "Horizon= " << getHorizon() << std::endl;


//    int nx = getNumPointsX();
//    double dx = (getMaxPoint().x()-getMinPoint().x())/(nx-1);
//    createStateBondsOfEachStatePoints(getHorizon(), dx);

    createStateBondsOfEachStatePoints(getHorizon());
    findPointsAreaState(); 
    calculateBondsVolumeState();
    calculatePointsBodyForceState();

    calcualtePmbConstantSpring();

    findZeroEnergyConstantTwo();

  }   

}

//************************************************************************************************************


PeriBox::~PeriBox()
{


}



void 
PeriBox::findPointsMassAndVolume()
{


  double sum = 0.0;
  int condition = 8;
  double density = Body::getDensity();
  double volume = 0;
  int size = 0;
  volume = (getMaxPoint().x() - getMinPoint().x())*(getMaxPoint().y() - getMinPoint().y())
                                                  *(getMaxPoint().z() - getMinPoint().z());
  double unit_mass = volume*density/(8*(getNumPointsX()-1)*(getNumPointsY()-1)*(getNumPointsZ()-1));
//  std::vector <PeriMaterialPointP> pointsArray = getMaterialPeriPoints();
  if (getPeriType() == Bond)
  {
    std::vector <PeriMaterialPointP> pointsArray = getMaterialPeriPoints();
    size = pointsArray.size();
    for (int i = 0; i < size; i++)
    {
      PeriMaterialPointP cur_point = pointsArray[i]; 
  //    StateMaterialPointP cur_point = pointsArray[i]; 
      if ((cur_point->getPosOld().x() == getMaxPoint().x()) || (cur_point->getPosOld().x() == getMinPoint().x()))
          condition /= 2;        
      if ((cur_point->getPosOld().y() == getMaxPoint().y()) || (cur_point->getPosOld().y() == getMinPoint().y()))
          condition /= 2;        
      if ((cur_point->getPosOld().z() == getMaxPoint().z()) || (cur_point->getPosOld().z() == getMinPoint().z()))
          condition /= 2;
      cur_point->setMass(condition*unit_mass);
      cur_point->setInitialVolume(condition*unit_mass/density);

      sum += condition*unit_mass;

//    std::cout << cur_point->getID() << "   " << "Condition: " << condition << "  Mass: " << cur_point->getMass() 
//              << std::endl;

//    std::cout << cur_point->getID() << "   " << "Condition: " << condition << "  Volume: " << cur_point->getVolume() 
//              << std::endl;


      condition = 8;
              
     }

  }
  else
  {
    std::vector <StateMaterialPointP> pointsStateArray;
    size = pointsStateArray.size();
    for (int i = 0; i < size; i++)
    {
//      PeriMaterialPointP cur_point = pointsArray[i]; 
      StateMaterialPointP cur_point = pointsStateArray[i]; 
      if ((cur_point->getPosOld().x() == getMaxPoint().x()) || (cur_point->getPosOld().x() == getMinPoint().x()))
          condition /= 2;        
      if ((cur_point->getPosOld().y() == getMaxPoint().y()) || (cur_point->getPosOld().y() == getMinPoint().y()))
          condition /= 2;        
      if ((cur_point->getPosOld().z() == getMaxPoint().z()) || (cur_point->getPosOld().z() == getMinPoint().z()))
          condition /= 2;
      cur_point->setMass(condition*unit_mass);
      cur_point->setInitialVolume(condition*unit_mass/density);

      sum += condition*unit_mass;

//    std::cout << cur_point->getID() << "   " << "Condition: " << condition << "  Mass: " << cur_point->getMass() 
//              << std::endl;

//    std::cout << cur_point->getID() << "   " << "Condition: " << condition << "  Volume: " << cur_point->getVolume() 
//              << std::endl;


      condition = 8;
              
     }

   }    
//  size = pointsArray.size(); 
//  for (int i = 0; i < size; i++)
//  {
//    PeriMaterialPointP cur_point = pointsArray[i]; 
//    StateMaterialPointP cur_point = pointsArray[i]; 
//    if ((cur_point->getPosOld().x() == getMaxPoint().x()) || (cur_point->getPosOld().x() == getMinPoint().x()))
//        condition /= 2;        
//    if ((cur_point->getPosOld().y() == getMaxPoint().y()) || (cur_point->getPosOld().y() == getMinPoint().y()))
//        condition /= 2;        
//    if ((cur_point->getPosOld().z() == getMaxPoint().z()) || (cur_point->getPosOld().z() == getMinPoint().z()))
//        condition /= 2;
//    cur_point->setMass(condition*unit_mass);
//    cur_point->setInitialVolume(condition*unit_mass/density);

//    sum += condition*unit_mass;

//    std::cout << cur_point->getID() << "   " << "Condition: " << condition << "  Mass: " << cur_point->getMass() 
//              << std::endl;

//    std::cout << cur_point->getID() << "   " << "Condition: " << condition << "  Volume: " << cur_point->getVolume() 
//              << std::endl;


//    condition = 8;
              
//  }

//  std::cout << "Total mass is: " << sum << std::endl;

}



void
PeriBox::findPointsAreaBond()
{
//  std::cout << "----------------findPointsArea PeriBox.cc--------------------" << std::endl;
  Vector3D force_vec(0.0, 0.0, 0.0);
//  std::vector <MaterialPointP> pointsArray = getMaterialPoints();
  std::vector <PeriMaterialPointP> pointsArray = getMaterialPeriPoints();

  int nx = getNumPointsX();
  int ny = getNumPointsY();
  int nz = getNumPointsZ();

  double dx = (getMaxPoint().x()-getMinPoint().x())/(nx-1);
  double dy = (getMaxPoint().y()-getMinPoint().y())/(ny-1);
  double dz = (getMaxPoint().z()-getMinPoint().z())/(nz-1);

  std::cout << "dx= " << dx << " dy= " << dy << " dz= " << dz << std::endl; // www

  ForceSurface  trac_surf = getTractionSurface();  
  double force = getTractionForce();

//  std::cout << "traction force= " << force << std::endl;

  int size = pointsArray.size();

//  std::cout << "Number of points= " << size << std::endl;

  for (int i = 0; i < size; i++)
  {
//    MaterialPointP cur_point = pointsArray[i];
    PeriMaterialPointP cur_point = pointsArray[i];
    bool surfaceOnX = ((cur_point->x() == getMinPoint().x()) || (cur_point->x() == getMaxPoint().x()));
    bool surfaceOnY = ((cur_point->y() == getMinPoint().y()) || (cur_point->y() == getMaxPoint().y()));
    bool surfaceOnZ = ((cur_point->z() == getMinPoint().z()) || (cur_point->z() == getMaxPoint().z()));

    if (surfaceOnX)
    {
      setAreaX(dy*dz);
      if (surfaceOnY)
      {
        setAreaX(getAreaX()/2);
        setAreaY(0.5*dx*dz);
        if (surfaceOnZ)
        {
          setAreaX(getAreaX()/2);
          setAreaY(getAreaY()/2);
          setAreaZ(0.25*dy*dx);
        }
      }
      else if (surfaceOnZ)
      {
        setAreaX(getAreaX()/2);
        setAreaZ(0.5*dx*dy);
      }        
    }
    else if (surfaceOnY)
    {
      setAreaY(dx*dz);
      if (surfaceOnZ)
      {
        setAreaY(getAreaY()/2);;
        setAreaZ(0.5*dx*dy);
      }
    }
    else if (surfaceOnZ)
    {
      setAreaZ(dx*dy);
    }

    

    
    switch (trac_surf)
    {
      case right:
      if ((surfaceOnX) && (cur_point->x() == getMaxPoint().x()))
      {
//      force_vec.x(-1*force);
        force_vec.x(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaX());
      }
      break;

      case left:
      if ((surfaceOnX) && (cur_point->x() == getMinPoint().x()))
      {
//      force_vec.x(force);
        force_vec.x(-1*force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaX());
      }
      break;

      case up:
      if ((surfaceOnY) && (cur_point->y() == getMaxPoint().y()))
      {
//      force_vec.y(-1*force);
        force_vec.y(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaY());
      }            
      break;

      case down:
      if ((surfaceOnY) && (cur_point->y() == getMinPoint().y()))
      {
        force_vec.y(-1*force);
//        force_vec.y(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaY());
      }            
      break;

      case front:
      if ((surfaceOnZ) && (cur_point->z() == getMaxPoint().z()))
      {
//      force_vec.z(-1*force);
        force_vec.z(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaZ());
      }
      break;

      case back:
      if ((surfaceOnZ) && (cur_point->z() == getMaxPoint().z()))
      {
        force_vec.z(-1*force);
//        force_vec.z(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaZ());
      }
      break;

      case upDownX:
      if ((surfaceOnX) && (cur_point->x() == getMinPoint().x()))
      {
//      force_vec.y(-1*force);
        force_vec.x(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaX());
      }            
      else if ((surfaceOnX) && (cur_point->x() == getMaxPoint().x()))
      {
        force_vec.x(-1*force);
//        force_vec.x(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaX());
      }
      break;


      case upDownY:
      if ((surfaceOnY) && (cur_point->y() == getMinPoint().y()))
      {
//      force_vec.y(-1*force);
        force_vec.y(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaY());
      }            
      else if ((surfaceOnY) && (cur_point->y() == getMaxPoint().y()))
      {
        force_vec.y(-1*force);
//        force_vec.y(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaY());
      }
      break;

      case upDownZ:
      if ((surfaceOnZ) && (cur_point->z() == getMinPoint().z()))
      {
//      force_vec.z(-1*force);
        force_vec.z(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaZ());
      }            
      else if ((surfaceOnZ) && (cur_point->z() == getMaxPoint().z()))
      {
        force_vec.z(-1*force);
//        force_vec.z(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaZ());
      }
      break;


      case none:
        cur_point->surfaceForce(force_vec);
      break;

      default:
        std::cout <<
        " The surface can be right, left, up, down, front, back, or none which means we do not have any traction force."
        << std::endl;
      break;

    } // end of switch

//    if (cur_point->getInitialArea() != 0.0)
//    std::cout << "Point " << cur_point->getID() << " " << cur_point->getPosOld() << // " surfaceOnX= " << surfaceOnX <<
//                                                   " surfaceOnY= " << surfaceOnY <<
//                                                   " surfaceOnZ= " << surfaceOnZ <<
//                                                   " surfaceForce= " << cur_point->surfaceForce();// << std::endl;
//                                                   " area= "  << cur_point->getInitialArea() << std::endl;

  } // end of for


//  setMaterialPoints(pointsArray);
  setMaterialPeriPoints(pointsArray);

  
} // end of function



void
PeriBox::findPointsAreaState()
{
//  std::cout << "----------------findPointsArea PeriBox.cc--------------------" << std::endl;
  Vector3D force_vec(0.0, 0.0, 0.0);
//  std::vector <MaterialPointP> pointsArray = getMaterialPoints();
  std::vector <StateMaterialPointP> pointsArray = getStatePoints();

  int nx = getNumPointsX();
  int ny = getNumPointsY();
  int nz = getNumPointsZ();

  double dx = (getMaxPoint().x()-getMinPoint().x())/(nx-1);
  double dy = (getMaxPoint().y()-getMinPoint().y())/(ny-1);
  double dz = (getMaxPoint().z()-getMinPoint().z())/(nz-1);

//  std::cout << "dx= " << dx << " dy= " << dy << " dz= " << dz << std::endl; // www

  ForceSurface  trac_surf = getTractionSurface();  
  double force = getTractionForce();

//  std::cout << "traction force= " << force << std::endl;

  int size = pointsArray.size();

//  std::cout << "Number of points= " << size << std::endl;

  for (int i = 0; i < size; i++)
  {
//    MaterialPointP cur_point = pointsArray[i];
    StateMaterialPointP cur_point = pointsArray[i];
    bool surfaceOnX = ((cur_point->x() == getMinPoint().x()) || (cur_point->x() == getMaxPoint().x()));
    bool surfaceOnY = ((cur_point->y() == getMinPoint().y()) || (cur_point->y() == getMaxPoint().y()));
    bool surfaceOnZ = ((cur_point->z() == getMinPoint().z()) || (cur_point->z() == getMaxPoint().z()));

    if (surfaceOnX)
    {
      setAreaX(dy*dz);
      if (surfaceOnY)
      {
        setAreaX(getAreaX()/2);
        setAreaY(0.5*dx*dz);
        if (surfaceOnZ)
        {
          setAreaX(getAreaX()/2);
          setAreaY(getAreaY()/2);
          setAreaZ(0.25*dy*dx);
        }
      }
      else if (surfaceOnZ)
      {
        setAreaX(getAreaX()/2);
        setAreaZ(0.5*dx*dy);
      }        
    }
    else if (surfaceOnY)
    {
      setAreaY(dx*dz);
      if (surfaceOnZ)
      {
        setAreaY(getAreaY()/2);;
        setAreaZ(0.5*dx*dy);
      }
    }
    else if (surfaceOnZ)
    {
      setAreaZ(dx*dy);
    }

    

    
    switch (trac_surf)
    {
      case right:
      if ((surfaceOnX) && (cur_point->x() == getMaxPoint().x()))
      {
//      force_vec.x(-1*force);
        force_vec.x(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaX());
      }
      break;

      case left:
      if ((surfaceOnX) && (cur_point->x() == getMinPoint().x()))
      {
//      force_vec.x(force);
        force_vec.x(-1*force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaX());
      }
      break;

      case up:
      if ((surfaceOnY) && (cur_point->y() == getMaxPoint().y()))
      {
//      force_vec.y(-1*force);
        force_vec.y(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaY());
      }            
      break;

      case down:
      if ((surfaceOnY) && (cur_point->y() == getMinPoint().y()))
      {
        force_vec.y(-1*force);
//        force_vec.y(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaY());
      }            
      break;

      case front:
      if ((surfaceOnZ) && (cur_point->z() == getMaxPoint().z()))
      {
//      force_vec.z(-1*force);
        force_vec.z(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaZ());
      }
      break;

      case back:
      if ((surfaceOnZ) && (cur_point->z() == getMaxPoint().z()))
      {
        force_vec.z(-1*force);
//        force_vec.z(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaZ());
      }
      break;

      case upDownX:
      if ((surfaceOnX) && (cur_point->x() == getMinPoint().x()))
      {
//      force_vec.y(-1*force);
        force_vec.x(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaX());
      }            
      else if ((surfaceOnX) && (cur_point->x() == getMaxPoint().x()))
      {
        force_vec.x(-1*force);
//        force_vec.x(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaX());
      }
      break;


      case upDownY:
      if ((surfaceOnY) && (cur_point->y() == getMinPoint().y()))
      {
//      force_vec.y(-1*force);
        force_vec.y(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaY());
      }            
      else if ((surfaceOnY) && (cur_point->y() == getMaxPoint().y()))
      {
        force_vec.y(-1*force);
//        force_vec.y(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaY());
      }
      break;

      case upDownZ:
      if ((surfaceOnZ) && (cur_point->z() == getMinPoint().z()))
      {
//      force_vec.z(-1*force);
        force_vec.z(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaZ());
      }            
      else if ((surfaceOnZ) && (cur_point->z() == getMaxPoint().z()))
      {
        force_vec.z(-1*force);
//        force_vec.z(force);
        cur_point->surfaceForce(force_vec);
        cur_point->setInitialArea(getAreaZ());
      }
      break;


      case none:
        cur_point->surfaceForce(force_vec);
      break;

      default:
        std::cout <<
        " The surface can be right, left, up, down, front, back, or none which means we do not have any traction force."
        << std::endl;
      break;

    } // end of switch

//    if (cur_point->getInitialArea() != 0.0)
//    std::cout << "Point " << cur_point->getID() << " " << cur_point->getPosOld() << // " surfaceOnX= " << surfaceOnX <<
//                                                   " surfaceOnY= " << surfaceOnY <<
//                                                   " surfaceOnZ= " << surfaceOnZ <<
//                                                   " surfaceForce= " << cur_point->surfaceForce();// << std::endl;
//                                                   " area= "  << cur_point->getInitialArea() << std::endl;

  } // end of for


//  setMaterialPoints(pointsArray);
  setStatePoints(pointsArray);

  
} // end of function



void
PeriBox::calculatePointsBodyForceBond()
{
  Vector3D zero(0.0, 0.0, 0.0);
  std::vector <PeriMaterialPointP> pointsArray = getMaterialPeriPoints();
  int size = pointsArray.size();
  for (int i = 0; i < size; i++)
  {
    PeriMaterialPointP cur_point = pointsArray[i];
    cur_point->setBodyForce(cur_point->getBodyForce() + 
                            cur_point->surfaceForce()*cur_point->getInitialArea()/cur_point->getInitialVolume());

  }

//  int num = 0;
//  for (int j = 0; j < size; j++)
//  {
//    if ((getMaterialPeriPoints()[j]->getBodyForce().x() != 0.0) ||
//        (getMaterialPeriPoints()[j]->getBodyForce().y() != 0.0) ||
//        (getMaterialPeriPoints()[j]->getBodyForce().z() != 0.0))
//        {
//          num++;
//          std::cout << "Num= " << num;
//          std::cout << " ID= " << getMaterialPeriPoints()[j]->getID() <<
//                       " Pos= " << getMaterialPeriPoints()[j]->getPosOld() <<
//                       " Force= " << getMaterialPeriPoints()[j]->getBodyForce() << std::endl;
//        }
//  }
  


}


void
PeriBox::calculatePointsBodyForceState()
{
  Vector3D zero(0.0, 0.0, 0.0);
  std::vector <StateMaterialPointP> pointsArray = getStatePoints();
  int size = pointsArray.size();
  for (int i = 0; i < size; i++)
  {
    StateMaterialPointP cur_point = pointsArray[i];
    cur_point->setBodyForce(cur_point->getBodyForce() + 
                            cur_point->surfaceForce()*cur_point->getInitialArea()/cur_point->getInitialVolume());

//    if (cur_point->getID() == 1055)
//       std::cout << "BodyForce= " << cur_point->getBodyForce() << std::endl;

  }

//  int num = 0;
//  for (int j = 0; j < size; j++)
//  {
//    if ((getMaterialPeriPoints()[j]->getBodyForce().x() != 0.0) ||
//        (getMaterialPeriPoints()[j]->getBodyForce().y() != 0.0) ||
//        (getMaterialPeriPoints()[j]->getBodyForce().z() != 0.0))
//        {
//          num++;
//          std::cout << "Num= " << num;
//          std::cout << " ID= " << getMaterialPeriPoints()[j]->getID() <<
//                       " Pos= " << getMaterialPeriPoints()[j]->getPosOld() <<
//                       " Force= " << getMaterialPeriPoints()[j]->getBodyForce() << std::endl;
//        }
//  }
  


}


//PeriBox::createPeriMaterialPointArray()
//{
  


//}


void
PeriBox::calculateBondsSurfaceCorrectionFactor()
{
  std::vector <PeriMaterialPointP> pointsArray = getMaterialPeriPoints();
  int size = pointsArray.size();
  for (int i = 0; i < size; i++)
  {
    PeriMaterialPointP cur_point = pointsArray[i];  
    PeriBondPArray bondsArray = cur_point->getHorizonBonds();
    for (int j = 0; j < bondsArray.size(); j++)
    {
      PeriBondP cur_bond = bondsArray[j];
      getPMB()->surfaceCorrectionFactor(cur_bond);
//      std::cout << "centerID= " << cur_bond->getCenterPoint()->getID() 
//                << " secondID= " << cur_bond->getSecondPoint()->getID()
//                << " surfaceCorrection= " << cur_bond->getSurfaceCorrectionFactor() << std::endl;
    }
  }

}
  


void
PeriBox::calculateMicromodulusFunction()
{
  double delX = (getMaxPoint().x()-getMinPoint().x())/(getNumPointsX()-1);
  int k = getHorizon()/delX;
  PeriPMBMaterialP pmbMaterial = getPMB();
  double microModulus = (pmbMaterial->getYoungsModulus())/(pow(delX, 4));
  if (k == 2)
  {
    microModulus *= 0.302942;
  }
  else if (k == 3)
  {
    microModulus *= 0.052385;
  }    
  else if (k == 4)
  {
    microModulus *= 0.017290;
  }    
  else 
  {
    microModulus *= 0.006819;
  }    
    
  pmbMaterial->setMicromodulusFunction(microModulus);

  std::cout << "delX= " << delX << " horizon/delX= " << getHorizon()/delX << std::endl;
  std::cout << "Young= " << pmbMaterial->getYoungsModulus() << " MicroFunction= " 
                         << pmbMaterial->getMicromodulusFunction() << std::endl;

}


void
PeriBox::calculateBondsVolumeBond()
{
//  std::cout << "----------------calculateBondsVolume PeriBox.cc--------------------" << std::endl;
  double dx = (getMaxPoint().x()-getMinPoint().x())/(getNumPointsX()-1);
  double horizon = getHorizon();
  std::vector <PeriMaterialPointP> pointsArray = getMaterialPeriPoints();
  int size = pointsArray.size();
  for (int i = 0; i < size; i++)
  {
    PeriMaterialPointP cur_point = pointsArray[i];
    cur_point->setInitialArea(dx*dx);
    cur_point->setInitialVolume(dx*dx*dx);  
    PeriBondPArray bondsArray = cur_point->getHorizonBonds();
    for (int j = 0; j < bondsArray.size(); j++)
    {
      PeriBondP cur_bond = bondsArray[j];
      cur_bond->calculateVolume(horizon, dx/2, dx);
//      getPMB()->surfaceCorrectionFactor(cur_bond);
//      std::cout << "centerID= " << cur_bond->getCenterPoint()->getID() 
//                << " secondID= " << cur_bond->getSecondPoint()->getID()
//                << " volume= " << cur_bond->getVolume() << std::endl;
    }
  }

}


void
PeriBox::calculateBondsVolumeState()
{
//  std::cout << "----------------calculateBondsVolume PeriBox.cc--------------------" << std::endl;
  double dx = (getMaxPoint().x() - getMinPoint().x())/(getNumPointsX() - 1);
  double horizon = getHorizon();
  std::vector <StateMaterialPointP> pointsArray = getStatePoints();
  int size = pointsArray.size();
  for (int i = 0; i < size; i++)
  {
    StateMaterialPointP cur_point = pointsArray[i];
    cur_point->setInitialArea(dx*dx);
    cur_point->setInitialVolume(dx*dx*dx);
  }


  for (int i = 0; i < size; i++)
  {
    StateMaterialPointP cur_point = pointsArray[i];
//    std::cout << "ID= " << cur_point->getID() << " Volume= " << cur_point->getInitialVolume()
//                                              << " Area= " << cur_point->getInitialArea() << std::endl;  
    StateBondPArray  bondsArray = cur_point->getStateBonds();
    for (int j = 0; j < bondsArray.size(); j++)
    {
      StateBondP cur_bond = bondsArray[j];
      cur_bond->calculateVolume(horizon, dx/2, dx);
//      getPMB()->surfaceCorrectionFactor(cur_bond);
//      std::cout << "centerID= " << cur_bond->getCenterPoint()->getID() 
//                << " secondID= " << cur_bond->getSecondPoint()->getID()
//                << " volume= " << cur_bond->getVolume() << std::endl;
    }
  }

}



void
PeriBox::applyBoundaryVelocity()
{

  Vector3D velocity_vec(0.0, 0.0, 0.0);

  std::vector <PeriMaterialPointP> pointsArray = getMaterialPeriPoints();

  ForceSurface  velocity_surf = getVelocityBoundarySurface();  
  Vector3D velocity = getBoundaryVelocity();

  int size = pointsArray.size();

  for (int i = 0; i < size; i++)
  {
    PeriMaterialPointP cur_point = pointsArray[i];
        
    switch (velocity_surf)
    {
      case right:
      if (cur_point->x() == getMaxPoint().x())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }
      break;

      case left:
      if (cur_point->x() == getMinPoint().x())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }
      break;

      case up:
      if (cur_point->y() == getMaxPoint().y())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }            
      break;

      case down:
      if (cur_point->y() == getMinPoint().y())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }            
      break;

      case front:
      if (cur_point->z() == getMaxPoint().z())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }
      break;

      case back:
      if (cur_point->z() == getMinPoint().z())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }
      break;

      case upDownX:
      if (cur_point->x() == getMinPoint().x())
      {
        cur_point->setVelOld(getBoundaryVelocity()*(-1));
        cur_point->setVelMid(getBoundaryVelocity()*(-1));
        cur_point->setVelNew(getBoundaryVelocity()*(-1));
        cur_point->setVelocityBoundaryFlag(true);
      }            
      else if (cur_point->x() == getMaxPoint().x())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }
      break;


      case upDownY:
      if (cur_point->y() == getMinPoint().y())
      {
        cur_point->setVelOld(getBoundaryVelocity()*(-1));
        cur_point->setVelMid(getBoundaryVelocity()*(-1));
        cur_point->setVelNew(getBoundaryVelocity()*(-1));
        cur_point->setVelocityBoundaryFlag(true);
      }            
      else if (cur_point->y() == getMaxPoint().y())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }
      break;

      case upDownZ:
      if (cur_point->z() == getMinPoint().z())
      {
        cur_point->setVelOld(getBoundaryVelocity()*(-1));
        cur_point->setVelMid(getBoundaryVelocity()*(-1));
        cur_point->setVelNew(getBoundaryVelocity()*(-1));
        cur_point->setVelocityBoundaryFlag(true);
      }            
      else if (cur_point->z() == getMaxPoint().z())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }
      break;


      case none:
//        cur_point->surfaceForce(force_vec);
      break;

      default:
        std::cout <<
        " The surface can be right, left, up, down, front, back, or none which means we do not have any traction force."
        << std::endl;
      break;

    } // end of switch

//    if (cur_point->getInitialArea() != 0.0)
//    std::cout << "Point " << cur_point->getID() << " " << cur_point->getPosOld() << // " surfaceOnX= " << surfaceOnX <<
//                                                   
//                                                   " velocity " << cur_point->getVelOld() << std::endl;
//                                                   //" newVel "  << cur_point->getVelNew() << std::endl;

  } // end of for

  setMaterialPeriPoints(pointsArray);
   
}




void
PeriBox::applyBoundaryVelocityState()
{

  Vector3D velocity_vec(0.0, 0.0, 0.0);

  std::vector <StateMaterialPointP> pointsArray = getStatePoints();

//  std::vector <StateMaterialPointP> pointsPositive;
//  std::vector <StateMaterialPointP> pointsNegative;

//  pointsPositive.clear();
//  pointsNegative.clear();



  ForceSurface  velocity_surf = getVelocityBoundarySurface();  
  Vector3D velocity = getBoundaryVelocity();

  int size = pointsArray.size();

  for (int i = 0; i < size; i++)
  {
    StateMaterialPointP cur_point = pointsArray[i];
        
    switch (velocity_surf)
    {
      case right:
      if (cur_point->x() == getMaxPoint().x())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
//        pointsPositive.emplace_back(cur_point);
        

      }
      break;

      case left:
      if (cur_point->x() == getMinPoint().x())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }
      break;

      case up:
      if (cur_point->y() == getMaxPoint().y())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }            
      break;

      case down:
      if (cur_point->y() == getMinPoint().y())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }            
      break;

      case front:
      if (cur_point->z() == getMaxPoint().z())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }
      break;

      case back:
      if (cur_point->z() == getMinPoint().z())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }
      break;

      case upDownX:
      if (cur_point->x() == getMinPoint().x())
      {
        cur_point->setVelOld(getBoundaryVelocity()*(-1));
        cur_point->setVelMid(getBoundaryVelocity()*(-1));
        cur_point->setVelNew(getBoundaryVelocity()*(-1));
        cur_point->setVelocityBoundaryFlag(true);
      }            
      else if (cur_point->x() == getMaxPoint().x())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }
      break;


      case upDownY:
      if (cur_point->y() == getMinPoint().y())
      {
        cur_point->setVelOld(getBoundaryVelocity()*(-1));
        cur_point->setVelMid(getBoundaryVelocity()*(-1));
        cur_point->setVelNew(getBoundaryVelocity()*(-1));
        cur_point->setVelocityBoundaryFlag(true);
      }            
      else if (cur_point->y() == getMaxPoint().y())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }
      break;

      case upDownZ:
      if (cur_point->z() == getMinPoint().z())
      {
        cur_point->setVelOld(getBoundaryVelocity()*(-1));
        cur_point->setVelMid(getBoundaryVelocity()*(-1));
        cur_point->setVelNew(getBoundaryVelocity()*(-1));
        cur_point->setVelocityBoundaryFlag(true);
      }            
      else if (cur_point->z() == getMaxPoint().z())
      {
        cur_point->setVelOld(getBoundaryVelocity());
        cur_point->setVelMid(getBoundaryVelocity());
        cur_point->setVelNew(getBoundaryVelocity());
        cur_point->setVelocityBoundaryFlag(true);
      }
      break;

      case tear:
      {
        int nnn = 3;
        double xxx = cur_point->x();
        double zzz = cur_point->z();

        double xLength = getMaxPoint().x() - getMinPoint().x();
        double dx = xLength/(getNumPointsX() - 1);
        double midXposition = getMinPoint().x() + xLength/2;
        double min = midXposition - nnn*dx;
        double max = midXposition + nnn*dx;
        bool midRight = (xxx <= max) && (xxx >= midXposition);
        bool midLeft =  (xxx < midXposition) && (xxx >= min);

        if ((zzz == getMaxPoint().z()) && midRight)
        {
          cur_point->setVelOld(getBoundaryVelocity()*(-1));
          cur_point->setVelMid(getBoundaryVelocity()*(-1));
          cur_point->setVelNew(getBoundaryVelocity()*(-1));
          cur_point->setVelocityBoundaryFlag(true);
          
        }            
        else if ((zzz == getMaxPoint().z()) && midLeft)
        {
          cur_point->setVelOld(getBoundaryVelocity());
          cur_point->setVelMid(getBoundaryVelocity());
          cur_point->setVelNew(getBoundaryVelocity());
          cur_point->setVelocityBoundaryFlag(true);
        }
      }
      break; 




      case none:
//        cur_point->surfaceForce(force_vec);
      break;

      default:
        std::cout <<
        " The surface can be right, left, up, down, front, back, or none which means we do not have any traction force."
        << std::endl;
      break;

    } // end of switch

//    if (cur_point->getInitialArea() != 0.0)
//    if (cur_point->getVelocityBoundaryFlag() == true)
//    {
//       std::cout << "MinPointBox= " << getMinPoint() << " MinPointBox= " 
//                                    << getMaxPoint() << std::endl;

//       std::cout << "Point " << cur_point->getID() << " "
//                              << cur_point->getPosOld() << " NewPos= " << cur_point->getPosNew() << std::endl
                                                         // << " surfaceOnX= " << surfaceOnX
//                                                   
//                                                   << "Velocity " << cur_point->getVelOld() <<
//                                                  " NewVel "  << cur_point->getVelNew() << std::endl << std::endl;
//    }

  } // end of for

  setStatePoints(pointsArray);
   
}



void
PeriBox::showBoundaryVelocityPoints()
{
  int num = 0;
  std::vector <StateMaterialPointP> pointsArray = getStatePoints();
  int size = pointsArray.size();
  
  for (int i = 0; i < size; i++)
  {
    StateMaterialPointP cur_point = pointsArray[i];
    if (cur_point->getVelocityBoundaryFlag() == true)
    {
      num += 1;
      std::cout << "num " << num;
      std::cout << " (" << cur_point->x() << ", " << cur_point->y() << ", " << cur_point->z() << ")" << std::endl;
      std::cout << "Velocity= (" << cur_point->getVelNew().x() << ", " << cur_point->getVelNew().y() << ", " << cur_point->getVelNew().z() << ")" << std::endl;
      std::cout << std::endl;
    }
  }

//  std::cout << "num= " << num << std::endl;

}




void
PeriBox::fillBoundaryPointsFields()
{

//  Vector3D velocity_vec(0.0, 0.0, 0.0);

  std::vector <StateMaterialPointP> pointsArray = getStatePoints();

  std::vector <StateMaterialPointP> pointsPositive;
  std::vector <StateMaterialPointP> pointsNegative;

  pointsPositive.clear();
  pointsNegative.clear();



  ForceSurface  velocity_surf = getVelocityBoundarySurface();  
//  Vector3D velocity = getBoundaryVelocity();

  int size = pointsArray.size();

   
        
  switch (velocity_surf)
  {
    case right:
    {
      for (int i = 0; i < size; i++)
      { 
        StateMaterialPointP cur_point = pointsArray[i]; 
        if (cur_point->x() == getMaxPoint().x())
        {
          pointsPositive.emplace_back(cur_point);
        }        
      }
      pointsNegative.clear();
      setPointsPositiveField(pointsPositive);
      setPointsNegativeField(pointsNegative);
    }
    break;

    case left:
    {
      for (int i = 0; i < size; i++)
      { 
        StateMaterialPointP cur_point = pointsArray[i]; 
        if (cur_point->x() == getMinPoint().x())
        {
          pointsPositive.emplace_back(cur_point);
        }        
      }
      pointsNegative.clear();
      setPointsPositiveField(pointsPositive);
      setPointsNegativeField(pointsNegative);
    }
    break;


    case up:
    {
      for (int i = 0; i < size; i++)
      { 
        StateMaterialPointP cur_point = pointsArray[i]; 
        if (cur_point->y() == getMaxPoint().y())
        {
          pointsPositive.emplace_back(cur_point);
        }        
      }
      pointsNegative.clear();
      setPointsPositiveField(pointsPositive);
      setPointsNegativeField(pointsNegative);
    }
    break;


    case down:
    {
      for (int i = 0; i < size; i++)
      { 
        StateMaterialPointP cur_point = pointsArray[i]; 
        if (cur_point->y() == getMinPoint().y())
        {
          pointsPositive.emplace_back(cur_point);
        }        
      }
      pointsNegative.clear();
      setPointsPositiveField(pointsPositive);
      setPointsNegativeField(pointsNegative);
    }
    break;


    case front:
    {
      for (int i = 0; i < size; i++)
      { 
        StateMaterialPointP cur_point = pointsArray[i]; 
        if (cur_point->z() == getMaxPoint().z())
        {
          pointsPositive.emplace_back(cur_point);
        }        
      }
      pointsNegative.clear();
      setPointsPositiveField(pointsPositive);
      setPointsNegativeField(pointsNegative);
    }
    break;


    case back:
    {
      for (int i = 0; i < size; i++)
      { 
        StateMaterialPointP cur_point = pointsArray[i]; 
        if (cur_point->z() == getMinPoint().z())
        {
          pointsPositive.emplace_back(cur_point);
        }        
      }
      pointsNegative.clear();
      setPointsPositiveField(pointsPositive);
      setPointsNegativeField(pointsNegative);
    }
    break;

    case upDownX:
    {
      for (int i = 0; i < size; i++)
      { 
        StateMaterialPointP cur_point = pointsArray[i]; 
        if (cur_point->x() == getMaxPoint().x())
        {
          pointsPositive.emplace_back(cur_point);
        }
        if (cur_point->x() == getMinPoint().x())
        {
          pointsNegative.emplace_back(cur_point);
        }
        
      }
      setPointsPositiveField(pointsPositive);
      setPointsNegativeField(pointsNegative);
    }
    break;


    case upDownY:
    {
      for (int i = 0; i < size; i++)
      { 
        StateMaterialPointP cur_point = pointsArray[i]; 
        if (cur_point->y() == getMaxPoint().y())
        {
          pointsPositive.emplace_back(cur_point);
        }
        if (cur_point->y() == getMinPoint().y())
        {
          pointsNegative.emplace_back(cur_point);
        }
        
      }
      setPointsPositiveField(pointsPositive);
      setPointsNegativeField(pointsNegative);
    }
    break;

 
    case upDownZ:
    {
      for (int i = 0; i < size; i++)
      { 
        StateMaterialPointP cur_point = pointsArray[i]; 
        if (cur_point->z() == getMaxPoint().z())
        {
          pointsPositive.emplace_back(cur_point);
        }
        if (cur_point->z() == getMinPoint().z())
        {
          pointsNegative.emplace_back(cur_point);
        }
        
      }
      setPointsPositiveField(pointsPositive);
      setPointsNegativeField(pointsNegative);
    }
    break;


    case tear:
    {

      int pointNum = 1055;

      int nnn = 2;
      int mmm = 3;
      for (int i = 0; i < size; i++)
      { 
        StateMaterialPointP cur_point = pointsArray[i]; 
        double xxx = cur_point->x();
        double zzz = cur_point->z();

        double xLength = getMaxPoint().x() - getMinPoint().x();
        double zLength = getMaxPoint().z() - getMinPoint().z();

        if (cur_point->getID() == pointNum)
            std::cout << "xLength= " << xLength << std::endl;

        double dx = xLength/(getNumPointsX() - 1);
        double dz = zLength/(getNumPointsZ() - 1);

        if (cur_point->getID() == pointNum)
            std::cout << "dx= " << dx << " dz= " << dz << std::endl;

        double midXposition = getMinPoint().x() + xLength/2;

        if (cur_point->getID() == pointNum)
            std::cout << "midXposition= " << midXposition << std::endl;

        double min = midXposition - nnn*dx;

        if (cur_point->getID() == pointNum)
            std::cout << "min= " << min << std::endl;

        double max = midXposition + nnn*dx;

        if (cur_point->getID() == pointNum)
            std::cout << "max= " << max << std::endl;

        double minZ = getMaxPoint().z() - mmm*dz;

        if (cur_point->getID() == pointNum)
            std::cout << "minZ= " << minZ << std::endl;            


        bool midRight = (xxx <= max) && (xxx >= midXposition);
        bool midLeft =  (xxx < midXposition) && (xxx >= min);
        bool midDown =  (zzz <= getMaxPoint().z()) && (zzz >= minZ);

        if (/*(zzz == getMaxPoint().z())*/midDown && midRight)
        {
          pointsNegative.emplace_back(cur_point);
        }            
        else if (/*(zzz == getMaxPoint().z())*/midDown && midLeft)
        {
          pointsPositive.emplace_back(cur_point);
        }
      }
      setPointsPositiveField(pointsPositive);
      setPointsNegativeField(pointsNegative);
   }
   break; 




   case none:
//        cur_point->surfaceForce(force_vec);
   break;


   default:
   std::cout <<
   " The surface can be right, left, up, down, front, back, or none which means we do not have any traction force."
   << std::endl;
   break;

 } // end of switch

//   std::cout << "Num of pos points= " << getPointsPositiveField().size() << ", Num of neg points= " << getPointsNegativeField().size() << std::endl;

//    if (cur_point->getInitialArea() != 0.0)
//    if (cur_point->getVelocityBoundaryFlag() == true)
//    {
//       std::cout << "MinPointBox= " << getMinPoint() << " MinPointBox= " 
//                                    << getMaxPoint() << std::endl;

//       std::cout << "Point " << cur_point->getID() << " "
//                              << cur_point->getPosOld() << " NewPos= " << cur_point->getPosNew() << std::endl
                                                         // << " surfaceOnX= " << surfaceOnX
//                                                   
//                                                   << "Velocity " << cur_point->getVelOld() <<
//                                                  " NewVel "  << cur_point->getVelNew() << std::endl << std::endl;
//    }

}



void
PeriBox::velocityBoundaryCondition()
{
  Vector3D velocity = getBoundaryVelocity();

//  std::cout << "velocityBoundary= " << getBoundaryVelocity() << std::endl;

  std::vector <StateMaterialPointP>  pointsPositive = getPointsPositiveField();
  int posSize = pointsPositive.size();
  for (int i = 0; i < posSize; i++)
  {
    StateMaterialPointP cur_point = pointsPositive[i];
    cur_point->setVelOld(velocity);
    cur_point->setVelMid(velocity);
    cur_point->setVelNew(velocity);
    cur_point->setVelocityBoundaryFlag(true);
  }
  
  std::vector <StateMaterialPointP>  pointsNegative = getPointsNegativeField();
  int negSize = pointsNegative.size();
  for (int j = 0; j < negSize; j++)
  {
    StateMaterialPointP cur_particle = pointsNegative[j];
    cur_particle->setVelOld(velocity*(-1));
    cur_particle->setVelMid(velocity*(-1));
    cur_particle->setVelNew(velocity*(-1));
    cur_particle->setVelocityBoundaryFlag(true);
  }  

  setPointsPositiveField(pointsPositive);
  setPointsNegativeField(pointsNegative);

}



void
PeriBox::findZeroEnergyConstantTwo()
{
  const double c = 5*pow(10,13);
  double E = getYoung();
  double dx = (getMaxPoint().x()-getMinPoint().x())/(getNumPointsX()-1); 
  setZeroEnergyConstantTwo((E*c)/(dx*dx));
//  std::cout << "zeroEnergyConstantTwo= " << getZeroEnergyConstantTwo() << std::endl;


}


void
PeriBox::calcualtePmbConstantSpring()
{

  double pi = M_PI;
  double E = getYoung();
  double nue = getPoission();
  double delta = getHorizon();
  double k = E/(3*(1 - 2*nue));
  setPmbConstantSpring((18.0*k)/(pi*pow(delta,4)));

} 

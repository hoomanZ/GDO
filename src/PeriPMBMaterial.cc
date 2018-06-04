#include <cmath>
#include <math.h>
#include <vector>

#include <math.h>

#include <PeriPMBMaterial.h>
#include <Body.h>
#include <PeriBondPArray.h>
#include <PeriMaterialPointP.h>
#include <PeriMaterialPoint.h>
#include <PeriBondP.h>
#include <PeriBond.h>
//#include <GeometryMath/Vector3D.h>
//#include <GeometryMath/Point3D.h>

using namespace FiniteElement;



PeriPMBMaterial::PeriPMBMaterial()
 : d_density(2440), d_horizon(1.0), d_fracture_energy(15000), d_youngs_modulus(72000000)//, d_body(new Body()) 
{
  calculateCriticalStretch();
  calculateMicromodulusFunction();
}


PeriPMBMaterial::PeriPMBMaterial(double fractureEnergy, double horizon)
 : d_density(2440), d_youngs_modulus(72000000), d_fracture_energy(fractureEnergy)
 , d_horizon(horizon)//, d_body(new Body()) 
{
  calculateCriticalStretch();
  calculateMicromodulusFunction();
}




PeriPMBMaterial::PeriPMBMaterial(double density, double horizon, double fractureEnergy, double youngModulus)
 : d_density(density), d_horizon(horizon), d_fracture_energy(fractureEnergy)
 , d_youngs_modulus(youngModulus)//, d_body(new Body()) 
{
  calculateCriticalStretch();
  calculateMicromodulusFunction();
}


PeriPMBMaterial::PeriPMBMaterial(double density, double horizon, double youngModulus, bool bondBased)
 : d_density(density), d_horizon(horizon), d_fracture_energy(15000), d_youngs_modulus(youngModulus)
{
  if (bondBased)
  {
    calculateCriticalStretch();
    calculateMicromodulusFunction();
  }

}

//PeriPMBMaterial::PeriPMBMaterial(double density, double horizon, double fractureEnergy,
//                                 double youngModulus, const BodySP& body)
// : d_density(density), d_horizon(horizon), d_fracture_energy(fractureEnergy), d_youngs_modulus(youngModulus)  
//{
//  setBody(body);
//  calculateCriticalStretch();
//  calculateMicromodulusFunction();
//}



void
PeriPMBMaterial::calculateCriticalStretch()
{
//  d_critical_stretch = std::sqrt((5*d_fracture_energy)/(6*d_youngs_modulus*d_horizon));
  d_critical_stretch = 0.01;
  std::cout << "critical stretch is s0= " << getCriticalStretch() << std::endl;
}



void
PeriPMBMaterial::calculateMicromodulusFunction()
{
  double c = (12*d_youngs_modulus)/(M_PI*std::pow(d_horizon, 4));
  setMicromodulusFunction(c);
  std::cout << "micromodulus function is c= " << getMicromodulusFunction() << std::endl;
}



void
PeriPMBMaterial::calculatePairwiseForceOfEachBond(const std::vector<PeriMaterialPointP>& pointsArray)
{
////  std::cout << "+++++++++++++++++++calculatePairwiseForceOfEachBond PeriPMB++++++++++++++++++++" << std::endl;
//  BodySP body = getBody();
//  std::vector<PeriMaterialPointP> pointsArray = body->getMaterialPeriPoints();
  int numBond = 0;
  int size = pointsArray.size();
////  std::cout << "PointsArraySize= " << size << std::endl;
  for (int i = 0; i < size; i++)
  { 
    PeriMaterialPointP cur_point = pointsArray[i];
    PeriBondPArray bondsArray = cur_point->getHorizonBonds();     
    int sizeBond = bondsArray.size();
    for (int j = 0; j < sizeBond; j++)
    {
      PeriBondP bond = bondsArray[j];
//      surfaceCorrectionFactor(bond); // find surface correction factor
      bool broken = bond->getBrokenBond();
//      if (cur_point->getID() == 1)
//         std::cout << "broken?= " << broken << std::endl;
      if (!broken)
      {      
        PeriMaterialPointP  center_point = bond->getCenterPoint();
        PeriMaterialPointP  second_point = bond->getSecondPoint();    
        Vector3D eta = second_point->getDisp() - center_point->getDisp();
        Vector3D sai = second_point->getPosOld() - center_point->getPosOld();
        double relative_elongation = ((eta + sai).length() - sai.length())/(sai.length());
//        if (relative_elongation != 0.0)
//        if (bond->getSurfaceCorrectionFactor() != 0.0)
//          std::cout << "center ID= " << center_point->getID() 
//                    << " second ID= " << second_point->getID()
//                   << " correction= " << bond->getSurfaceCorrectionFactor() << std::endl;
//                    << " relative elongation= " << relative_elongation << std::endl; 
        if (relative_elongation <= d_critical_stretch)
        {
          Vector3D unit_force = (eta + sai)/((eta + sai).length());
          pointsArray[i]->getHorizonBonds()[j]
                        ->setPairwiseForce((unit_force)*(relative_elongation*d_micromodulus_function));
          if ((center_point->getID() == size/4 + 1))// && (second_point->getID() == 100))
          {
            numBond ++;
/*            std::cout << "*******************numBond= " << numBond << "***********************" << std::endl;
            std::cout << "centBodyForce= " << center_point->getBodyForce() << std::endl;
            std::cout << "criticalStretch= " << d_critical_stretch;
            std::cout << " micromodulusFunction= " << d_micromodulus_function << std::endl;           
            std::cout << "centID= " << center_point->getID() << " secID()= " << second_point->getID() << std::endl;
            std::cout << "posCent= " << center_point->getPosOld() << std::endl;
            std::cout << "posSec= " << second_point->getPosOld() << std::endl;
            std::cout << "dispCent= " << center_point->getDisp() << std::endl; 
            std::cout << "dispSec= " << second_point->getDisp() << std::endl;
            std::cout << "sai= " << sai << " eta= " << eta << std::endl;
            std::cout << "s= " << relative_elongation << std::endl;
            std::cout << "unit_force= " << unit_force << std::endl;
            std::cout << "pairwiseForce= " << bond->getPairwiseForce() << std::endl;
            std::cout << "bondVolume= " << bond->getVolume() << std::endl;
            std::cout << "pairwiseForce*bondVolume= " << bond->getPairwiseForce()*bond->getVolume()
                                                      << std::endl << std::endl;*/
          }   
        }
        else
        {
          Vector3D zero(0.0, 0.0, 0.0);
//          bond->setPairwiseForce(zero);
          pointsArray[i]->getHorizonBonds()[j]->setPairwiseForce(zero);
//          bond->setBrokenBond(true);
          pointsArray[i]->getHorizonBonds()[j]->setBrokenBond(true);
//          std::cout << pointsArray[i]->getHorizonBonds()[j]->getBrokenBond() << std::endl;

        }
      } 
      
    }
    cur_point->calculateLocalDamage();
    

  }

}



Vector3D
PeriPMBMaterial::calculatePairwiseForce(Point3D pos_second, Point3D pos_center,
                                        Vector3D disp_second, Vector3D disp_center, bool broken)
{
  Vector3D zero(0.0, 0.0, 0.0);
  if (broken)
  {
//    return zero;
  }
  
  else
  {
    Vector3D eta = disp_second - disp_center;
    Vector3D sai = pos_second - pos_center;
    double relative_elongation = ((eta + sai).length() - sai.length())/(sai.length());
  //  if (relative_elongation != 0.0)
  //    std::cout << "center ID= " << center_point->getID() 
  //              << " second ID= " << second_point->getID() 
  //              << " relative elongation= " << relative_elongation << std::endl;
//    std::cout << (relative_elongation <= d_critical_stretch) <<
//               "  relative_elongation= " << relative_elongation << std::endl;  
//    if (relative_elongation <= d_critical_stretch)
//    {

      Vector3D unit_force = (eta + sai)/((eta + sai).length());
      return ((unit_force)*(relative_elongation*d_micromodulus_function));

//    }
//    else
//    {
//      return zero;
//    }
  }
}



Vector3D
PeriPMBMaterial::coefficientOfPoint(PeriMaterialPointP point)
{
  Vector3D force_x(0.0, 0.0, 0.0);
  Vector3D pos_force_x(0.0, 0.0, 0.0);
  Vector3D neg_force_x(0.0, 0.0, 0.0);

  Vector3D force_y(0.0, 0.0, 0.0);
  Vector3D pos_force_y(0.0, 0.0, 0.0);
  Vector3D neg_force_y(0.0, 0.0, 0.0);

  Vector3D force_z(0.0, 0.0, 0.0);
  Vector3D pos_force_z(0.0, 0.0, 0.0);
  Vector3D neg_force_z(0.0, 0.0, 0.0);

  PeriBondPArray bondsArray = point->getHorizonBonds();     
  int sizeBond = bondsArray.size();

  double delta_infinity = getHorizon()/50;

  for (int j = 0; j < sizeBond; j++)
  {
    PeriBondP bond = bondsArray[j];
    bool broken = bond->getBrokenBond();
//    std::cout << "broken= " << broken << std::endl;  
    PeriMaterialPointP  center_point = bond->getCenterPoint();
    PeriMaterialPointP  second_point = bond->getSecondPoint();
    double volume = second_point->getInitialVolume(); 
    if (second_point->x() >= center_point->x())
      pos_force_x = pos_force_x + calculatePairwiseForce(second_point->getPosOld(), center_point->getPosOld(),
                                                         second_point->getUniaxialTensionX(),
                                                         center_point->getUniaxialTensionX(), broken)*volume;
    else 
      neg_force_x = neg_force_x + calculatePairwiseForce(second_point->getPosOld(), center_point->getPosOld(),
                                                         second_point->getUniaxialTensionX(),
                                                         center_point->getUniaxialTensionX(), broken)*volume;

    if (second_point->y() >= center_point->y())
      pos_force_y = pos_force_y + calculatePairwiseForce(second_point->getPosOld(), center_point->getPosOld(),
                                                         second_point->getUniaxialTensionY(),
                                                         center_point->getUniaxialTensionY(), broken)*volume;
    else 
      neg_force_y = neg_force_y + calculatePairwiseForce(second_point->getPosOld(), center_point->getPosOld(),
                                                         second_point->getUniaxialTensionY(),
                                                         center_point->getUniaxialTensionY(), broken)*volume;

    if (second_point->z() >= center_point->z())
      pos_force_z = pos_force_z + calculatePairwiseForce(second_point->getPosOld(), center_point->getPosOld(),
                                                         second_point->getUniaxialTensionZ(),
                                                         center_point->getUniaxialTensionZ(), broken)*volume;
    else 
      neg_force_z = neg_force_z + calculatePairwiseForce(second_point->getPosOld(), center_point->getPosOld(),
                                                         second_point->getUniaxialTensionZ(),
                                                         center_point->getUniaxialTensionZ(), broken)*volume;    
  }
  
  force_x = pos_force_x - neg_force_x;
  force_y = pos_force_y - neg_force_y;
  force_z = pos_force_z - neg_force_z;

  

  Vector3D g(force_x.x(),//*pow(delta_infinity, 3),
             force_y.y(),//*pow(delta_infinity, 3),
             force_z.z());//*pow(delta_infinity, 3));
  
  
  double g_infinity = 0.0;
  Vector3D zero(0.0, 0.0, 0.0);
  Point3D point_zero(0.0, 0.0, 0.0);
  for (int i = 0; i < 51; i++)
    for (int j = 0; j < 51; j++)
      for (int k = 0; k < 51; k++)
      {
//        Point3D sec_point(i*delta_infinity, j*delta_infinity, k*delta_infinity);
        Point3D sec_point((i + 0.5)*delta_infinity, (j + 0.5)*delta_infinity, (k + 0.5)*delta_infinity);
        Vector3D disp(0.001*sec_point.x(), -0.00025*sec_point.y(), -0.00025*sec_point.z());
        if (sec_point != point_zero)
        { 
          g_infinity = g_infinity + 
            calculatePairwiseForce(sec_point, point_zero, disp, zero, false).x()
                                   *delta_infinity*delta_infinity*delta_infinity;
//          if (point->getID() == 1000)
//          {
//            std::cout << "PairwiseForce= " << calculatePairwiseForce(sec_point, point_zero, disp, zero, false) 
//                                           << std::endl;
//            std::cout << "g_infinity= " << g_infinity << std::endl;
//          }
        }
      }
  g_infinity *= 2;
  Vector3D g_c(g_infinity/g.x(), g_infinity/g.y(), g_infinity/g.z());
//  if (point->getID() == 1000)
//  {
//    std::cout << "force_x= " << force_x << std::endl;
//    std::cout << "force_y= " << force_y << std::endl;
//    std::cout << "force_z= " << force_z << std::endl;
//    std::cout << "g= " << g << std::endl;
//    std::cout << "g_infinity= " << g_infinity << std::endl;
//    std::cout << "g_c= " << g_c << std::endl << std::endl; 
//  } 
//  std::cout << point->getID() << " (" << g_c.x() << " " << g_c.y() << " " << g_c.z() << ")" << std::endl;
//  std::cout << point->getID() << " (" << g.x() << " " << g.y() << " " << g.z() << ")" << std::endl;  

  return g_c;
   
}

void
PeriPMBMaterial::surfaceCorrectionFactor(PeriBondP bond)
{
  PeriMaterialPointP  center_point = bond->getCenterPoint();
  PeriMaterialPointP  second_point = bond->getSecondPoint();   
  Vector3D g_c = (coefficientOfPoint(center_point) + coefficientOfPoint(second_point))*(0.5);
  Vector3D unit_vector = second_point->getPosOld() - center_point->getPosOld();
  unit_vector = unit_vector/unit_vector.length();
  double g_bond = (unit_vector.x()/g_c.x())*(unit_vector.x()/g_c.x()) +
                  (unit_vector.y()/g_c.y())*(unit_vector.y()/g_c.y()) +
                  (unit_vector.z()/g_c.z())*(unit_vector.z()/g_c.z());
  g_bond = pow(g_bond, -0.5);
  bond->setSurfaceCorrectionFactor(g_bond);

//        if ((bond->getSurfaceCorrectionFactor() != 0.0))// && (center_point->getID() == 1000))
//          std::cout << "center ID= " << center_point->getID() 
//                    << " second ID= " << second_point->getID()
//                    << " correction= " << bond->getSurfaceCorrectionFactor() << std::endl
//                    << "g_c= (" << g_c.x() << " " << g_c.y() << " " << g_c.z() << ")" << std::endl;
//                    << "g= (" << g.x() << " " << g.y() << " " << g.z() << ")" << std::endl;
//  return g_bond;

}

#include <Peridynamics.h>
#include <PeriPMBMaterial.h>
#include <PeriMaterialPointP.h>
#include <PeriMaterialPoint.h>
#include <PeriBondP.h>
#include <PeriBond.h>
#include <StateBondP.h>
#include <StateBond.h>
#include <StateMaterialPointP.h>
#include <StateMaterialPoint.h>
#include <Body.h>
#include <PeriBox.h>
#include <fstream>

#include <StateMaterialPointPArray.h>
#include <StateMaterialPointP.h>
 

using namespace FiniteElement;


Peridynamics::Peridynamics()
  : d_total_t(10.0), d_del_t(0.001), d_peri_box(new PeriBox()),
    d_alpha(0.0), d_beta(0.0), d_tolerance(0.0)
{

  setBasisPolynomialOrders(Peridynamics::NoPolynomial);
  setGoverningEquationType(Peridynamics::BondBasedPeridynamics);
  run();

}

Peridynamics::Peridynamics(double totalT, double delT)
  : d_total_t(totalT), d_del_t(delT), d_peri_box(new PeriBox())
{

  setBasisPolynomialOrders(Peridynamics::NoPolynomial);
  setGoverningEquationType(Peridynamics::BondBasedPeridynamics);
  run();

}


Peridynamics::Peridynamics(double totalT, double delT, PeriBoxP periBox)
  : d_total_t(totalT), d_del_t(delT), d_peri_box(periBox)
{
  if (d_peri_box->getPeriType() == PeriBox::PeriType::Bond)
  {
    run();
    setBasisPolynomialOrders(Peridynamics::NoPolynomial);
    setGoverningEquationType(Peridynamics::BondBasedPeridynamics);
  }
  else
  {
    setGoverningEquationType(Peridynamics::StateBasedPeridynamics);
    stateRun();
  }

}


Peridynamics::Peridynamics(double totalT, double delT, PeriBoxP periBox, const Matrix3D& initialVelocityGradient)
  : d_total_t(totalT), d_del_t(delT), d_peri_box(periBox)
{
  if (d_peri_box->getPeriType() == PeriBox::PeriType::Bond)
  {
    run();
    setBasisPolynomialOrders(Peridynamics::NoPolynomial);
    setGoverningEquationType(Peridynamics::BondBasedPeridynamics);
  }
  else
  {
    setBasisPolynomialOrders(Peridynamics::One);
    setGoverningEquationType(Peridynamics::StateBasedPeridynamics);
    stateRun(initialVelocityGradient);
  }
}


Peridynamics::Peridynamics(double totalT, double delT, PeriBoxP periBox, PeriCrackP periCrack)
  : d_total_t(totalT), d_del_t(delT), d_peri_box(periBox), d_peri_crack(periCrack)
{
  PeriMaterialPointPArray periPoints = getPeriBox()->getMaterialPeriPoints();
  getPeriCrack()->createInitialCrack(periPoints);
  setBasisPolynomialOrders(Peridynamics::NoPolynomial);
  setGoverningEquationType(Peridynamics::BondBasedPeridynamics);
  run();

}


Peridynamics::Peridynamics(double totalT, double delT, PeriBoxP periBox,
                           const Matrix3D& initialVelocityGradient, PeriCrackP periCrack)
  : d_total_t(totalT), d_del_t(delT), d_peri_box(periBox), d_peri_crack(periCrack)
{
  if (d_peri_box->getPeriType() == PeriBox::PeriType::Bond)
  {
    PeriMaterialPointPArray periPoints = getPeriBox()->getMaterialPeriPoints();
    getPeriCrack()->createInitialCrack(periPoints);
    setBasisPolynomialOrders(Peridynamics::NoPolynomial);
    setGoverningEquationType(Peridynamics::BondBasedPeridynamics);
    run();
  }
  else
  {
    setBasisPolynomialOrders(Peridynamics::One);
    setGoverningEquationType(Peridynamics::StateBasedPeridynamics);
    std::vector<StateMaterialPointP> statePoints = getPeriBox()->getStatePoints();
    getPeriCrack()->createInitialCrack(statePoints);
    stateRun(initialVelocityGradient);
    
  }

}


Peridynamics::Peridynamics(double totalT, double delT, PeriBoxP periBox,
                           const Matrix3D& initialVelocityGradient, 
                           PeriCrackP periCrack, const BasisPolynomialOrders& order)
  : d_total_t(totalT), d_del_t(delT), d_peri_box(periBox), d_peri_crack(periCrack)
{
  if (d_peri_box->getPeriType() == PeriBox::PeriType::Bond)
  {
    PeriMaterialPointPArray periPoints = getPeriBox()->getMaterialPeriPoints();
    getPeriCrack()->createInitialCrack(periPoints);
    setBasisPolynomialOrders(Peridynamics::NoPolynomial);
    setGoverningEquationType(Peridynamics::BondBasedPeridynamics);
    run();
  }
  else
  {
    setBasisPolynomialOrders(order);
    setGoverningEquationType(Peridynamics::StateBasedPeridynamics);
    std::vector<StateMaterialPointP> statePoints = getPeriBox()->getStatePoints();
    getPeriCrack()->createInitialCrack(statePoints);
    stateRun(initialVelocityGradient);
    
  }

}


//*************************************  A non-ordinary state based peridynamics method to model solid material deformation and fracture, Warren et al. 2009 
Peridynamics::Peridynamics(double totalT, double delT, PeriBoxP periBox,
                           const Matrix3D& initialVelocityGradient, 
                           PeriCrackP periCrack, const BasisPolynomialOrders& order, const GoverningEquationType& type)
  : d_total_t(totalT), d_del_t(delT), d_peri_box(periBox), d_peri_crack(periCrack)
{
  if (d_peri_box->getPeriType() == PeriBox::PeriType::Bond)
  {
    PeriMaterialPointPArray periPoints = getPeriBox()->getMaterialPeriPoints();
    getPeriCrack()->createInitialCrack(periPoints);
    setBasisPolynomialOrders(Peridynamics::NoPolynomial);
    setGoverningEquationType(Peridynamics::BondBasedPeridynamics);
    run();
  }
  else
  {
    setBasisPolynomialOrders(order);
    setGoverningEquationType(type);
    std::vector<StateMaterialPointP> statePoints = getPeriBox()->getStatePoints();
    getPeriCrack()->createInitialCrack(statePoints);
    stateRun(initialVelocityGradient);
    
  }

}
//************************************************************************************************************************************************************



Peridynamics::Peridynamics(double totalT, double delT, PeriBoxP periBox, 
                           PeriCrackP periCrack, const BasisPolynomialOrders& order, const GoverningEquationType& type)
  : d_total_t(totalT), d_del_t(delT), d_peri_box(periBox), d_peri_crack(periCrack)
{
  if (d_peri_box->getPeriType() == PeriBox::PeriType::Bond)
  {
    PeriMaterialPointPArray periPoints = getPeriBox()->getMaterialPeriPoints();
    getPeriCrack()->createInitialCrack(periPoints);
    setBasisPolynomialOrders(Peridynamics::NoPolynomial);
    setGoverningEquationType(Peridynamics::BondBasedPeridynamics);
    run();
  }
  else
  {
    setBasisPolynomialOrders(order);
    setGoverningEquationType(type);
    std::vector<StateMaterialPointP> statePoints = getPeriBox()->getStatePoints();
    getPeriCrack()->createInitialCrack(statePoints);

    std::cout << "Salam " << std::endl;

    linearElasticRun();


//    stateRun(initialVelocityGradient);
    
  }

}

//***************************************** Implicit Time Integration ***********************************


Peridynamics::Peridynamics(double totalT, double delT, PeriBoxP periBox, 
                           PeriCrackP periCrack, const BasisPolynomialOrders& order,
                           const GoverningEquationType& type, const TimeIntegration& timeIntegration,
                           double alpha, double beta, double tolerance)
  : d_total_t(totalT), d_del_t(delT), d_peri_box(periBox), d_peri_crack(periCrack)
{
  if (d_peri_box->getPeriType() == PeriBox::PeriType::Bond)
  {
    PeriMaterialPointPArray periPoints = getPeriBox()->getMaterialPeriPoints();
    getPeriCrack()->createInitialCrack(periPoints);
    setBasisPolynomialOrders(Peridynamics::NoPolynomial);
    setGoverningEquationType(Peridynamics::BondBasedPeridynamics);
    run();
  }
  else
  {
    setBasisPolynomialOrders(order);
    setGoverningEquationType(type);
    setTimeIntegration(timeIntegration);
    setAlpha(alpha);
    setBeta(beta);
    setTolerance(tolerance);
//    std::vector<StateMaterialPointP> statePoints = getPeriBox()->getStatePoints();
//    getPeriCrack()->createInitialCrack(statePoints);

    std::cout << "Salam " << std::endl;

//    testDivergance();

//    testDerivative();

    linearElasticRun();


//    stateRun(initialVelocityGradient);
    
  }

}


//**********************************************************************************************************



//Peridynamics::Peridynamics(PeriPMBMaterialP pmbMaterial)
//{
//  d_pmb_material->setDensity(pmbMaterial->getDensity());
//  d_pmb_material->setHorizon(pmbMaterial->getHorizon());
//  d_pmb_material->setFractureEnergy(pmbMaterial->getFractureEnergy());
//  d_pmb_material->setYoungsModulus(pmbMaterial->getYoungsModulus());
//  d_pmb_material->setCriticalStretch(pmbMaterial->getCriticalStretch());
//  d_pmb_material->setMicromodulusFunction(pmbMaterial->getMicromodulusFunction());
//  d_pmb_material->setBody(pmbMaterial->getBody());

//  d_body = d_pmb_material->getBody();
//  d_horizon = d_pmb_material->getHorizon();

//  createGeometry();
//  createBonds();
//  d_pmb_material->calculatePairwiseForceOfEachBond();

//}



Peridynamics::~Peridynamics()
{


}


void
Peridynamics::createGeometry()
{
//  BodySP body = d_pmb_material->getBody();
//  body->setPeriPoints();
//  d_pmb_material->getBody()->setMaterialPeriPoints(body->getMaterialPeriPoints());
//  d_body->setPeriPoints();

}


void
Peridynamics::createBonds()
{
//  double horizon = d_pmb_material->getHorizon();
//  d_body->createPeriBondsOfEachPeriPoints(d_horizon);

}


void
Peridynamics::findStressOfPoints()
{
  double dx = (d_peri_box->getMaxPoint().x() - d_peri_box->getMinPoint().x())/(d_peri_box->getNumPointsX() - 1);
  double dy = (d_peri_box->getMaxPoint().y() - d_peri_box->getMinPoint().y())/(d_peri_box->getNumPointsY() - 1);
  double dz = (d_peri_box->getMaxPoint().z() - d_peri_box->getMinPoint().z())/(d_peri_box->getNumPointsZ() - 1);
  Point3D delX(0.0, 0.5*dy, 0.5*dz);
  Point3D delY(0.5*dx, 0.0, 0.5*dz);
  Point3D delZ(0.5*dx, 0.5*dy, 0.0);
  PeriMaterialPointPArray periPoints = getPeriBox()->getMaterialPeriPoints();
  int pointsSize = periPoints.size();
//  std::cout << pointsSize << std::endl;
  for (int i = 0; i < pointsSize; i++)
  {
    PeriMaterialPointP cur_point = periPoints[i];
    Vector3D stress(0.0);
    for (int j = 0; j < pointsSize; j++)
    {
      PeriMaterialPointP point = periPoints[j];
      PeriBondPArray pointBonds = point->getHorizonBonds();
      int bondsSize = pointBonds.size();
      for (int k = 0; k < bondsSize; k++)
      {
        PeriBondP bond = pointBonds[k];
        Vector3D minX = cur_point->getPosOld() - delX;
        Vector3D maxX = cur_point->getPosOld() - delX*(-1);
        double xDiff = bond->getSecondPoint()->x() - bond->getCenterPoint()->x();
        if (intersectBondWithPlane(bond, minX, maxX) && xDiff > 0)
        {
          double fx = (bond->getPairwiseForce().length())*(xDiff/(bond->lengthOfSai()));
          stress.x(stress.x() + fx*pow(dx, 4));                   
        }
        Vector3D minY = cur_point->getPosOld() - delY;
        Vector3D maxY = cur_point->getPosOld() - delY*(-1);
        double yDiff = bond->getSecondPoint()->y() - bond->getCenterPoint()->y();
        if (intersectBondWithPlane(bond, minY, maxY) && yDiff > 0)
        {
          double fy = (bond->getPairwiseForce().length())*(yDiff/(bond->lengthOfSai()));
          stress.y(stress.y() + fy*pow(dy, 4));                   
        }
        Vector3D minZ = cur_point->getPosOld() - delZ;
        Vector3D maxZ = cur_point->getPosOld() - delZ*(-1);
        double zDiff = bond->getSecondPoint()->z() - bond->getCenterPoint()->z();
        if (intersectBondWithPlane(bond, minZ, maxZ) && zDiff > 0)
        {
          double fz = (bond->getPairwiseForce().length())*(zDiff/(bond->lengthOfSai()));
          stress.z(stress.z() + fz*pow(dy, 4));                   
        }       
      }
         
    }
    cur_point->setPeriStress(stress);
//    std::cout << "ID= " << cur_point->getID() << " stress= " << cur_point->getPeriStress() << std::endl;  

  }

}



void
Peridynamics::periPMBBoxVelocityVerlet()
{
////  std::cout << "--------------------------periPMBBoxVelocityVerlet------------------------" << std::endl;
  if (!(d_peri_box->getVelocityBoundarySurface() == PeriBox::ForceSurface::none))
    d_peri_box->applyBoundaryVelocity();

//  d_pmb_material->calculatePairwiseForceOfEachBond();
  d_peri_box->getPMB()->calculatePairwiseForceOfEachBond(d_peri_box->getMaterialPeriPoints());
//  d_body->calculateAccelerationOfEachPeriPoints();
//  d_body->calculateDisplacementAndVelocityOfPointsUsingVelocityVerlet(d_del_t);
  d_peri_box->calculateAccelerationOfEachPeriPoints();
  d_peri_box->calculateDisplacementAndVelocityOfPointsUsingVelocityVerlet(d_del_t);// v(n+1/2) = v(n)+a(n)*delt/2
//  d_pmb_material->calculatePairwiseForceOfEachBond();                            // u(n+1) = u(n) + v(n+1/2)*delt

  if (!(d_peri_box->getVelocityBoundarySurface() == PeriBox::ForceSurface::none))
    d_peri_box->applyBoundaryVelocity();

  d_peri_box->getPMB()->calculatePairwiseForceOfEachBond(d_peri_box->getMaterialPeriPoints()); 
//  d_body->calculateAccelerationOfEachPeriPoints();
  d_peri_box->calculateAccelerationOfEachPeriPoints();
//  std::vector<PeriMaterialPointP> periPointsArray = d_body->getMaterialPeriPoints();
  std::vector<PeriMaterialPointP> periPointsArray = d_peri_box->getMaterialPeriPoints();
  int size = periPointsArray.size();
  for (int i = 0; i < size; i++)
  {  
    PeriMaterialPointP cur_point = periPointsArray[i];
//    if (cur_point->getID() == size/4 + 1)
//    {
//      std::cout << "accelerationNew(n+1)= " << cur_point->getAccele() << std::endl;
//    }
    cur_point->setVelNew(cur_point->getVelMid() + cur_point->getAccele()*(0.5*d_del_t)); // v(n+1) = v(n+1/2)+a(n+1)*delt/2
//    if (cur_point->getID() == size/4 + 1)
//    {
//      std::cout << "velocityNew(n+1)= " << cur_point->getVelNew() << std::endl << std::endl;
//    }

    cur_point->setVelOld(cur_point->getVelNew());
 
  }

  if (!(d_peri_box->getVelocityBoundarySurface() == PeriBox::ForceSurface::none))
    d_peri_box->applyBoundaryVelocity();

  for (int k = 0; k < size; k++)
  {
    PeriMaterialPointP cur_point = periPointsArray[k];
    cur_point->setPosNew(cur_point->getPosOld() + cur_point->getDisp());    

  }
//  d_body->setMaterialPeriPoints(periPointsArray);
  d_peri_box->setMaterialPeriPoints(periPointsArray);

}


void
Peridynamics::run()
{
  int num_iter = 0;
  int number = 0;
  while (num_iter*d_del_t < d_total_t)
  {
    std::cout << std::endl << "num_iter= " << num_iter << std::endl;
    periPMBBoxVelocityVerlet();
    num_iter++;
    
    std::vector<PeriMaterialPointP> periPointsArray = d_peri_box->getMaterialPeriPoints();
    int size = periPointsArray.size();
    for (int i = 0; i < size; i++)
    {  
      PeriMaterialPointP cur_point = periPointsArray[i];
      if (cur_point->getID() == size/4 - 1)
        std::cout << "ID= " << cur_point->getID() << " Disp= (" <<
                               cur_point->getDisp().x() << ", " <<
                               cur_point->getDisp().y() << ", " <<
                               cur_point->getDisp().z() << ")" << " Damage= " << cur_point->getLocalDamage()<< std::endl;
//      if ((cur_point->getLocalDamage() != 1))// && (cur_point->getLocalDamage() < 1))
//        std::cout << "ID= " << cur_point->getID() << " Damage= " <<
//                               cur_point->getLocalDamage() << std::endl;  
    }

//   if ((num_iter == 1) || (num_iter == 20000) || (num_iter == 32000) || (num_iter == 35999))
//   {
//     findStressOfPoints();
//   }
    
   char buffer[2000];
   char *str = getcwd(buffer, 2000);
   std::string currentDirec = std::string(buffer);
    
    if (num_iter % 5000 == 1)
    {
      pointsFile(number, currentDirec);
      number++;            
    }


  }


}


  void
  Peridynamics::pointsFile(int& iter_num, const std::string& currentFolder)
  {
//    std::string pointFolderName = "pointsFolder";
//    int status2;
    int ret2;
//    if (iter_num == 0)
//    {

//      makeOutputFolder(status2, ret2, pointFolderName);
//    }
//    else
//    {
//      std::string direcLocation = currentFolder + "/" + pointFolderName;
//      ret2 = chdir(direcLocation.c_str());
//    }
    std::string fileName = "box_" + std::to_string(iter_num) + ".exdata";
    std::ofstream myfile(fileName);
    myfile << " Group name: box" // << "_" << std::to_string(iter_num) 
                                  << std::endl;
//    myfile << " #Fields=1" << std::endl;
    myfile << " #Fields=2" << std::endl;
    myfile << " 1) coordinates, coordinate, rectangular cartesian, #Components=3" << std::endl;
    myfile << "   x.  Value index= 1, #Derivatives= 0" << std::endl;
    myfile << "   y.  Value index= 2, #Derivatives= 0" << std::endl;
    myfile << "   z.  Value index= 3, #Derivatives= 0" << std::endl;
    myfile << " 2) localDamage, field, rectangular cartesian, #Components=1" << std::endl;
//    myfile << " 2) stressX, field, rectangular cartesian, #Components=1" << std::endl;
    myfile << "   localDamage. Value index= 4, #Derivatives= 0" << std::endl;
//    myfile << "   stressX. Value index= 4, #Derivatives= 0" << std::endl;


    std::vector<PeriMaterialPointP> periPointsArray = d_peri_box->getMaterialPeriPoints();
    int size = periPointsArray.size();
    for (int i = 0; i < size; i++)
    {  
      PeriMaterialPointP cur_point = periPointsArray[i];
      myfile << " Node:         " << cur_point->getID()  << std::endl;
      myfile << "  " << cur_point->getPosNew().x() << "  "
             << cur_point->getPosNew().y() << "  "
             << cur_point->getPosNew().z() << "  "
             << cur_point->getLocalDamage()
//             << cur_point->getPeriStress().x()  
             << std::endl;
    }
    
    myfile.close();
    ret2 = chdir(currentFolder.c_str());

  }


bool
Peridynamics::intersectBondWithPlane(const PeriBondP& bond, const Vector3D& minPoint, const Vector3D& maxPoint)
{
  bool intersect = false;
  PeriMaterialPointP center_point = bond->getCenterPoint();
  PeriMaterialPointP second_point = bond->getSecondPoint();   
  bool intersectX = interval(minPoint.x(), maxPoint.x(), center_point->x(), second_point->x());
  bool intersectY = interval(minPoint.y(), maxPoint.y(), center_point->y(), second_point->y());
  bool intersectZ = interval(minPoint.z(), maxPoint.z(), center_point->z(), second_point->z());


  intersect = intersectX && intersectY && intersectZ; 


  return intersect;
}



bool
Peridynamics::interval(double Min, double Max, double A, double B)
{
  bool intersection = false;
  if (A == B)
    intersection = ((Min <= B) && (B <= Max));
  else
  {
    double left = (Min - B)/(A - B);
    double right = (Max - B)/(A - B);
    intersection = (((left < 0) && (right < 0)) || ((left > 1) && (right > 1)));
    intersection = !(intersection);
  }
  return intersection; 
 
}



//********************************************* State-Based Peridynamics ******************************************


void
Peridynamics::stateLargeRotation()
{
  if (!(d_peri_box->getVelocityBoundarySurface() == PeriBox::ForceSurface::none))
  {
    d_peri_box->applyBoundaryVelocityState();
  }


  

  double E = d_peri_box->getYoung();
  double nue = d_peri_box->getPoission();
  double yield = d_peri_box->getYieldStrength();

  double firstLamme = (E*nue)/((1 + nue)*(1 - 2*nue));
  double secondLamme = E/(2*(1 + nue));


  StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
  int pointsSize = statePoints.size();


  int i = 0;
  
  for (i = 0; i < pointsSize; i++)
  {
    StateMaterialPointP cur_point = statePoints[i];
    cur_point->calculateForceVectorStateOfBonds(d_del_t, firstLamme, secondLamme, yield);
  
  }

  d_peri_box->calculateAccelerationOfEachStatePoints(d_del_t);  // Body.cc



}



void
Peridynamics::stateLargeRotation(const int& numIter, const Matrix3D& initialVelocityGradient)
{
  if (!(d_peri_box->getVelocityBoundarySurface() == PeriBox::ForceSurface::none))
  {
    d_peri_box->applyBoundaryVelocityState();
  }


  

  double E = d_peri_box->getYoung();
  double nue = d_peri_box->getPoission();
  double yield = d_peri_box->getYieldStrength();

  double firstLamme = (E*nue)/((1 + nue)*(1 - 2*nue));
  double secondLamme = E/(2*(1 + nue));

  double energyZeroConstantTwo = d_peri_box->getZeroEnergyConstantTwo();

//  std::cout << "zeroEnergyConstantTwo= " << energyZeroConstantTwo << std::endl;


  StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
  int pointsSize = statePoints.size();


  int i = 0;

  bool flag = false;

  do
  {
    flag = false;

    for (i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
      if (!(cur_point->getLocalDamage() == 0))
      {
        cur_point->LagrangianStrainTensor();
      }

    }


    for (i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
      if (!(cur_point->getLocalDamage() == 0))
      {
        cur_point->checkBrokenBonds(yield, flag);
      }
//      cur_point->calculateStateLocalDamage();
  
    }

  } while (flag);


  for (i = 0; i < pointsSize; i++)
  {
    StateMaterialPointP cur_point = statePoints[i];
    if (!(cur_point->getLocalDamage() == 0))
    {
      cur_point->calculateStateLocalDamage();
    }  

  }


  int orderFlag = 0;  
  if (getBasisPolynomialOrders()==Peridynamics::One)
  {
    orderFlag = 1;
  }
  else if (getBasisPolynomialOrders()==Peridynamics::Two)
  {
    orderFlag = 2;
  }
  else
  {
    orderFlag = 0;
  }


  if (getGoverningEquationType() != Divergance)
  {
    
    for (i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
      if (!(cur_point->getLocalDamage() == 0))
      {
        cur_point->calculateForceVectorStateOfBonds(d_del_t, firstLamme, 
                                                    secondLamme, yield,
                                                    numIter, initialVelocityGradient, orderFlag);
////        cur_point->calculateForceVectorStateOfBonds(d_del_t, firstLamme, 
////                                                    secondLamme, yield,
////                                                    numIter, initialVelocityGradient);
//        cur_point->calculateForceVectorStateOfBonds(d_del_t, firstLamme, 
//                                                    secondLamme, yield,
//                                                    numIter, initialVelocityGradient,
//                                                    energyZeroConstantTwo);

      }
//      else
//      {
//        cur_point->setForceBondsZero();
//        
//      }

    }

    d_peri_box->calculateAccelerationOfEachStatePoints(d_del_t);  // Body.cc
//    d_peri_box->calculateAccelerationOfEachStatePoints(d_del_t, springConstant); // Applying Penalty Approach to control the zero-energy modes

  }

  else
  {
    for (i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
      if (!(cur_point->getLocalDamage() == 0))
      {
        cur_point->calculateFirstPiolaKirchhoffStress(d_del_t, firstLamme, 
                                                      secondLamme, yield,
                                                      numIter, initialVelocityGradient, orderFlag);
      }
    }

    for (i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
      if (!(cur_point->getLocalDamage() == 0))
      {
        cur_point->forceVectorStateOfBonds();
      }
//      else
//      {
//        cur_point->setForceBondsZero();
//        
//      }

    }

    d_peri_box->accelerationOfEachStatePoints(d_del_t);  // Body.cc                
  }

//  double springConstant = d_peri_box->getPmbConstantSpring();

//  std::cout << "springConstant= " << springConstant << std::endl;


//  d_peri_box->accelerationOfEachStatePoints(d_del_t);  // Body.cc
//  d_peri_box->calculateAccelerationOfEachStatePoints(d_del_t, springConstant); // Applying Penalty Approach to control the zero-energy modes

}



void
Peridynamics::stateRun(const Matrix3D& initialVelocityGradient)
{
  StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
  int pointsSize = statePoints.size();

  int i = 0;
 
  for (i = 0; i < pointsSize; i++)
  {
    StateMaterialPointP cur_point = statePoints[i];
////    cur_point->shapeTensor();  //for simpleDeformationTensor
 //   std::cout << "After calculating shapeTensor of point number " << cur_point->getID() << std::endl;
    if (cur_point->getID() == 4229)
    {
      std::cout << "BasisPolynomialOrder is " << getBasisPolynomialOrders() << std::endl;
    }
    cur_point->setOrderFlag(getBasisPolynomialOrders());
    cur_point->deformationGradientUsingPolynomials();
//    std::cout << "After calculating deformationGradient of point number " << cur_point->getID() << std::endl; 
    cur_point->polarDecompositionOfDeformationGradient();
//    std::cout << "After polarDecompositionOfDeformationGradient of point number " << cur_point->getID() << std::endl; 
 
  } 

 

  int num_iter = 0;
  int number = 0;
  while (num_iter*d_del_t < d_total_t)
  {
    std::cout << std::endl << "num_iter= " << num_iter << std::endl;
    
//    int num_initial = 10;
//    int num_initial = 21;  
    int num_initial = 10;

//    if (num_iter == 0)
    if (num_iter < num_initial)
    {
      std::vector<StateMaterialPointP> statePointsArray = d_peri_box->getStatePoints();
      int size = statePointsArray.size();
      int i = 0;
      for (i = 0; i < size; i++)
      {  
        StateMaterialPointP cur_point = statePointsArray[i];
        if (!(cur_point->getLocalDamage() == 0))
        {
//          cur_point->setVelocityGradient(initialVelocityGradient*ramp(num_iter, num_initial));
          cur_point->setVelocityGradient(initialVelocityGradient*triangle(num_iter, num_initial));
        }
      }
    }

    stateLargeRotation(num_iter, initialVelocityGradient);
    num_iter ++;


    std::vector<StateMaterialPointP> statePointsArray = d_peri_box->getStatePoints();
    int size = statePointsArray.size();
    int i = 0;
    for (i = 0; i < size; i++)
    {  
      StateMaterialPointP cur_point = statePointsArray[i];
      if (cur_point->getID() == size/4 - 1)
        std::cout << "ID= " << cur_point->getID() << " Disp= (" <<
                               cur_point->getDisp().x() << ", " <<
                               cur_point->getDisp().y() << ", " <<
                               cur_point->getDisp().z() << ")" << " Damage= " <<
                               cur_point->getLocalDamage()<< std::endl;
//      if ((cur_point->getLocalDamage() != 1))// && (cur_point->getLocalDamage() < 1))
//        std::cout << "ID= " << cur_point->getID() << " Damage= " <<
//                               cur_point->getLocalDamage() << std::endl;  
    }


   char buffer[2000];
   char *str = getcwd(buffer, 2000);
   std::string currentDirec = std::string(buffer);
    
//    if (num_iter % 5 == 1)
    if (num_iter % 10 == 1)
    {
      pointsFileState(number, currentDirec);
      number++;            
    }


  }
}



void
Peridynamics::stateRun()
{
/*  StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
  int pointsSize = statePoints.size();

  int i = 0;
 
  for (i = 0; i < pointsSize; i++)
  {
    StateMaterialPointP cur_point = statePoints[i];
    cur_point->shapeTensor();
 //   std::cout << "After calculating shapeTensor of point number " << cur_point->getID() << std::endl;
    cur_point->deformationGradient();
//    std::cout << "After calculating deformationGradient of point number " << cur_point->getID() << std::endl; 
    cur_point->polarDecompositionOfDeformationGradient();
//    std::cout << "After polarDecompositionOfDeformationGradient of point number " << cur_point->getID() << std::endl; 
 
  } */

 

  int num_iter = 0;
  int number = 0;
  while (num_iter*d_del_t < d_total_t)
  {
    std::cout << std::endl << "num_iter= " << num_iter << std::endl;
    stateLargeRotation();
    num_iter ++;


/*    std::vector<StateMaterialPointP> statePointsArray = d_peri_box->getStatePoints();
    int size = statePointsArray.size();
    int i = 0;
    for (i = 0; i < size; i++)
    {  
      StateMaterialPointP cur_point = statePointsArray[i];
      if (cur_point->getID() == size/4 - 1)
        std::cout << "ID= " << cur_point->getID() << " Disp= (" <<
                               cur_point->getDisp().x() << ", " <<
                               cur_point->getDisp().y() << ", " <<
                               cur_point->getDisp().z() << ")" << " Damage= " <<
                               cur_point->getLocalDamage()<< std::endl;
//      if ((cur_point->getLocalDamage() != 1))// && (cur_point->getLocalDamage() < 1))
//        std::cout << "ID= " << cur_point->getID() << " Damage= " <<
//                               cur_point->getLocalDamage() << std::endl;  
    }


   char buffer[2000];
   char *str = getcwd(buffer, 2000);
   std::string currentDirec = std::string(buffer);
    
    if (num_iter % 500 == 1)
    {
      pointsFileState(number, currentDirec);
      number++;            
    } */


  }
}


void
Peridynamics::polynomialOrderPoint()
{

  int orderFlag = 0;  
  if (getBasisPolynomialOrders()==Peridynamics::One)
  {
    orderFlag = 1;
  }
  else if (getBasisPolynomialOrders()==Peridynamics::Two)
  {
    orderFlag = 2;
  }
  else
  {
    orderFlag = 0;
  }



  StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
  int pointsSize = statePoints.size();

  int i = 0;
 
  for (i = 0; i < pointsSize; i++)
  {
    StateMaterialPointP cur_point = statePoints[i];
    cur_point->setOrderFlag(orderFlag);
  }

}


/*************************

void
Peridynamics::linearElasticRun()
{

  polynomialOrderPoint();
  

  double E = d_peri_box->getYoung();
  double nue = d_peri_box->getPoission();
  double yield = d_peri_box->getYieldStrength();

  double firstLamme = (E*nue)/((1 + nue)*(1 - 2*nue));
  double secondLamme = E/(2*(1 + nue));

  double alpha = getAlpha();
  double beta = getBeta();

  double density = getPeriBox()->getDensity();

  int num_iter = 0;
  int number = 0;



   char buffer[2000];
   char *str = getcwd(buffer, 2000);
   std::string currentDirec = std::string(buffer);


  if (getTimeIntegration() == Implicit)
  {
    StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
    int pointsSize = statePoints.size();
    for (int i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
      cur_point->updateBroydenMatrix(alpha);
    }

  }


  if (!(d_peri_box->getVelocityBoundarySurface() == PeriBox::ForceSurface::none))
  {
    d_peri_box->fillBoundaryPointsFields();
  }


  while (num_iter*d_del_t < d_total_t)
  {

    if (!(d_peri_box->getVelocityBoundarySurface() == PeriBox::ForceSurface::none))
    {
//      d_peri_box->applyBoundaryVelocityState();
      d_peri_box->velocityBoundaryCondition();
//      d_peri_box->showBoundaryVelocityPoints();
    }

//    d_peri_box->showBoundaryVelocityPoints();

    std::cout << std::endl << "num_iter= " << num_iter << std::endl;

//    damageModelingUsingCriticalStrain(yield);
    damageModelingUsingCriticalStrain(yield, num_iter);


    StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
    int pointsSize = statePoints.size();

    int i = 0;

    for (i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
//      if (cur_point->getID() == 1055)
//         std::cout << cur_point->getDisp() << std::endl;  
      cur_point->updateNumIter(num_iter);               
    }




    if (getTimeIntegration() == Explicit)
    {
 
      for (i = 0; i < pointsSize; i++)
      {
        StateMaterialPointP cur_point = statePoints[i];
  //      if (cur_point->getID() == 1055)
  //         std::cout << cur_point->getDisp() << std::endl;  
        cur_point->cauchyStress(firstLamme, secondLamme);               
      }



      for (i = 0; i < pointsSize; i++)
      {
        StateMaterialPointP cur_point = statePoints[i];
        if (getGoverningEquationType() == Divergance)
        {  
          cur_point->updatePosition(density, d_del_t);
        }
        if (getGoverningEquationType() == StateBasedPeridynamics)
        {
          cur_point->updatePositionUsingStateBasedPeridynamics(density, d_del_t);
        }

      }
                    
    }

    else
    {
 
      double tolerance = getTolerance();
      for (i = 0; i < pointsSize; i++)
      {
        int governingType = 0;
        if (getGoverningEquationType() == StateBasedPeridynamics) governingType = 1;
        StateMaterialPointP cur_point = statePoints[i];
        int id = 3762; 


        cur_point->updateAccelerationDisplacementImplicitly(firstLamme, secondLamme,
                                                            alpha, beta, density, d_del_t, 
                                                            tolerance, governingType);


//        cur_point->updateAccelerationDisplacementImplicitly(firstLamme, secondLamme,
//                                                            alpha, beta, density, d_del_t, 
//                                                            tolerance, governingType,
//                                                            num_iter, currentDirec);



        Vector3D averageAccele = (cur_point->getAcceleNew() + cur_point->getAccele())*(d_del_t/2.0);
        cur_point->setVelNew(cur_point->getVelOld() + averageAccele);

//        if (cur_point->getID() == 1050)
//        {
//          std::cout << "accele= " << cur_point->getAccele() << std::endl;
//          std::cout << "acceleNew= " << cur_point->getAcceleNew() << std::endl;
//          std::cout << "velOld= " << cur_point->getVelOld() << std::endl;
//          std::cout << "velNew= " << cur_point->getVelNew() << std::endl;
//        }

        cur_point->setDispOld(cur_point->getDispNew());
        cur_point->setVelOld(cur_point->getVelNew());
        cur_point->setAccele(cur_point->getAcceleNew());
        cur_point->setPosNew(cur_point->getPosOld() + cur_point->getDisp());

      }

    }




    std::vector<StateMaterialPointP> statePointsArray = d_peri_box->getStatePoints();
    int size = statePointsArray.size();
    int id = size/4 - 1;
    for (i = 0; i < size; i++)
    {  
      StateMaterialPointP cur_point = statePointsArray[i];
      if (cur_point->getID() == id) 
//      size/4 - 1)
      {
        std::cout << "ID= " << cur_point->getID() << " Disp= (" <<
                               cur_point->getDisp().x() << ", " <<
                               cur_point->getDisp().y() << ", " <<
                               cur_point->getDisp().z() << ")" << " Damage= " <<
                               cur_point->getLocalDamage()<< std::endl;
//        std::cout  << "DispNew= " << cur_point->getDispNew() << std::endl;
//      if ((cur_point->getLocalDamage() != 1))// && (cur_point->getLocalDamage() < 1))
//        std::cout << "ID= " << cur_point->getID() << " Damage= " <<
//                               cur_point->getLocalDamage() << std::endl;
      }  
    }
    
//    if (num_iter % 5 == 1)
    if (num_iter % 10 == 1)
    {
         pointsFileState(number, currentDirec);

//      if (num_iter > 80 && num_iter <90)
//      {
//         printDetails(number, currentDirec);    
//      }
      number++;            
    }

    num_iter++;  

  }

}


************************/




void
Peridynamics::linearElasticRun()
{

  polynomialOrderPoint();
  

  double E = d_peri_box->getYoung();
  double nue = d_peri_box->getPoission();
  double yield = d_peri_box->getYieldStrength();

  double firstLamme = (E*nue)/((1 + nue)*(1 - 2*nue));
  double secondLamme = E/(2*(1 + nue));

  double alpha = getAlpha();
  double beta = getBeta();

  double density = getPeriBox()->getDensity();

  int num_iter = 0;
  int number = 0;

  int id = 3909;
  int id2 = 4184;

   char buffer[2000];
   char *str = getcwd(buffer, 2000);
   std::string currentDirec = std::string(buffer);


  if (getTimeIntegration() == Implicit)
  {
    StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
    int pointsSize = statePoints.size();
    for (int i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
      cur_point->updateBroydenMatrix(alpha);
    }

  }


  if (!(d_peri_box->getVelocityBoundarySurface() == PeriBox::ForceSurface::none))
  {
    d_peri_box->fillBoundaryPointsFields();
  }


  while (num_iter*d_del_t < d_total_t)
  {

    if (!(d_peri_box->getVelocityBoundarySurface() == PeriBox::ForceSurface::none))
    {
//      d_peri_box->applyBoundaryVelocityState();
      d_peri_box->velocityBoundaryCondition();
//      d_peri_box->showBoundaryVelocityPoints();
    }

//    d_peri_box->showBoundaryVelocityPoints();

    std::cout << std::endl << "num_iter= " << num_iter << std::endl;

//    damageModelingUsingCriticalStrain(yield);
    damageModelingUsingCriticalStrain(yield, num_iter);


    StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
    int pointsSize = statePoints.size();

    int i = 0;

    for (i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
//      if (cur_point->getID() == 1055)
//         std::cout << cur_point->getDisp() << std::endl;  
      cur_point->updateNumIter(num_iter);               
    }




    if (getTimeIntegration() == Explicit)
    {
 
      for (i = 0; i < pointsSize; i++)
      {
        StateMaterialPointP cur_point = statePoints[i];
  //      if (cur_point->getID() == 1055)
  //         std::cout << cur_point->getDisp() << std::endl;  
        cur_point->cauchyStress(firstLamme, secondLamme);               
      }



      for (i = 0; i < pointsSize; i++)
      {
        StateMaterialPointP cur_point = statePoints[i];
        if (getGoverningEquationType() == Divergance)
        {  
          cur_point->updatePosition(density, d_del_t);
        }
        if (getGoverningEquationType() == StateBasedPeridynamics)
        {
          cur_point->updatePositionUsingStateBasedPeridynamics(density, d_del_t);
        }

      }
                    
    }

    else
    {
 
//      bool num_condition = (num_iter > 80 && num_iter < 90);
      bool num_condition = (num_iter < 490);
//      bool num_condition = (num_iter > 0 && num_iter < 4);
//      std::string fileName = "details_" + std::to_string(num_iter) + ".com";
      int ret;


      double tolerance = getTolerance();
      for (i = 0; i < pointsSize; i++)
      {
        int governingType = 0;
        if (getGoverningEquationType() == StateBasedPeridynamics) governingType = 1;
        StateMaterialPointP cur_point = statePoints[i];
 
        std::string fileName = "details_" + std::to_string(cur_point->getID()) + "_" + std::to_string(num_iter) + ".com";

        bool writing_condition = ((cur_point->getID() == id) && num_condition);
        writing_condition = writing_condition || (cur_point->getID() == id2);

//        cur_point->updateAccelerationDisplacementImplicitly(firstLamme, secondLamme,
//                                                            alpha, beta, density, d_del_t, 
//                                                            tolerance, governingType);


        cur_point->updateAccelerationDisplacementImplicitly(firstLamme, secondLamme,
                                                            alpha, beta, density, d_del_t, 
                                                            tolerance, governingType,
                                                            num_iter, currentDirec);



        Vector3D averageAccele = (cur_point->getAcceleNew() + cur_point->getAccele())*(d_del_t/2.0);
        cur_point->setVelNew(cur_point->getVelOld() + averageAccele);

//        if (cur_point->getID() == 1050)
//        {
//          std::cout << "accele= " << cur_point->getAccele() << std::endl;
//          std::cout << "acceleNew= " << cur_point->getAcceleNew() << std::endl;
//          std::cout << "velOld= " << cur_point->getVelOld() << std::endl;
//          std::cout << "velNew= " << cur_point->getVelNew() << std::endl;
//        }

        cur_point->setDispOld(cur_point->getDispNew());
        cur_point->setVelOld(cur_point->getVelNew());
        cur_point->setAccele(cur_point->getAcceleNew());
        cur_point->setPosNew(cur_point->getPosOld() + cur_point->getDisp());


        if (writing_condition)
        {
          std::ofstream myfile;
          myfile.open(fileName, std::ofstream::app);
          myfile << "{AcceleNew[" << num_iter + 1 << "]- Accele[" << num_iter << "]}*0.5= " << averageAccele << std::endl;
          myfile << "VelNew[" << num_iter + 1 << "]= " << cur_point->getVelNew() << std::endl;
          myfile << "PosNew[" << num_iter + 1 << "]= " << cur_point->getPosNew() << std::endl;
          myfile.close();
          ret = chdir(currentDirec.c_str());

        }  


      }

    }




    std::vector<StateMaterialPointP> statePointsArray = d_peri_box->getStatePoints();
    int size = statePointsArray.size();
////    int id = size/4 - 1;
    for (i = 0; i < size; i++)
    {  
      StateMaterialPointP cur_point = statePointsArray[i];
      if (cur_point->getID() == id) 
//       size/4 - 1)
      {
        std::cout << "ID= " << cur_point->getID() << " Disp= (" <<
                               cur_point->getDisp().x() << ", " <<
                               cur_point->getDisp().y() << ", " <<
                               cur_point->getDisp().z() << ")" << " Damage= " <<
                               cur_point->getLocalDamage()<< std::endl;
//        std::cout  << "DispNew= " << cur_point->getDispNew() << std::endl;
//      if ((cur_point->getLocalDamage() != 1))// && (cur_point->getLocalDamage() < 1))
//        std::cout << "ID= " << cur_point->getID() << " Damage= " <<
//                               cur_point->getLocalDamage() << std::endl;
      }  
    }


////   char buffer[2000];
////   char *str = getcwd(buffer, 2000);
////   std::string currentDirec = std::string(buffer);
    
//    if (num_iter % 5 == 1)
//    if (num_iter % 10 == 1)
//    {
         pointsFileState(number, currentDirec);

//      if (num_iter > 80 && num_iter <90)
//      {
//         printDetails(number, currentDirec);    
//      }
      number++;            
//    }

    num_iter++;  

  }

}




/*******
void
Peridynamics::damageModelingUsingCriticalStrain(const double& yield, const int& num_iter)
{

  StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
  int pointsSize = statePoints.size();
  int i = 0;
  bool flag = false;
//  bool cond = false;
  do
  {
    flag = false;

    for (i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
      if (!(cur_point->getLocalDamage() == 0))
      {
        cur_point->LagrangianStrainTensor();
      }
      
    }


    for (i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
      if (!(cur_point->getLocalDamage() == 0))
      {
        cur_point->checkBrokenBonds(yield, flag);
      }

      cur_point->calculateStateLocalDamage();
  
    }

  } while (flag);


  for (i = 0; i < pointsSize; i++)
  {
    StateMaterialPointP cur_point = statePoints[i];
    if (!(cur_point->getLocalDamage() == 0))
    {
      cur_point->calculateStateLocalDamage();
    }  

  }

}
*************/


void
Peridynamics::damageModelingUsingCriticalStrain(const double& yield, const int& num_iter)
{

  char buffer[2000];
  char *str = getcwd(buffer, 2000);
  std::string currentDirec = std::string(buffer);


  int ret;
//  bool writing_condition = ((num_iter > 20) && (num_iter < 45));
  bool writing_condition = (num_iter < 490);
  bool nodeID = true;
//  std::string fileName = "damage_modeling_detail_" + std::to_string(num_iter) + ".com";

  StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
  int pointsSize = statePoints.size();
  int i = 0;
  int k = 0;
  int l = 0;

  int id = 3909;
  int id2 = 4184;

  bool flag = false;
//  bool cond = false;
  do
  {
    flag = false;

    for (i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
      std::string fileName = "damage_modeling_detail_" + std::to_string(cur_point->getID()) + "_" + std::to_string(num_iter) + ".com";
      if (!(cur_point->getLocalDamage() == 0))
      {
        cur_point->LagrangianStrainTensor();
      }

      nodeID = (cur_point->getID() == id) || (cur_point->getID() == id2);
//      writing_condition = nodeID && writing_condition;
//      if (nodeID && writing_condition)
//      {
//        std::cout << "************************ writing_condition= " << writing_condition << std::endl;
//        std::cout << "************************ nodeID= " << nodeID << std::endl;
//        std::cout << "************************ pointID= " << cur_point->getID() << std::endl;
//      }

      
//      if ((cur_point->getID() == id) && writing_condition)
      if (nodeID && writing_condition)
//      if (cond)
      {
        std::ofstream myfile;
        myfile.open(fileName, std::ofstream::app);

        myfile << "Node ID= " << cur_point->getID() << std::endl << std::endl;
        myfile << "Deformation Gradient= F= " << std::endl << cur_point->getDeformationGradient();
        myfile << "rightCauchyGreenDeformation= F.Transpose*F= " << std::endl << cur_point->getRightCauchyGreenDeformation();
        myfile << "Lagrangian Strain Tensor C= (F.Transpose*F - I)*0.5= " << std::endl << cur_point->getLagrangianStrainTensor();
        myfile << "Yield= " << yield << std::endl;

        myfile << std::endl;
        myfile.close();
        ret = chdir(currentDirec.c_str());
      } 
    }


    for (k = 0; k < pointsSize; k++)
    {
      StateMaterialPointP cur_point = statePoints[k];
      std::string fileName = "damage_modeling_detail_" + std::to_string(cur_point->getID()) + "_" + std::to_string(num_iter) + ".com";
      nodeID = (cur_point->getID() == id) || (cur_point->getID() == id2);

      if (!(cur_point->getLocalDamage() == 0))
      {
        cur_point->checkBrokenBonds(yield, flag);
      }

//      if ((cur_point->getID() == id) && writing_condition)
      if (nodeID && writing_condition)
//      if (cond)
      {
        std::ofstream myfile;
        myfile.open(fileName, std::ofstream::app);

        StateBondPArray bondsArray = cur_point->getStateBonds();
        int size = bondsArray.size();
        int j = 0;
        for (j = 0; j < size; j++)
        {
          StateBondP bond = bondsArray[j];
          StateMaterialPointP centerPoint = bond->getCenterPoint();
          StateMaterialPointP secondPoint = bond->getSecondPoint();
          myfile << "********* (" << j + 1 << ")  Bond(" << bond->getCenterPoint()->getID() << ", " << bond->getSecondPoint()->getID() << ") *********  Broken? "
                 << bond->getBrokenBond() << "  ******" << std::endl;
          myfile << "Lagrangian Strain Tensor of second point " << bond->getSecondPoint()->getID() 
                 << "= " << std::endl << bond->getSecondPoint()->getLagrangianStrainTensor();
          myfile << "Average Strain Tensor= " << std::endl << bond->getAverageStrainTensor();
          myfile << "Equivalent Strain= " << bond->getEquivalentStrain() << std::endl;
          myfile << "Yield= " << yield << std::endl;
//          myfile << "Bond Froce Vector= " << bond->getForceVectorState() << std::endl; 
          myfile << "Bond Influence Function= " << bond->getInfluenceFunction() << std::endl;
          myfile << "Bond Center Point Reference Position= " << centerPoint->getPosOld() << std::endl;
          myfile << "Bond Second Point Reference Position= " << secondPoint->getPosOld() << std::endl;  
          myfile << "Bond Sai= " << bond->getSai() << std::endl;
//          myfile << "Bond Polynomial Basis Order 2= " << bond->getPolynomialBasisOrder2() << std::endl;
          myfile << "Bond Volume= " << bond->getVolume() << std::endl;    
          myfile << "Bond Y<Sai>= Deformation Vector= " << bond->getDeformationVector() << std::endl;

//          myfile << std::endl;        

        }

        myfile << "Point Damage Index= " << cur_point->getLocalDamage() << std::endl;
        myfile << std::endl;

        myfile.close();
        ret = chdir(currentDirec.c_str());
      }

//      cur_point->calculateStateLocalDamage();
  
    }

  } while (flag);


  for (l = 0; l < pointsSize; l++)
  {
    StateMaterialPointP cur_point = statePoints[l];
//    if (!(cur_point->getLocalDamage() == 0))
//    {
      cur_point->calculateStateLocalDamage();
//    }

    if (nodeID && writing_condition)
//      if (cond)
    {
      std::string fileName = "damage_modeling_detail_" + std::to_string(cur_point->getID()) + "_" + std::to_string(num_iter) + ".com";
      std::ofstream myfile;
      myfile.open(fileName, std::ofstream::app); 
      nodeID = (cur_point->getID() == id) || (cur_point->getID() == id2);
      myfile << "Point Damage Index= " << cur_point->getLocalDamage() << std::endl;
      myfile << std::endl;
      myfile.close();
      ret = chdir(currentDirec.c_str());
    }
 

  }

}






void
Peridynamics::damageModelingUsingCriticalStrain(const double& yield)
{

  StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
  int pointsSize = statePoints.size();
  int i = 0;
  bool flag = false;
  do
  {
    flag = false;

    for (i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
      if (!(cur_point->getLocalDamage() == 0))
      {
        cur_point->LagrangianStrainTensor();
      }

    }


    for (i = 0; i < pointsSize; i++)
    {
      StateMaterialPointP cur_point = statePoints[i];
      if (!(cur_point->getLocalDamage() == 0))
      {
        cur_point->checkBrokenBonds(yield, flag);
      }
//      cur_point->calculateStateLocalDamage();
  
    }

  } while (flag);


  for (i = 0; i < pointsSize; i++)
  {
    StateMaterialPointP cur_point = statePoints[i];
    if (!(cur_point->getLocalDamage() == 0))
    {
      cur_point->calculateStateLocalDamage();
    }  

  }

}







  void
  Peridynamics::pointsFileState(int& iter_num, const std::string& currentFolder)
  {
//    std::string pointFolderName = "pointsFolder";
//    int status2;
    int ret2;
//    if (iter_num == 0)
//    {

//      makeOutputFolder(status2, ret2, pointFolderName);
//    }
//    else
//    {
//      std::string direcLocation = currentFolder + "/" + pointFolderName;
//      ret2 = chdir(direcLocation.c_str());
//    }
    std::string fileName = "box_" + std::to_string(iter_num) + ".exdata";
    std::ofstream myfile(fileName);
    myfile << " Group name: box" // << "_" << std::to_string(iter_num) 
                                  << std::endl;
//    myfile << " #Fields=1" << std::endl;
    myfile << " #Fields=2" << std::endl;
    myfile << " 1) coordinates, coordinate, rectangular cartesian, #Components=3" << std::endl;
    myfile << "   x.  Value index= 1, #Derivatives= 0" << std::endl;
    myfile << "   y.  Value index= 2, #Derivatives= 0" << std::endl;
    myfile << "   z.  Value index= 3, #Derivatives= 0" << std::endl;
    myfile << " 2) localDamage, field, rectangular cartesian, #Components=1" << std::endl;
//    myfile << " 3) equivalentStrain, field, rectangular cartesian, #Components=1" << std::endl;
//    myfile << " 2) stressX, field, rectangular cartesian, #Components=1" << std::endl;
    myfile << "   localDamage. Value index= 4, #Derivatives= 0" << std::endl;
//    myfile << "   equivalentStrain. Value index= 5, #Derivatives= 0" << std::endl;
//    myfile << "   stressX. Value index= 4, #Derivatives= 0" << std::endl;


    std::vector<StateMaterialPointP> statePointsArray = d_peri_box->getStatePoints();
    int size = statePointsArray.size();
    for (int i = 0; i < size; i++)
    {  
      StateMaterialPointP cur_point = statePointsArray[i];
      myfile << " Node:         " << cur_point->getID()  << std::endl;
      myfile << "  " << cur_point->getPosNew().x() << "  "
             << cur_point->getPosNew().y() << "  "
             << cur_point->getPosNew().z() << "  "
             << cur_point->getLocalDamage() << "  "
//             << cur_point->getPeriStress().x()  
             << std::endl;
    }
    
    myfile.close();
    ret2 = chdir(currentFolder.c_str());

  }




void
Peridynamics::printDetails(int& iter_num, const std::string& currentFolder)
{
  int id = 3762;
  std::vector<StateMaterialPointP> statePointsArray = d_peri_box->getStatePoints();
  int size = statePointsArray.size();
  for (int i = 0; i < size; i++)
  {  
    StateMaterialPointP cur_point = statePointsArray[i];
    if (cur_point->getID() == id)
    {
      if (getTimeIntegration() == Implicit)
      {
        Matrix3D BroydenCurrent(0.0);
        Matrix3D BroydenNew(0.0);
        Vector3D residualVectorCurrent(0.0);
        Vector3D residualVectorNew(0.0);
        cur_point->getInformation(BroydenCurrent,
                                  BroydenNew,
                                  residualVectorCurrent,
                                  residualVectorNew);


      }

    }

  }

}




void
Peridynamics::testDerivative()
{

  bool secondOrder = true;
//  bool secondOrder = false;

  StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
  int pointsSize = statePoints.size();
  int i = 0;
  for (i = 0; i < pointsSize; i++)
  {
    StateMaterialPointP cur_point = statePoints[i];
    if (cur_point->getID() == 3920)
    {    
      Vector3D pos = cur_point->newPositionVector(); 
      std::cout << "pos= " << pos << std::endl;
      cur_point->derivative(secondOrder);
    }
  }

}




/*void
Peridynamics::testDivergance()
{

//  bool secondOrder = true;
  bool secondOrder = false;

  bool diverganceMethod = true;
//  bool diverganceMethod = false;


  StateMaterialPointPArray statePoints = d_peri_box->getStatePoints();
  int pointsSize = statePoints.size();
  int i = 0;
  for (i = 0; i < pointsSize; i++)
  {
    StateMaterialPointP cur_point = statePoints[i];
    if (cur_point->getID() == 1440)
    {    
      Vector3D pos = cur_point->newPositionVector(); 
      std::cout << "pos= " << pos << std::endl;
      cur_point->divergance(diverganceMethod, secondOrder);
    }
  }

}*/





double
Peridynamics::ramp(const int& num_iter, const int& num_initial)
{
  
  if (num_initial == 1)
  {
    return 1.0;
  }
  else
  {
    double denominator = (double) (num_initial - 1);    
    double coeff =  1/denominator;
//    if (getID() == 4224)
//    {
//      std::cout << "coeff= " << coeff << std::endl;
//    } 
    return coeff*num_iter;
  }


}


double
Peridynamics::triangle(const int& num_iter, const int& num_initial)
{
  double result = 0.0;
  if (num_initial == 1)
  {
    result = 1.0;
  }
  else
  {
    double denominator = (double) (num_initial - 1);    
    double coeff =  (2.0)/denominator;
//    if (getID() == 4224)
//    {
//      std::cout << "coeff= " << coeff << std::endl;
//    } 
    double half = denominator/(2.0);
    double num = (double) (num_iter);
    if ((0.0 <= num) && (num <= half))
    {
      result = coeff*num;
    }
    else if ((half <= num) && (num <= denominator))
    {
      result = 2 - coeff*num;
    }
    else
    {
      result = 0.0;
    }
  }

  return result;
  
}






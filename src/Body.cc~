#include <Body.h>
#include <PeriMaterialPointP.h>
#include <PeriMaterialPoint.h>
#include <StateMaterialPointP.h>
#include <StateMaterialPoint.h>
#include <StateBond.h>
#include <StateBondP.h>
#include <PeriBond.h>
#include <PeriBondP.h>
#include <PeriBondPArray.h>
#include <PeriPMBMaterial.h>
#include <GeometryMath/Vector3D.h>
#include <GeometryMath/Matrix3D.h>

#include <fstream>

#include <cmath>
#include <memory>

using namespace FiniteElement;

  Body::Body() 
  : d_density(0.0), d_poission_ratio(0.0), d_young_modulus(0.0)
  {
    d_elements_array.clear();
    d_nodes_array.clear();
    d_global_stiffness.clear();
    d_global_mass.clear();
    d_global_lumped_mass.clear();
    d_external_force.clear();
    d_glob_old_disp.clear();
    d_glob_disp.clear();
    d_glob_new_disp.clear();
    d_elements_array.reserve(1500);
    d_nodes_array.reserve(2000);
    d_global_stiffness.reserve(6000);
    d_global_mass.reserve(6000);
    d_global_lumped_mass.reserve(6000);
    d_external_force.reserve(6000);
    d_glob_old_disp.reserve(6000);
    d_glob_disp.reserve(6000);
    d_glob_new_disp.reserve(6000);
  }

  Body::~Body()
  {
  }


  void
  Body::clearStiff()
  {
    d_global_stiffness.clear();
  }


  void
  Body::clearMass()
  {
    d_global_mass.clear();
  }

  void
  Body::clearLumpedMass()
  {
    d_global_lumped_mass.clear();
  }

  void
  Body::clearExternalForce()
  {
    d_external_force.clear();
  }

  void
  Body::clearOldGlobDisp()
  {
    d_glob_old_disp.clear();
  }

  void
  Body::clearGlobDisp()
  {
    d_glob_old_disp.clear();
  }

  void
  Body::clearNewGlobDisp()
  {
    d_glob_new_disp.clear();
  }

  void
  Body::clearNodesArray()
  {
    d_nodes_array.clear();
  }

  void
  Body::clearElementsArray()
  {
    d_elements_array.clear();
  }


  void
  Body::initialMatrix(std::vector <std::vector <double> >& Matrix, int size)
  {
//    for (int i = 0; i < size; i++)
//      for (int j = 0; j < size; j++)
//        Matrix[i][j] = 0.0;
    std::vector <double> row;
    for (int j = 0; j < size; j++)
    {
      for (int i = 0; i < size; i++)
        row.push_back(0.0);
      Matrix.push_back(row);
      row.clear();
    }
//    for (int ii = 0; ii < size; ii++)
//      for (int jj = 0; jj < size; jj++)
//      std::cout << Matrix[ii][jj] << std::endl;
    
  }

  void
  Body::initialVector(std::vector <double>& Vector, int size)
  {
    for (int i = 0; i < size; i++)
    {
      Vector.push_back(0.0);
    }
  }

  void
  Body::setProperties(double density, double young, double poission)
  {
    d_elements_array.clear();
    d_nodes_array.clear();
    d_global_stiffness.clear();
    d_global_mass.clear();
    d_glob_old_disp.clear();
    d_glob_disp.clear();
    d_glob_new_disp.clear();
 
    setDensity(density);
    setPoission(poission);
    setYoung(young);
  }


  void
  Body::setProperties(double density, double young, double poission, double yield)
  {
    d_elements_array.clear();
    d_nodes_array.clear();
    d_global_stiffness.clear();
    d_global_mass.clear();
    d_glob_old_disp.clear();
    d_glob_disp.clear();
    d_glob_new_disp.clear();
 
    setDensity(density);
    setPoission(poission);
    setYoung(young);
    setYieldStrength(yield);
  }



  bool
  Body::isNodeInElement(NodeP node, SymmetricMaterialTrilinearElementSP element)
  {
    bool belong = false;
    TrilinearVolumeElementSP elem = element->triElement();
    std::array <NodeP, 8> nodes = elem->getArray();
    for (int i = 0; i < 8; i++)
    {
      if (node->getID() == nodes[i]->getID())
      {
        belong = true;
      break;
      } 
    }
    return belong;
  }


  int
  Body::globalIndex(SymmetricMaterialTrilinearElementSP element, int i)
  {
    int index = i/3;
    std::array <NodeP, 8> nodes = element->triElement()->getArray();
    int node_id = nodes[index]->getID();
//      if (element->id() == 1)
//       std::cout << "nodeID= " << node_id << ", ";
    return 3*node_id-3+(i%3);
//    return 3*(element->triElement()->getArray())[i/3]->getID()-3+i%3;
  }



  void
  Body::makeGlobalMatrices()
  {
    std::vector<std::vector<double> > globStiff;
    std::vector<std::vector<double> > globMass;
    std::vector<std::vector<double> > globLumpedMass;

    std::vector<double> globExternalForce;
    std::vector<double> globOldDisp;
    std::vector<double> globDisp;
    std::vector<double> globNewDisp;

    int size = 3*getNodes().size();
//    std::cout << size << std::endl;
    initialMatrix(globStiff, size);
    initialMatrix(globMass, size);
    initialMatrix(globLumpedMass, size);
    initialVector(globExternalForce, size);
    initialVector(globOldDisp, size);
    initialVector(globDisp, size);
    initialVector(globNewDisp, size);
//    std::cout << "size= " << size << std::endl;
//    std::cout << "globMass size= " << globMass.size() << std::endl;
//    for (int ii = 0; ii < globMass.size(); ii++)
//    {
//      for (int jj = 0; jj < globMass.size(); jj++)
//      {
//        std::cout << "globMass(" << ii << ", " << jj << ")= " << 
//                      globMass[ii][jj] << "   ";
//      }
//    }    

    int globalJ = 0;
    int globalK = 0;
//    std::cout << d_elements_array.size() << std::endl;   
//    for (auto element_iter = d_elements_array.begin(); element_iter != d_elements_array.end(); element_iter++)
    for (int element_iter = 0; element_iter < d_elements_array.size(); element_iter++)
    {
//      SymmetricMaterialTrilinearElementSP cur_element = *element_iter;
      SymmetricMaterialTrilinearElementSP cur_element = d_elements_array[element_iter];
      cur_element->stiffnessMatrix(6);
      std::cout << cur_element->id() << " StiffnessMatrix of elements are created." << std::endl;
      cur_element->massMatrix(6);
      cur_element->surfaceForceVector(6);
//      std::cout << cur_element->id() << " SurfaceForceVector of elements are created." << std::endl;
      cur_element->bodyForceVector(6);
//      cur_element->triElement()->oldDisplacementArray();
//      cur_element->triElement()->displacementArray();
//      cur_element->triElement()->newDisplacementArray();
      cur_element->createDispVectors();
      for (int j = 0; j < 24 ; j++)
      { 
        globalJ = globalIndex(cur_element, j);
        globExternalForce[globalJ] += cur_element->surfaceVec()(j) + cur_element->bodyVec()(j);
        globOldDisp[globalJ] = cur_element->oldDispVec()(j);
        globDisp[globalJ] = cur_element->dispVec()(j);
        globNewDisp[globalJ] = cur_element->newDispVec()(j);
//        std::cout << "surface force(" << j << ")=" << cur_element->surfaceVec()(j) << "  ";
        for (int k = 0; k < 24; k++)
          { 
            globalK = globalIndex(cur_element, k);
            
//            if (element_iter == 0) 
//            std::cout << "j= " << j << ", k= " << k
//                      << ", globalJ= " << globalJ
//                      << ", globalK= " << globalK << std::endl; 
            globStiff[globalJ][globalK] += cur_element->stiffMat()(j, k);
            globMass[globalJ][globalK] += cur_element->massMat()(j, k);
            globLumpedMass[globalJ][globalK] += cur_element->lumpedMassMat()(j, k);
//            double test = cur_element->stiffMat()(j, k);
          }

      }
    }
   setGlobalStiffness(globStiff);
   setGlobalMass(globMass);
   setGlobalLumpedMass(globLumpedMass);
   setExternalForce(globExternalForce);
   setOldGlobDisp(globOldDisp);
   setGlobDisp(globDisp);
   setNewGlobDisp(globNewDisp);
   std::cout << getExternalForce().size() << std::endl;

//   int sizee = getExternalForce().size();
//   for (int ll = 0; ll < sizee; ll++)
//   {
//      if ((ll%3 == 0) /*&& (ll != 0)*/) std::cout << std::endl << ll/3+1 << "   ";
//      std::cout << "(" << ll << ")= " << getExternalForce()[ll] << "  ";
//   }
//   std::cout << std::endl;


    

}


//**************************************************** Peridynamics ***********************************************

void
Body::setPeriPoints()
{
  std::vector<MaterialPointP> pointsArray = getMaterialPoints();
  std::vector<PeriMaterialPointP> periPointsArray;    
  int size = pointsArray.size();
//  std::cout << size << std::endl;
  for (int i = 0; i < size; i++)
  {
    MaterialPointP cur_point = pointsArray[i];
    PeriMaterialPointP peri_point(new PeriMaterialPoint(cur_point));
    periPointsArray.emplace_back(peri_point);
    
  }
  setMaterialPeriPoints(periPointsArray);

//  for (int i = 0; i < size; i++)
//  {
//    std::cout << getMaterialPeriPoints()[i]->getID() << "  ";
//  }
}




void
Body::createPeriBondsOfEachPeriPoints(const double& horizon)
{
  std::vector<PeriMaterialPointP> periPointsArray = getMaterialPeriPoints();
  int size = periPointsArray.size();
//  std::cout << size << std::endl;
  for (int i = 0; i < size; i++)
  {
    PeriMaterialPointP cur_peri_point = periPointsArray[i];
//    std::cout << cur_peri_point->getID() << std::endl;
    PeriBondPArray bondsArray;
//    PeriBondP bond;
//    PeriBond bond;
//    std::shared_ptr<PeriBond> bond;
    bondsArray.clear();
    for (int j = 0; j < size; j++)
    {
//      std::cout << j << std::endl;
      PeriMaterialPointP neighb_point = periPointsArray[j]; 
      Vector3D sai = neighb_point->getPosOld() - cur_peri_point->getPosOld();
//      std::cout << sai.x() << " " << sai.y() << " " << sai.z() << " " << std::endl;
      if ((j != i) && (sai.length() <= horizon))
      {
        PeriBondP bond(new PeriBond(cur_peri_point, neighb_point, false));
//        bond->setBrokenBond(false);
//        bond.setBrokenBond(false);

//        bond->setCenterPoint(cur_peri_point);
//        bond.setCenterPoint(cur_peri_point);
//        bond->setSecondPoint(neighb_point);
//        bond.setSecondPoint(neighb_point);
//        std::cout << "MMMM" << std::endl;  
        bondsArray.emplace_back(bond);
//        continue;
      }
    }
//    std::cout << "PointID= " << cur_peri_point->getID() << " BondSize= " << bondsArray.size() << std::endl;
    cur_peri_point->setHorizonBonds(bondsArray);    
  }

  setMaterialPeriPoints(periPointsArray);
//  std::cout << "OK" << std::endl;
}  





void
Body::calculateAccelerationOfEachPeriPoints()
{
  std::vector<PeriMaterialPointP> periPointsArray = getMaterialPeriPoints();
  int size = periPointsArray.size();
  for (int i = 0; i < size; i++)
  {
    PeriMaterialPointP cur_peri_point = periPointsArray[i];
    Vector3D internal_force(0.0, 0.0, 0.0);
    PeriBondPArray bondsArray = cur_peri_point->getHorizonBonds();
    int boundSize = bondsArray.size();
    for (int j = 0; j < boundSize; j++)
    { 
      PeriBondP bond = bondsArray[j];
      PeriMaterialPointP second_point = bond->getSecondPoint();
      double volume = second_point->getInitialVolume();
//      internal_force = internal_force + bond->getPairwiseForce()*volume*bond->getSurfaceCorrectionFactor();
      internal_force = internal_force + bond->getPairwiseForce()*bond->getVolume();

    }
    
//    Vector3D zero(0.0, 0.0, 0.0);
//    if (!(cur_peri_point->getBodyForce() == zero))
//    if (cur_peri_point->getID() == size/3 + 1)    
//      std::cout << "ID= " << cur_peri_point->getID() << " External= " << cur_peri_point->getBodyForce()
//                          << " Internal= " << internal_force << std::endl << std::endl;    
//      std::cout << "internalForce= " << internal_force << std::endl << std::endl;

    Vector3D acceleration = (internal_force + cur_peri_point->getBodyForce())/(d_density);
    periPointsArray[i]->setAccele(acceleration);
//    if (cur_peri_point->getID() == 1)
//      std::cout << "ID= " << cur_peri_point->getID() << " Acceleration= " << cur_peri_point->getAccele() << std::endl; 
  }

  setMaterialPeriPoints(periPointsArray);

}






void
Body::calculateDisplacementAndVelocityOfPointsUsingVelocityVerlet(const double& del_t)
{

//  std::cout << "+++++++++++++++++++++++++Hello Velocity Verlet++++++++++++++++++++++++++" << std::endl;
  std::vector<PeriMaterialPointP> periPointsArray = getMaterialPeriPoints();
  int size = periPointsArray.size();
  for (int i = 0; i < size; i++)
  {
    PeriMaterialPointP cur_point = periPointsArray[i];
//    if (cur_point->getID() == size/4 + 1)
//    {
//      std::cout << "ID= " << cur_point->getID() << " velOld(n)= " << cur_point->getVelOld() << std::endl;
//      std::cout << "acceleration(n)= " << cur_point->getAccele() << " delT= " << del_t << std::endl;
//    }
    cur_point->setVelMid(cur_point->getVelOld() + cur_point->getAccele()*(0.5*del_t)); // v(n+1/2) = v(n) + a(n)*del_t/2
    if (cur_point->getVelocityBoundaryFlag() == true)
    {
      cur_point->setVelMid(cur_point->getVelOld());
    }
//    if (cur_point->getID() == size/4 + 1)
//    {
//      std::cout << "velMid(n+0.5)= " << cur_point->getVelMid() << std::endl;
//      std::cout << "disp(n)= " << cur_point->getDisp() << std::endl;
//    }
    cur_point->setDispNew(cur_point->getDisp() + cur_point->getVelMid()*(del_t)); // u(n+1) = u(n) + v(n+1/2)*del_t
//    if (cur_point->getID() == size/4 + 1)
//    {
//      std::cout << "dispNew(n+1)= " << cur_point->getDispNew() << std::endl;
//    }
    cur_point->setDisp(cur_point->getDispNew());
//    FiniteElement::PeriPMBMaterial::calculatePairwiseForceOfEachBond();



  }
  setMaterialPeriPoints(periPointsArray);
  
}




//******************************************** StateBased Peridynamics **************************************


void
Body::setStatePoints()
{
  std::vector<PeriMaterialPointP> pointsArray = getMaterialPeriPoints();
  std::vector<StateMaterialPointP> statePointsArray;    
  int size = pointsArray.size();
//  std::cout << size << std::endl;
  for (int i = 0; i < size; i++)
  {
    PeriMaterialPointP cur_peri_point = pointsArray[i];
    StateMaterialPointP peri_point(new StateMaterialPoint(cur_peri_point));
    statePointsArray.emplace_back(peri_point);
    
  }
  setStatePoints(statePointsArray);
//  std::cout << "Number of state points= " << getStatePoints().size() << std::endl;

//  for (int i = 0; i < size; i++)
//  {
//    std::cout << getStatePoints()[i]->getID() << "  ";
//  }

}



void
Body::setStatePoints(const Matrix3D& initialVelocityGradient)
{
  std::vector<PeriMaterialPointP> pointsArray = getMaterialPeriPoints();
  std::vector<StateMaterialPointP> statePointsArray;    
  int size = pointsArray.size();
//  std::cout << size << std::endl;
  for (int i = 0; i < size; i++)
  {
    PeriMaterialPointP cur_peri_point = pointsArray[i];
    StateMaterialPointP peri_point(new StateMaterialPoint(cur_peri_point, initialVelocityGradient));
    statePointsArray.emplace_back(peri_point);
    
  }
  setStatePoints(statePointsArray);
//  std::cout << "Number of state points= " << getStatePoints().size() << std::endl;

//  for (int i = 0; i < size; i++)
//  {
//    std::cout << getStatePoints()[i]->getID() << "  ";
//  }

}



/***************************************



void
Body::createStateBondsOfEachStatePoints(const double& horizon)
{
 
  std::vector<StateMaterialPointP> statePointsArray = getStatePoints();
  int size = statePointsArray.size();
  std::cout << size << std::endl;
  for (int i = 0; i < size; i++)
  {
    StateMaterialPointP cur_state_point = statePointsArray[i];

//    std::cout << cur_state_point->getID() << "  " ; //<< std::endl;
    StateBondPArray bondsArray;
//    PeriBondP bond;
//    PeriBond bond;
//    std::shared_ptr<PeriBond> bond;
    bondsArray.clear();
    for (int j = 0; j < size; j++)
    {
//      std::cout << j << std::endl;
      StateMaterialPointP neighb_point = statePointsArray[j]; 
      Vector3D sai = neighb_point->getPosOld() - cur_state_point->getPosOld();
//      std::cout << sai.x() << " " << sai.y() << " " << sai.z() << " " << std::endl;
      if ((j != i) && (sai.length() <= horizon))
      {
//        StateBondP bond(new StateBond(cur_state_point, neighb_point, false));



         StateBondP bond(new StateBond(cur_state_point, neighb_point, horizon, false));
        
//         bond->setBrokenBond(false);
//         bond.setBrokenBond(false);
//         bond->setCenterPoint(cur_peri_point);
//         bond.setCenterPoint(cur_peri_point);
//         bond->setSecondPoint(neighb_point);
//         bond.setSecondPoint(neighb_point);
//         std::cout << "MMMM" << std::endl;  
         bondsArray.emplace_back(bond);
       
//          continue;
       }
    
    }
//    std::cout << "PointID= " << cur_state_point->getID() << " BondSize= " << bondsArray.size() << std::endl;
    cur_state_point->setStateBonds(bondsArray);    
  }

//  setMaterialPeriPoints(periPointsArray);
//  std::cout << "OK" << std::endl;

}

***************************/



void
Body::createStateBondsOfEachStatePoints(const double& horizon)
{


  int id = 3909;
  int id2 = 4184;

  char buffer[2000];
  char *str = getcwd(buffer, 2000);
  std::string currentDirec = std::string(buffer);
  int ret;

  int numBond = 0;

  std::vector<StateMaterialPointP> statePointsArray = getStatePoints();
  int size = statePointsArray.size();
  std::cout << size << std::endl;
  for (int i = 0; i < size; i++)
  {
    StateMaterialPointP cur_state_point = statePointsArray[i];

//    std::cout << cur_state_point->getID() << "  " ; //<< std::endl;
    StateBondPArray bondsArray;
//    PeriBondP bond;
//    PeriBond bond;
//    std::shared_ptr<PeriBond> bond;
    bondsArray.clear();
    for (int j = 0; j < size; j++)
    {
//      std::cout << j << std::endl;
      StateMaterialPointP neighb_point = statePointsArray[j]; 
      Vector3D sai = neighb_point->getPosOld() - cur_state_point->getPosOld();

      if ((cur_state_point->getID() == id) || (cur_state_point->getID() == id2))
      {
        std::string fileName = "create_bonds_" + std::to_string(cur_state_point->getID()) + ".com";
        std::ofstream myfile;
        myfile.open(fileName, std::ofstream::app);
        myfile << "Point " << cur_state_point->getID() << ": " << std::endl;
        myfile << "(" << cur_state_point->getID() << ", " << neighb_point->getID() << "):" << std::endl;
        myfile << "Center Point Reference Position= " << cur_state_point->getPosOld() << std::endl;
        myfile << "Second Point Reference Position= " << neighb_point->getPosOld() << std::endl;
        myfile << "Sai= SecondPos - CenterPos= " << sai << std::endl;
        myfile << "SaiLength= " << sai.length() << std::endl;
        myfile << "Horizon= " << horizon << std::endl;
        myfile.close();
        ret = chdir(currentDirec.c_str());
      }



//      std::cout << sai.x() << " " << sai.y() << " " << sai.z() << " " << std::endl;
      if ((j != i) && (sai.length() <= horizon))
      {
//        StateBondP bond(new StateBond(cur_state_point, neighb_point, false));



         StateBondP bond(new StateBond(cur_state_point, neighb_point, horizon, false));
        
//         bond->setBrokenBond(false);
//         bond.setBrokenBond(false);
//         bond->setCenterPoint(cur_peri_point);
//         bond.setCenterPoint(cur_peri_point);
//         bond->setSecondPoint(neighb_point);
//         bond.setSecondPoint(neighb_point);
//         std::cout << "MMMM" << std::endl;  
         bondsArray.emplace_back(bond);



         if (cur_state_point->getID() == id)
         {
           std::string fileName = "create_bonds_" + std::to_string(cur_state_point->getID()) + ".com";
           std::ofstream myfile;
           myfile.open(fileName, std::ofstream::app);
           numBond ++;          
           myfile << "Bond number " << numBond << " =(" << bond->getCenterPoint()->getID() << ", " << bond->getSecondPoint()->getID() << ") is created." << std::endl;
           myfile << std::endl;
           myfile.close();
           ret = chdir(currentDirec.c_str());
         }          
//          continue;
       }
    
       else
       {
         if (cur_state_point->getID() == id)
         {
           std::string fileName = "create_bonds_" + std::to_string(cur_state_point->getID()) + ".com";
           std::ofstream myfile;
           myfile.open(fileName, std::ofstream::app);
           myfile << "No bond is created." << std::endl;
           myfile << std::endl;
           myfile.close();
           ret = chdir(currentDirec.c_str()); 
          } 

        }

    }
//    std::cout << "PointID= " << cur_state_point->getID() << " BondSize= " << bondsArray.size() << std::endl;
    cur_state_point->setStateBonds(bondsArray);    
  }

//  setMaterialPeriPoints(periPointsArray);
//  std::cout << "OK" << std::endl;
}  


void
Body::createStateBondsOfEachStatePoints(const double& delX, const double& horizon)
{
  std::vector<StateMaterialPointP> statePointsArray = getStatePoints();
  int size = statePointsArray.size();
  std::cout << size << std::endl;
  for (int i = 0; i < size; i++)
  {
    StateMaterialPointP cur_state_point = statePointsArray[i];
//    std::cout << cur_state_point->getID() << "  " ; //<< std::endl;
    StateBondPArray bondsArray;
//    PeriBondP bond;
//    PeriBond bond;
//    std::shared_ptr<PeriBond> bond;
    bondsArray.clear();
    for (int j = 0; j < size; j++)
    {
//      std::cout << j << std::endl;
      StateMaterialPointP neighb_point = statePointsArray[j]; 
      Vector3D sai = neighb_point->getPosOld() - cur_state_point->getPosOld();
//      std::cout << sai.x() << " " << sai.y() << " " << sai.z() << " " << std::endl;
      if ((j != i) && (sai.length() <= horizon))
      {

//        StateBondP bond(new StateBond(cur_state_point, neighb_point, false));
//        StateBondP bond(new StateBond(cur_state_point, neighb_point, horizon, false));
        StateBondP bond(new StateBond(cur_state_point, neighb_point, delX, false));


//        bond->setBrokenBond(false);
//        bond.setBrokenBond(false);

//        bond->setCenterPoint(cur_peri_point);
//        bond.setCenterPoint(cur_peri_point);
//        bond->setSecondPoint(neighb_point);
//        bond.setSecondPoint(neighb_point);
//        std::cout << "MMMM" << std::endl;  
        bondsArray.emplace_back(bond);
//        continue;
      }
    }
//    std::cout << "PointID= " << cur_state_point->getID() << " BondSize= " << bondsArray.size() << std::endl;
    cur_state_point->setStateBonds(bondsArray);    
  }

//  setMaterialPeriPoints(periPointsArray);
//  std::cout << "OK" << std::endl;
}  




void
Body::accelerationOfEachStatePoints(const double& delT)
{
  Vector3D zero(0.0, 0.0, 0.0);
  std::vector<StateMaterialPointP> statePointsArray = getStatePoints();
  int size = statePointsArray.size();
  int i = 0;
  for (i = 0; i < size; i++)
  {
    StateMaterialPointP cur_state_point = statePointsArray[i];
    Vector3D internal_force(0.0, 0.0, 0.0);
    StateBondPArray bondsArray = cur_state_point->getStateBonds();
    int boundSize = bondsArray.size();
    int j = 0;
    for (j = 0; j < boundSize; j++)
    { 
      StateBondP bond = bondsArray[j];
//      if (bond->getBrokenBond() == false)
//      {
//        StateMaterialPointP second_point = bond->getSecondPoint();
//        double volume = second_point->getInitialVolume();
//        internal_force = internal_force + bond->getPairwiseForce()*volume*bond->getSurfaceCorrectionFactor();
//        StateBondP twin = twinBond(bond);
        internal_force = internal_force + (bond->getFirstPiolaDivergenceIntegrand())
                                           *bond->getVolume();
//        bool condition1 = ((bond->getCenterPoint()->getID() == 700) &&      
//                           (bond->getSecondPoint()->getID() == 613));

//        if (condition1)
//        {
//          std::cout << "Bond(" << bond->getCenterPoint()->getID() << "," 
//                               << bond->getSecondPoint()->getID() << ")" << std::endl;
//          std::cout << "ForceVectorState= " << bond->getForceVectorState() << std::endl;

//          std::cout << "TwinBond(" << twin->getCenterPoint()->getID() << "," 
//                               << twin->getSecondPoint()->getID() << ")" << std::endl;
//          std::cout << "TwinForceVectorState= " << twin->getForceVectorState() << std::endl;

//          std::cout << "Force exerted from point " << bond->getSecondPoint()->getID() 
//                    << " to the point " << bond->getCenterPoint()->getID() << " is= " << std::endl
//                    << bond->getForceVectorState() - twin->getForceVectorState() << std::endl;
//          std::cout << "Current deformation vector of the bond= " << bond->getDeformationVector() << std::endl;
//          std::cout << "Sai of the bond= " << bond->getSai() << std::endl;


//        }          
//      }

    }
    
//    Vector3D zero(0.0, 0.0, 0.0);
//    if (!(cur_peri_point->getBodyForce() == zero))
//    if (cur_peri_point->getID() == size/3 + 1)
//    if (cur_state_point->getID() == 700)    
//      std::cout << "ID= " << cur_peri_point->getID() << " External= " << cur_peri_point->getBodyForce()
//                          << " Internal= " << internal_force << std::endl << std::endl;
//      std::cout << "ID= " << cur_state_point->getID() << " External= " << cur_state_point->getBodyForce()
//                          << " Internal= " << internal_force << std::endl << std::endl;    
//      std::cout << "internalForce= " << internal_force << std::endl << std::endl;

    Vector3D acceleration = (internal_force + cur_state_point->getBodyForce())/(d_density);
    if (cur_state_point->getVelocityBoundaryFlag() == true)
    {
      statePointsArray[i]->setAccele(zero);
    }
    else
    {
      statePointsArray[i]->setAccele(acceleration);
    }


    cur_state_point->setVelNew(cur_state_point->getVelNew() + (cur_state_point->getAccele())*delT);
    cur_state_point->setDisp(cur_state_point->getDisp() + (cur_state_point->getVelNew())*delT);
//    cur_state_point->setPosNew(cur_state_point->getPosNew() + cur_state_point->getDisp());
    cur_state_point->setPosNew(cur_state_point->getPosOld() + cur_state_point->getDisp());

//    if (cur_peri_point->getID() == 1)
//    if (cur_state_point->getID() == 700)
//    {
//      std::cout << "ID= " << cur_state_point->getID() << " Acceleration= " << cur_state_point->getAccele() << std::endl;
//      std::cout << "VelNew= " << cur_state_point->getVelNew() << std::endl;
//      std::cout << "Disp= " << cur_state_point->getDisp() << std::endl;
//      std::cout << "PosNew= " << cur_state_point->getPosNew() << std::endl;
//      std::cout << "PosOld= " << cur_state_point->getPosOld() << std::endl;
//    }
     
  }

  setStatePoints(statePointsArray);

}



void
Body::calculateAccelerationOfEachStatePoints(const double& delT)
{
  Vector3D zero(0.0, 0.0, 0.0);
  std::vector<StateMaterialPointP> statePointsArray = getStatePoints();
  int size = statePointsArray.size();
  int i = 0;
  for (i = 0; i < size; i++)
  {
    StateMaterialPointP cur_state_point = statePointsArray[i];
    Vector3D internal_force(0.0, 0.0, 0.0);
    StateBondPArray bondsArray = cur_state_point->getStateBonds();
    int boundSize = bondsArray.size();
    int j = 0;
    for (j = 0; j < boundSize; j++)
    { 
      StateBondP bond = bondsArray[j];
//      if (bond->getBrokenBond() == false)
//      {
        StateMaterialPointP second_point = bond->getSecondPoint();
//        double volume = second_point->getInitialVolume();
//        internal_force = internal_force + bond->getPairwiseForce()*volume*bond->getSurfaceCorrectionFactor();
        StateBondP twin = twinBond(bond);
        internal_force = internal_force + (bond->getForceVectorState() - twin->getForceVectorState())
                                           *bond->getVolume();
//        bool condition1 = ((bond->getCenterPoint()->getID() == 700) &&      
//                           (bond->getSecondPoint()->getID() == 613));

//        if (condition1)
//        {
//          std::cout << "Bond(" << bond->getCenterPoint()->getID() << "," 
//                               << bond->getSecondPoint()->getID() << ")" << std::endl;
//          std::cout << "ForceVectorState= " << bond->getForceVectorState() << std::endl;

//          std::cout << "TwinBond(" << twin->getCenterPoint()->getID() << "," 
//                               << twin->getSecondPoint()->getID() << ")" << std::endl;
//          std::cout << "TwinForceVectorState= " << twin->getForceVectorState() << std::endl;

//          std::cout << "Force exerted from point " << bond->getSecondPoint()->getID() 
//                    << " to the point " << bond->getCenterPoint()->getID() << " is= " << std::endl
//                    << bond->getForceVectorState() - twin->getForceVectorState() << std::endl;
//          std::cout << "Current deformation vector of the bond= " << bond->getDeformationVector() << std::endl;
//          std::cout << "Sai of the bond= " << bond->getSai() << std::endl;


//        }          
//      }

    }
    
//    Vector3D zero(0.0, 0.0, 0.0);
//    if (!(cur_peri_point->getBodyForce() == zero))
//    if (cur_peri_point->getID() == size/3 + 1)
//    if (cur_state_point->getID() == 700)    
//      std::cout << "ID= " << cur_peri_point->getID() << " External= " << cur_peri_point->getBodyForce()
//                          << " Internal= " << internal_force << std::endl << std::endl;
//      std::cout << "ID= " << cur_state_point->getID() << " External= " << cur_state_point->getBodyForce()
//                          << " Internal= " << internal_force << std::endl << std::endl;    
//      std::cout << "internalForce= " << internal_force << std::endl << std::endl;

    Vector3D acceleration = (internal_force + cur_state_point->getBodyForce())/(d_density);
    if (cur_state_point->getVelocityBoundaryFlag() == true)
    {
      statePointsArray[i]->setAccele(zero);
    }
    else
    {
      statePointsArray[i]->setAccele(acceleration);
    }


    cur_state_point->setVelNew(cur_state_point->getVelNew() + (cur_state_point->getAccele())*delT);
    cur_state_point->setDisp(cur_state_point->getDisp() + (cur_state_point->getVelNew())*delT);
//    cur_state_point->setPosNew(cur_state_point->getPosNew() + cur_state_point->getDisp());
    cur_state_point->setPosNew(cur_state_point->getPosOld() + cur_state_point->getDisp());

//    if (cur_peri_point->getID() == 1)
//    if (cur_state_point->getID() == 700)
//    {
//      std::cout << "ID= " << cur_state_point->getID() << " Acceleration= " << cur_state_point->getAccele() << std::endl;
//      std::cout << "VelNew= " << cur_state_point->getVelNew() << std::endl;
//      std::cout << "Disp= " << cur_state_point->getDisp() << std::endl;
//      std::cout << "PosNew= " << cur_state_point->getPosNew() << std::endl;
//      std::cout << "PosOld= " << cur_state_point->getPosOld() << std::endl;
//    }
     
  }

  setStatePoints(statePointsArray);

}




void
Body::calculateAccelerationOfEachStatePoints(const double& delT, const double& springConstant)
{
  double c_hg = 0.015;
  Vector3D zero(0.0, 0.0, 0.0);
  std::vector<StateMaterialPointP> statePointsArray = getStatePoints();
  int size = statePointsArray.size();
  int i = 0;
  for (i = 0; i < size; i++)
  {
    StateMaterialPointP cur_state_point = statePointsArray[i];
    Matrix3D deformation = cur_state_point->getDeformationGradient();
    Vector3D internal_force(0.0, 0.0, 0.0);
    StateBondPArray bondsArray = cur_state_point->getStateBonds();
    int boundSize = bondsArray.size();
    int j = 0;
    for (j = 0; j < boundSize; j++)
    { 
      StateBondP bond = bondsArray[j];
//      if (bond->getBrokenBond() == false)
//      {
        StateMaterialPointP second_point = bond->getSecondPoint();
//        double volume = second_point->getInitialVolume();
//        internal_force = internal_force + bond->getPairwiseForce()*volume*bond->getSurfaceCorrectionFactor();
        StateBondP twin = twinBond(bond);
        Matrix3D twinDeformation = twin->getCenterPoint()->getDeformationGradient();
        Vector3D bondStateForce = bond->getForceVectorState() - (bond->hourglassForceDensity(springConstant, deformation))*c_hg;
        Vector3D twinStateForce = twin->getForceVectorState() - (twin->hourglassForceDensity(springConstant, twinDeformation))*c_hg;
        internal_force = internal_force + (bondStateForce - twinStateForce)
                                           *bond->getVolume();

//        internal_force = internal_force + (bond->getForceVectorState() - twin->getForceVectorState())
//                                           *bond->getVolume();
//        bool condition1 = ((bond->getCenterPoint()->getID() == 700) &&      
//                           (bond->getSecondPoint()->getID() == 613));

//        if (condition1)
//        {
//          std::cout << "Bond(" << bond->getCenterPoint()->getID() << "," 
//                               << bond->getSecondPoint()->getID() << ")" << std::endl;
//          std::cout << "ForceVectorState= " << bond->getForceVectorState() << std::endl;

//          std::cout << "TwinBond(" << twin->getCenterPoint()->getID() << "," 
//                               << twin->getSecondPoint()->getID() << ")" << std::endl;
//          std::cout << "TwinForceVectorState= " << twin->getForceVectorState() << std::endl;

//          std::cout << "Force exerted from point " << bond->getSecondPoint()->getID() 
//                    << " to the point " << bond->getCenterPoint()->getID() << " is= " << std::endl
//                    << bond->getForceVectorState() - twin->getForceVectorState() << std::endl;
//          std::cout << "Current deformation vector of the bond= " << bond->getDeformationVector() << std::endl;
//          std::cout << "Sai of the bond= " << bond->getSai() << std::endl;


//        }          
//      }

    }
    
//    Vector3D zero(0.0, 0.0, 0.0);
//    if (!(cur_peri_point->getBodyForce() == zero))
//    if (cur_peri_point->getID() == size/3 + 1)
//    if (cur_state_point->getID() == 700)    
//      std::cout << "ID= " << cur_peri_point->getID() << " External= " << cur_peri_point->getBodyForce()
//                          << " Internal= " << internal_force << std::endl << std::endl;
//      std::cout << "ID= " << cur_state_point->getID() << " External= " << cur_state_point->getBodyForce()
//                          << " Internal= " << internal_force << std::endl << std::endl;    
//      std::cout << "internalForce= " << internal_force << std::endl << std::endl;

    Vector3D acceleration = (internal_force + cur_state_point->getBodyForce())/(d_density);
    if (cur_state_point->getVelocityBoundaryFlag() == true)
    {
      statePointsArray[i]->setAccele(zero);
    }
    else
    {
      statePointsArray[i]->setAccele(acceleration);
    }


    cur_state_point->setVelNew(cur_state_point->getVelNew() + (cur_state_point->getAccele())*delT);
    cur_state_point->setDisp(cur_state_point->getDisp() + (cur_state_point->getVelNew())*delT);
//    cur_state_point->setPosNew(cur_state_point->getPosNew() + cur_state_point->getDisp());
    cur_state_point->setPosNew(cur_state_point->getPosOld() + cur_state_point->getDisp());

//    if (cur_peri_point->getID() == 1)
//    if (cur_state_point->getID() == 700)
//    {
//      std::cout << "ID= " << cur_state_point->getID() << " Acceleration= " << cur_state_point->getAccele() << std::endl;
//      std::cout << "VelNew= " << cur_state_point->getVelNew() << std::endl;
//      std::cout << "Disp= " << cur_state_point->getDisp() << std::endl;
//      std::cout << "PosNew= " << cur_state_point->getPosNew() << std::endl;
//      std::cout << "PosOld= " << cur_state_point->getPosOld() << std::endl;
//    }
     
  }

  setStatePoints(statePointsArray);

}


StateBondP
Body::twinBond(const StateBondP& bond)
{
  StateMaterialPointP second_point = bond->getSecondPoint();
  StateMaterialPointP center_point = bond->getCenterPoint();
  StateBondPArray bondsArray = second_point->getStateBonds();    
  int boundSize = bondsArray.size();
  int j = 0;
  for (j = 0; j < boundSize; j++)
  { 
    StateBondP twinBond = bondsArray[j];
    if (twinBond->getSecondPoint() == center_point)
    {
      return twinBond;
    }

  }    

}



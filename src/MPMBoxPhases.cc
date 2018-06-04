#include <MPMBoxPhases.h>
#include <BoxSP.h>
#include <Body.h>
#include <ComplicatedGeometrySP.h>
#include <ComplicatedGeometry.h>
#include <MaterialPointP.h>
#include <MaterialPoint.h>
#include <Solver.h>

#include <iostream>
#include <fstream>
#include <string>

#include <sys/stat.h>
#include <unistd.h>

#include <eigen/Eigen/Dense>


using namespace FiniteElement;


MPMBoxPhases::MPMBoxPhases()
  : d_mpm_box(new MPMBox()), d_del_t(0.001), d_time(1.0), d_solver_type(Lumped), d_update_stress_type(USL), d_implicit_parameter(0.0)
{
  d_consistent_mass.clear();
  d_inverse_mass.clear();
  d_lumped_mass.clear();
  d_velocity.clear();
  d_internal_force.clear();
  d_external_force.clear();

  d_consistent_mass.reserve(6000);
  d_inverse_mass.reserve(6000);
  d_lumped_mass.reserve(6000);
  d_velocity.reserve(6000);
  d_internal_force.reserve(6000);
  d_external_force.reserve(6000);

  mpmSolver();

}



MPMBoxPhases::MPMBoxPhases(MPMBoxP box, solverType solver, double time, double totalTime)
  : d_mpm_box(box), d_del_t(time), d_time(totalTime), d_solver_type(solver), d_update_stress_type(USL), d_implicit_parameter(0.0)
{
  d_consistent_mass.clear();
  d_inverse_mass.reserve(6000);
  d_inverse_mass.clear();
  d_lumped_mass.clear();
  d_velocity.clear();
  d_internal_force.clear();
  d_external_force.clear();

  d_consistent_mass.reserve(6000);
  d_lumped_mass.reserve(6000);
  d_velocity.reserve(6000);
  d_internal_force.reserve(6000);
  d_external_force.reserve(6000);

  mpmSolver();

}


MPMBoxPhases::MPMBoxPhases(MPMBoxP box, double parameter, double time, double totalTime)
  : d_mpm_box(box), d_del_t(time), d_time(totalTime), d_solver_type(Lumped), d_update_stress_type(USL), d_implicit_parameter(parameter)
{
  d_consistent_mass.clear();
  d_inverse_mass.reserve(6000);
  d_inverse_mass.clear();
  d_lumped_mass.clear();
  d_velocity.clear();
  d_internal_force.clear();
  d_external_force.clear();

  d_consistent_mass.reserve(6000);
  d_lumped_mass.reserve(6000);
  d_velocity.reserve(6000);
  d_internal_force.reserve(6000);
  d_external_force.reserve(6000);

  mpmSolver();

}


MPMBoxPhases::MPMBoxPhases(MPMBoxP box, solverType solver, double parameter, double time, double totalTime)
  : d_mpm_box(box), d_del_t(time), d_time(totalTime), d_solver_type(solver), d_update_stress_type(USL), d_implicit_parameter(parameter)
{
  d_consistent_mass.clear();
  d_inverse_mass.reserve(6000);
  d_inverse_mass.clear();
  d_lumped_mass.clear();
  d_velocity.clear();
  d_internal_force.clear();
  d_external_force.clear();

  d_consistent_mass.reserve(6000);
  d_lumped_mass.reserve(6000);
  d_velocity.reserve(6000);
  d_internal_force.reserve(6000);
  d_external_force.reserve(6000);

  mpmSolver();

}


MPMBoxPhases::MPMBoxPhases(MPMBoxP box, solverType solver, updateStressType update, double parameter, double time, double totalTime)
  : d_mpm_box(box), d_del_t(time), d_time(totalTime), d_solver_type(solver), d_update_stress_type(update), d_implicit_parameter(parameter)
{
  d_consistent_mass.clear();
  d_inverse_mass.reserve(6000);
  d_inverse_mass.clear();
  d_lumped_mass.clear();
  d_velocity.clear();
  d_internal_force.clear();
  d_external_force.clear();

  d_consistent_mass.reserve(6000);
  d_lumped_mass.reserve(6000);
  d_velocity.reserve(6000);
  d_internal_force.reserve(6000);
  d_external_force.reserve(6000);

  std::cout << "Successfull constructor before mpmSolver();" << std::endl;

  mpmSolver();

}



MPMBoxPhases::~MPMBoxPhases()
{
}

//***********************************************************************************************************

void
MPMBoxPhases::initializationPhaseLumped()
{
//  createNodesConsistentMass();
//  createNodesInverseMass();
  createNodesLumpedMass();
//  std::cout << "createNodesLumpedMass()" << std::endl;
  createNodesVelocityLumped();
//  createNodesVelocityConsistent();
  createNodesMomentum();
    if (d_update_stress_type == USF)
    {
      findPointsStrainIncrementAndUpdateStrain();
      findPointsStressIncrementAndUpdateStress();
    }
  createNodesInternalForce();
  createNodesExternalForce();
  if (d_implicit_parameter != 0.0)
    createNodesIncrementForce();  

//  mpmPrintOld();
//  mpmPrintNew();

}


void
MPMBoxPhases::initializationPhaseConsistent()
{
  createNodesConsistentMass();
//  createNodesInverseMass();
//  createNodesLumpedMass();
//  createNodesVelocityLumped();
  createNodesVelocityConsistent();
  if (d_implicit_parameter != 0.0)
    createNodesMomentum();
    if (d_update_stress_type == USF)
    {
      findPointsStrainIncrementAndUpdateStrain();
      findPointsStressIncrementAndUpdateStress();
    }
  createNodesInternalForce();
  createNodesExternalForce();
  if (d_implicit_parameter != 0.0)
    createNodesIncrementForce();  


//  mpmPrintOld();
//  mpmPrintNew();

}


void
MPMBoxPhases::initializationPhaseLumpedUsingMomentum()
{
//  createNodesConsistentMass();
//  createNodesInverseMass();
  createNodesLumpedMass();
  createNodesVelocityLumped();
//  createNodesVelocityConsistent();
  createNodesMomentum();
    if (d_update_stress_type == USF)
    {
      findPointsStrainIncrementAndUpdateStrain();
      findPointsStressIncrementAndUpdateStress();
    }
  createNodesInternalForce();
  createNodesExternalForce();
  if (d_implicit_parameter != 0.0)
    createNodesIncrementForce();  


//  mpmPrintOld();
//  mpmPrintNew();

}


void
MPMBoxPhases::initializationPhaseConsistentUsingMomentum()
{
  createNodesConsistentMass();
  createNodesInverseMass();
//  createNodesLumpedMass();
//  createNodesVelocityLumped();
  createNodesVelocityConsistent();
  createNodesMomentum();
    if (d_update_stress_type == USF)
    {
      findPointsStrainIncrementAndUpdateStrain();
      findPointsStressIncrementAndUpdateStress();
    }
  createNodesInternalForce();
  createNodesExternalForce();
  if (d_implicit_parameter != 0.0)
    createNodesIncrementForce();  


//  mpmPrintOld();
//  mpmPrintNew();

}


BoxSP
MPMBoxPhases::getGridBox()
{
  if (d_mpm_box->getBoxFlag())
  {
    BoxSP box = d_mpm_box->getBox(); 
    return box;
  }
  else
  {
    ComplicatedGeometrySP geometry = d_mpm_box->getGeometry();
//    std::cout << "GetGridBox() " << std::endl;
    BoxSP box = geometry->getBox();
//    std::cout << "GetGridBox() " << std::endl;
    return box;
  }
}    


std::vector<NodeP>
MPMBoxPhases::getGridNodes()
{
  if (d_mpm_box->getBoxFlag())
  {
    BoxSP box = d_mpm_box->getBox();
    std::vector<NodeP> nodesArray = box->getNodes();
    return nodesArray;
  }  
  else
  {
    ComplicatedGeometrySP geometry = d_mpm_box->getGeometry();
    BoxSP box = geometry->getBox();
    std::vector<NodeP> nodesArray = box->getNodes();
    return nodesArray;
  }  
}


std::vector<SymmetricMaterialTrilinearElementSP>
MPMBoxPhases::getGridElements()
{
  if (d_mpm_box->getBoxFlag())
  {
    BoxSP box = d_mpm_box->getBox();
    std::vector<SymmetricMaterialTrilinearElementSP> elementsArray = box->getElements();
    return elementsArray;
  }  
  else
  {
    ComplicatedGeometrySP geometry = d_mpm_box->getGeometry();
    BoxSP box = geometry->getBox();
    std::vector<SymmetricMaterialTrilinearElementSP> elementsArray = box->getElements();
    return elementsArray;
  }  
}


std::vector<MaterialPointP>
MPMBoxPhases::getGeometryPoints()
{
  if (d_mpm_box->getBoxFlag())
  {
    BoxSP box = d_mpm_box->getBox();
    std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
    return pointsArray;
  }  
  else
  {
    ComplicatedGeometrySP geometry = d_mpm_box->getGeometry();
//    BoxSP box = geometry->getBox();
//    std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
    std::vector<MaterialPointP> pointsArray = geometry->getMaterialPointsArray();
    return pointsArray;
  }  
}



void
MPMBoxPhases::createNodesConsistentMass()
{
//  double sum = 0.0;
  d_consistent_mass.clear();
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();
  int size = nodesArray.size();

  initialMatrix(d_consistent_mass, size);

//  std::cout << size << std::endl;

  for (int i = 0; i < size; i++)
  {
    NodeP cur_node1 = nodesArray[i];
//    std::cout << cur_node1->getID() << std::endl;
    for (int j = 0; j < size; j++)
    {
      NodeP cur_node2 = nodesArray[j];
//      std::cout << cur_node2->getID() << std::endl;
      d_consistent_mass[i][j] = d_mpm_box->consistentMassComponent(cur_node1, cur_node2);
//      sum += d_consistent_mass[i][j];
//      std::cout << d_consistent_mass[i][j] << "  "; 
//      d_consistent_mass[i][j] = 0.0;
    }
   
//    std::cout << std::endl;

  }
  
//  std::cout << sum << std::endl << std::endl;

}


void
MPMBoxPhases::createNodesInverseMass()
{
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();

  int size = nodesArray.size();

//  std::vector<std::vector<double> > inv_mass_matrix;
  inverse(d_consistent_mass, d_inverse_mass, size);

//  for(int m = 0; m < size; m++)
//  {
//    for(int n = 0; n < size; n++)
//      std::cout << d_inverse_mass[m][n] << "  ";
//    std::cout << std::endl;
//  }

//  std::cout << std::endl;

/*  MatrixXd consistMass(size, size);
  MatrixXd invMass(size, size);
  for (int nn = 0; nn < size; nn++)
    for (int mm = 0; mm < size; mm++)
    {
      consistMass(nn, mm) = d_consistent_mass[nn][mm];
      invMass(nn, mm) = d_inverse_mass[nn][mm];
    }

//  std::cout << consistMass*invMass << std::endl << std::endl;
//  std::cout << invMass*consistMass << std::endl << std::endl;*/
}



void
MPMBoxPhases::createNodesLumpedMass()
{
  d_lumped_mass.clear();
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();
  int size = nodesArray.size();

  initialVector(d_lumped_mass, size);

  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];
//    std::cout << "Node's ID: " << cur_node->getID() << "  " << std::endl;
    d_lumped_mass[i] = d_mpm_box->lumpedMassComponent(cur_node);
//    d_lumped_mass.emplace_back(d_mpm_box->lumpedMassComponent(cur_node));
//    std::cout << d_lumped_mass[i] << std::endl;
  }
//    std::cout << "The size of lumped matrix: " << d_lumped_mass.size() << std::endl;
}


void
MPMBoxPhases::createNodesVelocityLumped()
{
  int number = 0;
  d_velocity.clear();
  createNodesLumpedMass();
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();

  int size = nodesArray.size();

//  std::cout << "Size of nodes= " << size << std::endl;

//    initialVector(d_velocity, size);

  Vector3D zero(0.0);
  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];

//    if (d_lumped_mass[i] == 0.0)
//    {
//      number++;
//      std::cout << "ID= " << cur_node->getID() << std::endl;
//    }
   
//    d_velocity[i] = (d_mpm_box->momentumComponent(cur_node))*(1/d_lumped_mass[i]);
//    Vector3D vec = d_mpm_box->momentumComponent(cur_node)*(1/d_lumped_mass[i]);
    if ((nodesBoundaryCondition(cur_node)) || (d_lumped_mass[i] <= MTOL))
    {
      cur_node->velocity(zero);
    }
    else
    {
      d_velocity.emplace_back(d_mpm_box->momentumComponent(cur_node)*(1/d_lumped_mass[i]));
      cur_node->velocity(d_velocity[i]);
    }
//    d_velocity.emplace_back(d_mpm_box->momentumComponent(cur_node)*(1/d_mpm_box->lumpedMassComponent(cur_node)));
//    std::cout << d_velocity[i] << std::endl;
//    if ((d_solver_type == Lumped) || (d_solver_type == LumpedUsingMomentum))
//    if (nodesBoundaryCondition(cur_node))
//    {
//      cur_node->velocity(zero);
//    }
//    else
//    {
//      cur_node->velocity(d_velocity[i]);
//    }
//    std::cout << "NodeID= " << cur_node->getID() << " Node velocity= " << cur_node->getID() << "  " << cur_node->velocity() << std::endl;
   } 

//   std::cout << "Number of nodes with lumped mass equal to zero= " << number << std::endl;
}


void
MPMBoxPhases::createNodesVelocityConsistent()
{
  std::vector<std::vector<double> > massMatrix;
  std::vector<double> momentumNodesVector;
  std::vector<double> nodesVelocity;
  
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();
  int size = nodesArray.size();

  initialMatrix(massMatrix, 3*size);
  initialVector(momentumNodesVector, 3*size);

  Vector3D component(0.0);

  for (int i = 0; i < size; i++)
  {
    int a = 3*i + 1;
//    Vector3D component = d_velocity[i]*(d_lumped_mass[i]);
    component = d_mpm_box->momentumComponent(nodesArray[i]);
    momentumNodesVector[a-1] = component.x();
    momentumNodesVector[a] = component.y();
    momentumNodesVector[a+1] = component.z();
    for (int j = 0; j < size; j++)
    {
      int b = 3*j + 1;
      massMatrix[a-1][b-1] = d_consistent_mass[i][j];
      massMatrix[a][b] = d_consistent_mass[i][j];
      massMatrix[a+1][b+1] = d_consistent_mass[i][j];
    }
    component.set(0.0);
  }
 
   MUMPS(massMatrix, momentumNodesVector, nodesVelocity);

  Vector3D vel(0.0);
  Vector3D zero(0.0);
  for (int k = 0; k < size; k++)
  {
    NodeP cur_node = nodesArray[k];
    vel.x(nodesVelocity[3*k]);
    vel.y(nodesVelocity[3*k+1]);
    vel.z(nodesVelocity[3*k+2]);

//    std::cout << vel << std::endl;
    if (nodesBoundaryCondition(cur_node))
    {
      cur_node->velocity(zero);
    }
    else
    {
      cur_node->velocity(vel);
    }

    vel.set(0.0);
//    std::cout << "Node" << cur_node->getID() << " Velocity: " << cur_node->velocity() << std::endl;
  }  

//  for (int ii = 0; ii < 3*size; ii++)
//  {
//    std::cout << nodesVelocity[ii] << "  ";      
//    std::cout << momentumNodesVector[ii] << "  ";
//    if (ii%3 == 2) std::cout << std::endl;     
//    for (int jj = 0; jj < 3*size; jj++)
//    {
//      std::cout << massMatrix[ii][jj] << "  ";
//    }
//    std::cout //<< "               " << ii%3 
//              << std::endl << std::endl;
//  }


//  for (int jj; jj < size; jj++)
//     std::cout << d_velocity[jj] << std::endl;

}


void
MPMBoxPhases::createNodesMomentum()
{
  
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();

  int size = nodesArray.size();

//    initialVector(d_velocity, size);

Vector3D zero(0.0);

  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];   
//    d_velocity[i] = (d_mpm_box->momentumComponent(cur_node))*(1/d_lumped_mass[i]);
//    Vector3D vec = d_mpm_box->momentumComponent(cur_node)*(1/d_lumped_mass[i]);
//    d_velocity.emplace_back(d_mpm_box->momentumComponent(cur_node)*(1/d_lumped_mass[i]));
//   std::cout << d_velocity[i] << std::endl;

    if (nodesBoundaryCondition(cur_node))
    {
      cur_node->momentum(zero);
    }
    else
    {
    cur_node->momentum(d_mpm_box->momentumComponent(cur_node));
    }

//    std::cout << cur_node->momentum() << std::endl;
  }
}


void
MPMBoxPhases::createNodesInternalForce()
{
  d_internal_force.clear();
  d_mpm_box->updatePointsVolume();
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();

  int size = nodesArray.size();
//  std::cout << "NodeSize= " << size << std::endl;
  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];   
//    d_internal_force[i] = d_mpm_box->internalForce(cur_node);
    d_internal_force.emplace_back(d_mpm_box->internalForce(cur_node));
//    std::cout << d_internal_force[i] << std::endl;
    cur_node->internalForce(d_internal_force[i]);
//    std::cout << cur_node->internalForce() << std::endl;
  }
}


void
MPMBoxPhases::createNodesExternalForce()
{
  std::vector<Vector3D> body_force;
  std::vector<Vector3D> point_force;
  d_external_force.clear();
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();
  int size = nodesArray.size();
  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];   
//    body_force[i] = d_mpm_box->bodyForceComponent(cur_node);
    body_force.emplace_back(d_mpm_box->bodyForceComponent(cur_node));
//    point_force[i] = d_mpm_box->pointForceComponent(cur_node);
    point_force.emplace_back(d_mpm_box->pointForceComponent(cur_node));
    cur_node->bodyForce(body_force[i]);
    cur_node->pointForce(point_force[i]);
//    std::cout << d_internal_force[i] << "   " << body_force[i] + point_force[i] << std::endl;
    cur_node->externalForce(cur_node->bodyForce() + cur_node->pointForce());
//    std::cout << cur_node->externalForce() << std::endl;
//    std::cout << cur_node->pointForce() << std::endl;     
  }

}



void
MPMBoxPhases::createNodesIncrementForce()         //  Impilict
{
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();
  int size = nodesArray.size();
  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];
    cur_node->incrementForce(d_mpm_box->forceIncrementComponent(cur_node)*d_del_t);
//    std::cout << std::endl; 
//    Matrix3D mat = cur_node->incrementForce();
//    mat.printMatrix();
  }
}

//***********************************************************************************************************

void
MPMBoxPhases::lagrangianPhaseLumped()
{
//    if (d_update_stress_type == USF)
//    {
//      findPointsStrainIncrementAndUpdateStrain();
//      findPointsStressIncrementAndUpdateStress();
//    }
//  findNodesVelocityConsistentMethod();
//  if (d_implicit_parameter == 0.0)
//  {

    findNodesLumpedAcceleration();
//    std::cout << "findNodesLumpedAcceleration()" << std::endl;
//    findNodesConsistentAcceleration();  
    updateNodesVelocity();
//    std::cout << "updateNodesVelocity()" << std::endl;
    updateNodesMomentum();
//    std::cout << "updateNodesMomentum()" << std::endl;
//  }
  if (d_implicit_parameter != 0.0)
  {
    findNodesVelocityImplicitMethod();
    updateAcceleration();
  }  
//  updateNodesMomentum();
  updatePointsVelocityAndPosition();
//  std::cout << "updatePointsVelocityAndPosition()" << std::endl;
//  updatePointsVelocityAndPositionLumpedUsingMomentum();
//  updatePointsVelocityAndPositionConsistentUsingMomentum();
    if (d_update_stress_type != USF)
    {
      findPointsStrainIncrementAndUpdateStrain();
      findPointsStressIncrementAndUpdateStress();
    }

//  mpmPrintOld();
//  mpmPrintNew();

}


void
MPMBoxPhases::lagrangianPhaseConsistent()
{
//    if (d_update_stress_type == USF)
//    {
//      findPointsStrainIncrementAndUpdateStrain();
//      findPointsStressIncrementAndUpdateStress();
//    }
//  findNodesLumpedAcceleration();
  findNodesConsistentAcceleration();  
  updateNodesVelocity();
  if (d_implicit_parameter != 0.0)
  {
    updateNodesMomentum();
    findNodesVelocityImplicitMethod();
    updateAcceleration();
  }  

//  updateNodesMomentum();
  updatePointsVelocityAndPosition();
//  updatePointsVelocityAndPositionLumpedUsingMomentum();
//  updatePointsVelocityAndPositionConsistentUsingMomentum();
    if (d_update_stress_type != USF)
    {
      findPointsStrainIncrementAndUpdateStrain();
      findPointsStressIncrementAndUpdateStress();
    }

//  mpmPrintOld();
//  mpmPrintNew();

}


void
MPMBoxPhases::lagrangianPhaseLumpedUsingMomentum()
{
//  findNodesLumpedAcceleration();
//  findNodesConsistentAcceleration();  
//  updateNodesVelocity();
  if (d_implicit_parameter == 0.0)
  {
//    if (d_update_stress_type == USF)
//    {
//      findPointsStrainIncrementAndUpdateStrain();
//      findPointsStressIncrementAndUpdateStress();
//    }
    updateNodesMomentum();
//    updatePointsVelocityAndPosition();
    updatePointsVelocityAndPositionLumpedUsingMomentum();
//    updatePointsVelocityAndPositionConsistentUsingMomentum();
    updateNodesVelocityLumpedUsingMomentum();
    if (d_update_stress_type != USF)
    {
      findPointsStrainIncrementAndUpdateStrain();
      findPointsStressIncrementAndUpdateStress();
    }
  }
  else
  {
//    if (d_update_stress_type == USF)
//    {
//      findPointsStrainIncrementAndUpdateStrain();
//      findPointsStressIncrementAndUpdateStress();
//    }
    findNodesLumpedAcceleration();
    updateNodesVelocity();
    updateNodesMomentum();
    findNodesVelocityImplicitMethod();
    findNodesTotalForce();
    updateNodesMomentumImplicit();
    updatePointsVelocityAndPositionLumpedUsingMomentum();
    updateNodesVelocityLumpedUsingMomentum();
    if (d_update_stress_type != USF)
    {
      findPointsStrainIncrementAndUpdateStrain();
      findPointsStressIncrementAndUpdateStress();
    }
  }
//  mpmPrintOld();
//  mpmPrintNew();

}


void
MPMBoxPhases::lagrangianPhaseConsistentUsingMomentum()
{
  if (d_implicit_parameter == 0)
  {
//    if (d_update_stress_type == USF)
//    {
//      findPointsStrainIncrementAndUpdateStrain();
//      findPointsStressIncrementAndUpdateStress();
//    }
//    findNodesLumpedAcceleration();
//    findNodesConsistentAcceleration();  
//    updateNodesVelocity();
    updateNodesMomentum();
//    updatePointsVelocityAndPosition();
//    updatePointsVelocityAndPositionLumpedUsingMomentum();
    updatePointsVelocityAndPositionConsistentUsingMomentum();
    updateNodesVelocityConsistentUsingMomentum();
    if (d_update_stress_type != USF)
    {
      findPointsStrainIncrementAndUpdateStrain();
      findPointsStressIncrementAndUpdateStress();
    }
  }
  else
  {
//    if (d_update_stress_type == USF)
//    {
//      findPointsStrainIncrementAndUpdateStrain();
//      findPointsStressIncrementAndUpdateStress();
//    }
    findNodesConsistentAcceleration();
    updateNodesVelocity();
    updateNodesMomentum();
    findNodesVelocityImplicitMethod();
    findNodesTotalForce();
    updateNodesMomentumImplicit();
    updatePointsVelocityAndPositionConsistentUsingMomentum();
    updateNodesVelocityConsistentUsingMomentum();
    if (d_update_stress_type != USF)
    {
      findPointsStrainIncrementAndUpdateStrain();
      findPointsStressIncrementAndUpdateStress();
    }
  }

//  mpmPrintOld();
//  mpmPrintNew();

}


void
MPMBoxPhases::findNodesLumpedAcceleration()
{
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();

  int size = nodesArray.size();
  Vector3D zero(0.0);
  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];

//    std::cout << "d_lumped_mass[" << i << "]= " << d_lumped_mass[i] << std::endl;
//    if (d_implicit_parameter == 0)
//    { 

    if ((nodesBoundaryCondition(cur_node)) || (d_lumped_mass[i] <= MTOL))
    {
      cur_node->acceleration(zero);
    }
    else
    {
      Vector3D intForce = cur_node->internalForce();
      Vector3D extForce = cur_node->externalForce();
      Vector3D force = intForce + extForce;
//      std::cout << "Force= " << force << std::endl;
      Vector3D accele = force/d_lumped_mass[i];
      cur_node->acceleration(accele);
//      cur_node->acceleration((cur_node->internalForce() + cur_node->externalForce())/(d_lumped_mass[i]));
    }    

//    }
//    else
//    {
//      cur_node->acceleration((cur_node->newVelocityImplicit() - cur_node->velocity())*(1/d_del_t));
//    }
//      std::cout << "NodeID= " << cur_node->getID() << " Accele= " 
//                << cur_node->acceleration().x() << ", " << cur_node->acceleration().y() << ", "
//                << cur_node->acceleration().z() << std::endl;
  }

}


void
MPMBoxPhases::findNodesConsistentAcceleration()
{
  std::vector<std::vector<double> > massMatrix;
  std::vector<double> forceNodesVector;
  std::vector<double> nodesAcceleration;
  
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();

  int size = nodesArray.size();

  initialMatrix(massMatrix, 3*size);
  initialVector(forceNodesVector, 3*size);

  Vector3D component(0.0);

  for (int i = 0; i < size; i++)
  {
    int a = 3*i + 1;
//    Vector3D component = d_velocity[i]*(d_lumped_mass[i]);
    component = nodesArray[i]->internalForce() + nodesArray[i]->externalForce();
    forceNodesVector[a-1] = component.x();
    forceNodesVector[a] = component.y();
    forceNodesVector[a+1] = component.z();
    for (int j = 0; j < size; j++)
    {
      int b = 3*j + 1;
      massMatrix[a-1][b-1] = d_consistent_mass[i][j];
      massMatrix[a][b] = d_consistent_mass[i][j];
      massMatrix[a+1][b+1] = d_consistent_mass[i][j];
    }
    component.set(0.0);
  }
 
   MUMPS(massMatrix, forceNodesVector, nodesAcceleration);

  Vector3D acc(1.0);
  Vector3D zero(0.0);
  for (int k = 0; k < size; k++)
  {
    NodeP cur_node = nodesArray[k];
    acc.x(nodesAcceleration[3*k]);
    acc.y(nodesAcceleration[3*k+1]);
    acc.z(nodesAcceleration[3*k+2]);

//    std::cout << vel << std::endl;
    if (nodesBoundaryCondition(cur_node))
    {
      cur_node->acceleration(zero);
    }
    else
    {
      cur_node->acceleration(acc);
    } 

//    std::cout << "Node" << cur_node->getID() << " Acceleration: " << cur_node->acceleration() << std::endl;
  }  

//  for (int ii = 0; ii < 3*size; ii++)
//  {
//       std::cout << nodesVelocity[ii] << "  ";      
//      std::cout << momentumNodesVector[ii] << "  ";
//      if (ii%3 == 2) std::cout << std::endl;     
//    for (int jj = 0; jj < 3*size; jj++)
//    {
//      std::cout << massMatrix[ii][jj] << "  ";
//    }
//  std::cout << std::endl;
//  }


//  for (int jj; jj < size; jj++)
//     std::cout << d_velocity[jj] << std::endl;

}



void
MPMBoxPhases::updateNodesVelocity()
{
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();

  int size = nodesArray.size();
  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i]; 
    cur_node->newVelocity(cur_node->velocity() + cur_node->acceleration()*d_del_t);
//    std::cout << "ID= " << cur_node->getID() << "  " << cur_node->velocity() << "   " << cur_node->newVelocity() << std::endl;
  }
    
}


void
MPMBoxPhases::updateNodesMomentum()
{
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();

  int size = nodesArray.size();
  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];
//    if (d_implicit_parameter == 0)
//    { 
    if (nodesBoundaryCondition(cur_node))
    {
      cur_node->newMomentum(cur_node->momentum());
    }
    else
    {
      cur_node->newMomentum(cur_node->momentum() + (cur_node->internalForce()+cur_node->externalForce())*d_del_t);
    } 

//    }
//    else
//    {
//      cur_node->newMomentum(cur_node->momentum() + (cur_node->totalForceSemiImplicit())*d_del_t);
//    }
//    std::cout << cur_node->momentum() << "   " << cur_node->newMomentum() << std::endl;
//    std::cout << cur_node->velocity() << "   " << cur_node->newVelocity() << std::endl;
  }
    
}


void
MPMBoxPhases::updateNodesMomentumImplicit()
{
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
//  std::vector<NodeP> nodesArray = getGridNodes();

  int size = nodesArray.size();
  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];
//    if (d_implicit_parameter == 0)
//    { 
//      cur_node->newMomentum(cur_node->momentum() + (cur_node->internalForce()+cur_node->externalForce())*d_del_t);
//    }
//    else
//    {
    if (nodesBoundaryCondition(cur_node))
    {
      cur_node->newMomentum(cur_node->momentum());
    }
    else
    {
     cur_node->newMomentum(cur_node->momentum() + (cur_node->totalForceSemiImplicit())*d_del_t);
    } 
 
//    }
//    std::cout << cur_node->momentum() << "   " << cur_node->newMomentum() << std::endl;
//    std::cout << cur_node->velocity() << "   " << cur_node->newVelocity() << std::endl;
  }
    
}




void
MPMBoxPhases::updatePointsVelocityAndPosition()
{

  Vector3D accSum(0.0);
  Vector3D velSum(0.0);
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
//  std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
  std::vector<MaterialPointP> pointsArray = getGeometryPoints();
  int pointSize = pointsArray.size();
//  std::cout << "pointSize= " << pointSize << std::endl;
  for (int k = 0; k < pointSize; k++)
  {
    MaterialPointP cur_point = pointsArray[k];
//    std::cout << "ID= " << cur_point->getID() << " Volume= " << cur_point->getInitialVolume() << std::endl;
    std::vector<NodeP> nodesArray = box->getNodes(); 
    int size = nodesArray.size();
//    std::cout << "nodeSize= " << size << std::endl;
    accSum.reset();
    velSum.reset();
    for (int i = 0; i < size; i++)
    { 
      NodeP cur_node = nodesArray[i];
//      if (k == 0)
//      std::cout << "acceleration " << cur_node->acceleration() << "  velocity " << cur_node->newVelocity() << std::endl;
//      if (d_solver_type == ConsistentUsingMomentum)
//      {
//        for (int j = 0; j < size; j++)
//        {
//          NodeP cur_node2 = nodesArray[j];
//          accSum = accSum + (cur_node2->internalForce()+cur_node2->externalForce())*
//                            (d_lumped_mass[j]*d_mpm_box->interpolate(cur_node2, cur_point)); 
//          velSum = velSum + cur_node2->newMomentum()*
//                            (d_lumped_mass[j]*d_mpm_box->interpolate(cur_node2, cur_point));
//        }
        
//      }
//      else
//      { 
//       std::cout << "nodeAcceleration= " << cur_node->acceleration() << std::endl;
       accSum = accSum + cur_node->acceleration()*d_mpm_box->interpolate(cur_node, cur_point);
       velSum = velSum + cur_node->newVelocity()*d_mpm_box->interpolate(cur_node, cur_point);
//       std::cout << "accSum= " << accSum << std::endl; 
//      }      
     }

//    std::cout << "accSum= " << accSum << "  velSum= " << velSum << std::endl;
    updateFunctions(cur_point, accSum, velSum);

    nodesArray.clear();

   }

}




void
MPMBoxPhases::updatePointsVelocityAndPositionConsistentUsingMomentum()
{

  Vector3D accSum(0.0);
  Vector3D velSum(0.0);
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<MaterialPointP> pointsArray = getGeometryPoints();
//  std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
  int pointSize = pointsArray.size();
  for (int k = 0; k < pointSize; k++)
  {
    MaterialPointP cur_point = pointsArray[k];
    std::vector<NodeP> nodesArray = box->getNodes(); 
    int size = nodesArray.size();
    accSum.reset();
    velSum.reset();
    for (int i = 0; i < size; i++)
    { 
      NodeP cur_node1 = nodesArray[i];
//      if (d_solver_type == ConsistentUsingMomentum)
//      {
        for (int j = 0; j < size; j++)
        {
          NodeP cur_node2 = nodesArray[j];
          accSum = accSum + (cur_node2->internalForce()+cur_node2->externalForce())*
                            //(d_inverse_mass[i][j]*d_mpm_box->interpolate(cur_node2, cur_point));
                            (d_inverse_mass[i][j]*d_mpm_box->interpolate(cur_node1, cur_point)); 
          velSum = velSum + cur_node2->newMomentum()*
                            //(d_inverse_mass[i][j]*d_mpm_box->interpolate(cur_node2, cur_point));
                            (d_inverse_mass[i][j]*d_mpm_box->interpolate(cur_node1, cur_point));
        }
        
//      }
//      else
//      { 
//        accSum = accSum + cur_node->acceleration()*d_mpm_box->interpolate(cur_node, cur_point);
//        velSum = velSum + cur_node->newVelocity()*d_mpm_box->interpolate(cur_node, cur_point);
//      }      
     }

    updateFunctions(cur_point, accSum, velSum);

    nodesArray.clear();

  }

}


void
MPMBoxPhases::updatePointsVelocityAndPositionLumpedUsingMomentum()
{

  Vector3D accSum(0.0);
  Vector3D velSum(0.0);
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<MaterialPointP> pointsArray = getGeometryPoints();
//  std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
  int pointSize = pointsArray.size();
  for (int k = 0; k < pointSize; k++)
  {
    MaterialPointP cur_point = pointsArray[k];
    std::vector<NodeP> nodesArray = box->getNodes(); 
    int size = nodesArray.size();
    accSum.reset();
    velSum.reset();
//    for (int i = 0; i < size; i++)
//    { 
//      NodeP cur_node = nodesArray[i];
//      if (d_solver_type == ConsistentUsingMomentum)
//      {
        for (int j = 0; j < size; j++)
        {
          NodeP cur_node2 = nodesArray[j];
          if (d_implicit_parameter == 0)
          {
            accSum = accSum + (cur_node2->internalForce()+cur_node2->externalForce())*
                              ((1/d_lumped_mass[j])*d_mpm_box->interpolate(cur_node2, cur_point)); 
            velSum = velSum + cur_node2->newMomentum()*
                              ((1/d_lumped_mass[j])*d_mpm_box->interpolate(cur_node2, cur_point));
          }
          else
          {
            accSum = accSum + (cur_node2->totalForceSemiImplicit())*
                              ((1/d_lumped_mass[j])*d_mpm_box->interpolate(cur_node2, cur_point)); 
            velSum = velSum + cur_node2->newMomentum()*
                              ((1/d_lumped_mass[j])*d_mpm_box->interpolate(cur_node2, cur_point));
          }
        }
        
//      }
//      else
//      { 
//        accSum = accSum + cur_node->acceleration()*d_mpm_box->interpolate(cur_node, cur_point);
//        velSum = velSum + cur_node->newVelocity()*d_mpm_box->interpolate(cur_node, cur_point);
//      }      
//     }

    updateFunctions(cur_point, accSum, velSum);

    nodesArray.clear();

   }

}


void
MPMBoxPhases::updateFunctions(MaterialPointP& cur_point, const Vector3D& accSum, const Vector3D& velSum)
{
    Vector3D zero(0.0);
    cur_point->setVelNew(cur_point->getVelOld() + accSum*d_del_t);

//    std::cout << "Node " << cur_point->getID() << " Vel: " << cur_point->getVelOld() << "  " 
//                                      << cur_point->getVelNew() << std::endl;

    if (!(cur_point->getHasBoundaryCondition()))
    {
//      cur_point->setPosNew(cur_point->getPosOld() + cur_point->getVelNew()*d_del_t);
      cur_point->setPosNew(cur_point->getPosOld() + velSum*d_del_t);
//      std::cout << "PointID= " << cur_point->getID() << " Pos= " << cur_point->getPosNew() << std::endl;
      cur_point->setHasBoundaryCondition(false);
    }
    else
    {
      cur_point->setPosNew(cur_point->getPosOld());
      cur_point->setHasBoundaryCondition(true);
//      std::cout << cur_point->getPosOld() << std::endl;
    }


    cur_point->setDisp(cur_point->getPosNew() - cur_point->getPosOld());


    if (cur_point->getNoneZeroExtForce())
    {
      if (d_mpm_box->getBoxFlag())
      {
        cur_point->setPointForce(d_mpm_box->getBox()->getPointForce());
      }
      else
      {
        cur_point->setPointForce(d_mpm_box->getGeometry()->getPointForce());  
      }      
      cur_point->setNoneZeroExtForce(true);
    }   
    else
    {
      cur_point->setPointForce(zero);
      cur_point->setNoneZeroExtForce(false);
    }      

//    std::cout << "Pos: " << cur_point->getPosOld() << "  " << cur_point->getPosNew() << std::endl << std::endl;
}


void
MPMBoxPhases::updateNodesVelocityLumpedUsingMomentum()
{
  d_velocity.clear();
  createNodesLumpedMass();
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
  int size = nodesArray.size();

//    initialVector(d_velocity, size);

  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];   
//    d_velocity[i] = (d_mpm_box->momentumComponent(cur_node))*(1/d_lumped_mass[i]);
//    Vector3D vec = d_mpm_box->momentumComponent(cur_node)*(1/d_lumped_mass[i]);
//    d_velocity.emplace_back(d_mpm_box->momentumComponent(cur_node)*(1/d_lumped_mass[i]));
//    d_velocity.emplace_back(d_mpm_box->momentumComponent(cur_node)*(1/d_mpm_box->lumpedMassComponent(cur_node)));
//    std::cout << d_velocity[i] << std::endl;
//    if ((d_solver_type == Lumped) || (d_solver_type == LumpedUsingMomentum))
    Vector3D zero(0.0);
    if (nodesBoundaryCondition(cur_node))
    {
      cur_node->newVelocity(zero);
    }
    else
    {
      cur_node->newVelocity(nodesArray[i]->newMomentum()*(1/d_lumped_mass[i]));
    }  

  }
}



void 
MPMBoxPhases::updateNodesVelocityConsistentUsingMomentum()
{
  std::vector<std::vector<double> > massMatrix;
  std::vector<double> momentumNodesVector;
  std::vector<double> nodesVelocity;
  
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
  int size = nodesArray.size();

  initialMatrix(massMatrix, 3*size);
  initialVector(momentumNodesVector, 3*size);

  Vector3D component(0.0);

  for (int i = 0; i < size; i++)
  {
    int a = 3*i + 1;
//    Vector3D component = d_velocity[i]*(d_lumped_mass[i]);
    component = d_mpm_box->newMomentumComponent(nodesArray[i]);
    momentumNodesVector[a-1] = component.x();
    momentumNodesVector[a] = component.y();
    momentumNodesVector[a+1] = component.z();
    for (int j = 0; j < size; j++)
    {
      int b = 3*j + 1;
      massMatrix[a-1][b-1] = d_consistent_mass[i][j];
      massMatrix[a][b] = d_consistent_mass[i][j];
      massMatrix[a+1][b+1] = d_consistent_mass[i][j];
    }
    component.set(0.0);
    a = 0;
  }
 
   MUMPS(massMatrix, momentumNodesVector, nodesVelocity);

  Vector3D vel(0.0);
  Vector3D zero(0.0);
  for (int k = 0; k < size; k++)
  {
    NodeP cur_node = nodesArray[k];
    vel.x(nodesVelocity[3*k]);
    vel.y(nodesVelocity[3*k+1]);
    vel.z(nodesVelocity[3*k+2]);

//    std::cout << vel << std::endl;
    if (nodesBoundaryCondition(cur_node))
    {
      cur_node->newVelocity(zero);
    }
    else
    {
      cur_node->newVelocity(vel);
    }  


    vel.set(0.0);
//    std::cout << "Node" << cur_node->getID() << " Velocity: " << cur_node->velocity() << " New Velocity: " << cur_node->newVelocity() << std::endl;
  }  

//  for (int ii = 0; ii < 3*size; ii++)
//  {
//       std::cout << nodesVelocity[ii] << "  ";      
//      std::cout << momentumNodesVector[ii] << "  ";
//      if (ii%3 == 2) std::cout << std::endl;     
//    for (int jj = 0; jj < 3*size; jj++)
//    {
//      std::cout << massMatrix[ii][jj] << "  ";
//    }
//  std::cout << std::endl;
//  }


//  for (int jj; jj < size; jj++)
//     std::cout << d_velocity[jj] << std::endl;

}



void
MPMBoxPhases::findPointsStrainIncrementAndUpdateStrain()
{
  Point3D tempPosition(0.0, 0.0, 0.0);
  Matrix3D sum(0.0);
  Matrix3D gradientDeformation(0.0);
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<MaterialPointP> pointsArray = getGeometryPoints();
//  std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();

  int pointSize = pointsArray.size();
  for (int k = 0; k < pointSize; k++)
  {
    MaterialPointP cur_point = pointsArray[k];
    if (d_update_stress_type == MUSL)
    {
      tempPosition = cur_point->getPosOld();
      cur_point->setPosOld(cur_point->getPosNew());
    }
   
    std::vector<NodeP> nodesArray = box->getNodes(); 
    int size = nodesArray.size();
    sum.set(0.0);
//    Matrix3D tens;
    Vector3D velocity(0.0);
    for (int i = 0; i < size; i++)
    {
      NodeP cur_node = nodesArray[i];
      velocity.set(0.0);
      if (d_update_stress_type == USF)
        velocity = cur_node->velocity();
//         Matrix3D tens(d_mpm_box->gradientInterpolate(cur_node, cur_point), cur_node->oldVelocity());
      else
        velocity = cur_node->newVelocity(); 
//      Matrix3D tens(d_mpm_box->gradientInterpolate(cur_node, cur_point), cur_node->newVelocity());
      Matrix3D tens(d_mpm_box->gradientInterpolate(cur_node, cur_point), velocity);
      gradientDeformation = gradientDeformation + tens; 
      sum = sum + (tens + tens.Transpose());
      tens.set(0.0);
    }
    cur_point->setDeformationGradientNew(cur_point->getDeformationGradientOld() + 
                                        (gradientDeformation*cur_point->getDeformationGradientOld())*d_del_t);
    cur_point->setStrainIncrement(sum*(d_del_t/2));
    cur_point->setStrainNew(cur_point->getStrainOld() + cur_point->getStrainIncrement());
/*    if (k == pointSize/2)
    {
      std::cout << "Point " << cur_point->getID() << " " << std::endl;
      cur_point->getDeformationGradientOld().printMatrix();
      std::cout << cur_point->getDeformationGradientOld().Determinant();
      std::cout << std::endl;  
      cur_point->getDeformationGradientNew().printMatrix();
      std::cout << cur_point->getDeformationGradientNew().Determinant();   
//      cur_point->getStrainOld().printMatrix();
//      std::cout << std::endl;
//      cur_point->getStrainIncrement().printMatrix();
      std::cout << std::endl;
//      cur_point->getStrainNew().printMatrix();
    } */ 


    nodesArray.clear();
    if (d_update_stress_type == MUSL)
    {
      cur_point->setPosOld(tempPosition);
      tempPosition.set(0.0);
    }
    
  }

}



void
MPMBoxPhases::findPointsStressIncrementAndUpdateStress()     // constitutive law
{
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<MaterialPointP> pointsArray = getGeometryPoints();
//  std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
  int pointSize = pointsArray.size();
  for (int k = 0; k < pointSize; k++)
  {
    MaterialPointP cur_point = pointsArray[k];
    cur_point->setStressIncrement(box->getElements()[0]->tensor().doubleContract(cur_point->getStrainIncrement()));
    cur_point->setStressNew(cur_point->getStressOld() + cur_point->getStressIncrement());

/*    if (k == pointSize/2)
    {
      std::cout << "Point " << cur_point->getID() << " " << std::endl;    
      cur_point->getStressOld().printMatrix();
      std::cout << std::endl;
      cur_point->getStrainIncrement().printMatrix();
      std::cout << std::endl;
      cur_point->getStressIncrement().printMatrix();
      std::cout << std::endl;
      cur_point->getStressNew().printMatrix();
    }  */
  }


}



void
MPMBoxPhases::findNodesVelocityImplicitMethod()
{
  std::vector<std::vector<double> > massMatrixImplicit;
  std::vector<double> momentumNodesLumped;
  std::vector<double> nodesVelocityImplicit;
  
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
  int size = nodesArray.size();

  initialMatrix(massMatrixImplicit, 3*size);
  initialVector(momentumNodesLumped, 3*size);

  Vector3D component(0.0);

  for (int i = 0; i < size; i++)
  {
//    int a = 3*i + 1;
    int a = 3*i;
//    Vector3D component = d_velocity[i]*(d_lumped_mass[i]);
//    component = nodesArray[i]->newVelocity()*(d_lumped_mass[i]);
    component = nodesArray[i]->newMomentum();
//    std::cout << component << std::endl;
    Matrix3D deltaForce = nodesArray[i]->incrementForce();
//    deltaForce.printMatrix();
    momentumNodesLumped[a] = component.x();
    momentumNodesLumped[a+1] = component.y();
    momentumNodesLumped[a+2] = component.z();
    if ((d_solver_type == Lumped) || (d_solver_type == LumpedUsingMomentum))
    {
      for (int j = 0; j < 3; j++)
      {
        for (int k = 0; k < 3; k++)
//        int b = 3*j + 1;
          massMatrixImplicit[a+j][a+k] = delta(j, k)*d_lumped_mass[i] - d_implicit_parameter*d_del_t*deltaForce.get(j,k)/2;
//          massMatrixImplicit[a][b] = d_consistent_mass[i][j];
//          massMatrixImplicit[a+1][b+1] = d_consistent_mass[i][j];
      }
    }
    else
    {
      for (int jj = 0; jj < size; jj++)
      {
        int b = 3*jj;
        for (int j = 0; j < 3; j++)
        {
          for (int k = 0; k < 3; k++)
            massMatrixImplicit[a+j][b+k] = delta(j, k)*
                                           d_consistent_mass[i][jj] - d_implicit_parameter*d_del_t*deltaForce.get(j,k)/2;
        }
      }


    }
    component.set(0.0);
    a = 0;
  }

//  int massSize = massMatrixImplicit.size();
//  for (int mm = 0; mm < massSize; mm++)
//    for (int nn = 0; nn < massSize; nn++)
//    {
//      std::cout  << massMatrixImplicit[mm][nn] << " ";
//      if (nn == massSize-1)
//         std::cout << std::endl;
//    } 


//  int momentSize = momentumNodesLumped.size();
//  for (int mm = 0; mm < momentSize; mm++)
//    {
//      std::cout  << momentumNodesLumped[mm] << " ";
//    } 
//  std::cout << std::endl;

   MUMPS(massMatrixImplicit, momentumNodesLumped, nodesVelocityImplicit);

  Vector3D vel(0.0);
  Vector3D zero(0.0);
  for (int k = 0; k < size; k++)
  {
    NodeP cur_node = nodesArray[k];
    vel.x(nodesVelocityImplicit[3*k]);
    vel.y(nodesVelocityImplicit[3*k+1]);
    vel.z(nodesVelocityImplicit[3*k+2]);

//    std::cout << vel << std::endl;
    if (nodesBoundaryCondition(cur_node))
    {
      cur_node->newVelocityImplicit(zero);
    }
    else
    {
      cur_node->newVelocityImplicit(vel);
    }  

    vel.set(0.0);
//    std::cout << "Node" << cur_node->getID() << 
//                 " New Velocity: " << cur_node->newVelocity() <<
//                 " New Velocity Implicit: " << cur_node->newVelocityImplicit() << std::endl;
//  }

  

//  for (int ii = 0; ii < 3*size; ii++)
//  {
//       std::cout << nodesVelocity[ii] << "  ";      
//      std::cout << momentumNodesVector[ii] << "  ";
//      if (ii%3 == 2) std::cout << std::endl;     
//    for (int jj = 0; jj < 3*size; jj++)
//    {
//      std::cout << massMatrix[ii][jj] << "  ";
//    }
//  std::cout << std::endl;
  }


//  for (int jj; jj < size; jj++)
//     std::cout << d_velocity[jj] << std::endl;


}





void
MPMBoxPhases::updateAcceleration()
{
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
  int size = nodesArray.size();
  Vector3D zero(0.0);
  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];
    if (nodesBoundaryCondition(cur_node))
    {
      cur_node->acceleration(zero);
      cur_node->newVelocity(zero);
    }
    else
    {
      cur_node->acceleration((cur_node->newVelocityImplicit() - cur_node->velocity())*(1/d_del_t));
      cur_node->newVelocity(cur_node->newVelocityImplicit());
    }  
  
  }
}


int
MPMBoxPhases::delta(int i, int j)
{
  if (i == j) return 1;
  else return 0;
}



void
MPMBoxPhases::findNodesTotalForce()
{
//  BoxSP box = d_mpm_box->getBox();
  BoxSP box = getGridBox();
  std::vector<NodeP> nodesArray = box->getNodes();
  int size = nodesArray.size();
  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];
    Vector3D updatedForce(0.0);
    updatedForce = d_mpm_box->forceIncrementComponent(cur_node)*cur_node->newVelocityImplicit();
//    std::cout << cur_node->newVelocityImplicit() << std::endl;
    updatedForce = updatedForce*(d_del_t/2);
//    std::cout << updatedForce << std::endl;
    updatedForce = updatedForce + (cur_node->internalForce() + cur_node->externalForce());
//    std::cout << updatedForce << std::endl;
//    std::cout << cur_node->internalForce() << std::endl;
    cur_node->totalForceSemiImplicit(updatedForce*d_implicit_parameter + (cur_node->internalForce() + cur_node->externalForce())*(1 - d_implicit_parameter));
//    Vector3D force(0.0);
//    force = (cur_node->newVelocityImplicit() - cur_node->velocity())/d_del_t;
//    std::cout << "Acceleration " << force << std::endl;
//    force = force*d_lumped_mass[i];
//    std::cout << cur_node->totalForceSemiImplicit() << "  " << force << std::endl;
//    std::cout << d_lumped_mass[i] << std::endl;
//    cur_node->newVelocity(cur_node->newVelocityImplicit());
      
  }

}

//***********************************************************************************************************

void
MPMBoxPhases::convectivePhase()
{
  updatePointsProperties();
  redefineBackgroundGrid();

}



void
MPMBoxPhases::updatePointsProperties()
{
//  BoxSP box = d_mpm_box->getBox();
//  std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
  std::vector<MaterialPointP> pointsArray = getGeometryPoints();
  int pointSize = pointsArray.size();
   
//  int i = 20;
//  std::cout << pointsArray[i]->getPosOld().x() << ",  "
//            << pointsArray[i]->getPosOld().y() << ",  "
//            << pointsArray[i]->getPosOld().z() << std::endl; 

//  std::cout << pointsArray[i]->getMass() << std::endl; 

//  updateBoundaryConditionAndPointsForce();

  for (int k = 0; k < pointSize; k++)
  {
    MaterialPointP cur_point = pointsArray[k];
    cur_point->setPosOld(cur_point->getPosNew());  
    cur_point->setVelOld(cur_point->getVelNew());

    cur_point->setDeformationGradientOld(cur_point->getDeformationGradientNew());
  
    cur_point->setStrainOld(cur_point->getStrainNew());  
    cur_point->setStressOld(cur_point->getStressNew());
/*    if (k == pointSize/2)
    {
      std::cout << "Point " << cur_point->getID() << std::endl; 
      std::cout << cur_point->getVelOld() << "  " << cur_point->getVelNew() << std::endl;
      std::cout << cur_point->getPosOld() << "  " << cur_point->getPosNew() << std::endl;
    }  */ 

  }


//  int i = 1;
//  std::cout << pointsArray[i]->getPosOld().x() << ",  "
//            << pointsArray[i]->getPosOld().y() << ",  "
//            << pointsArray[i]->getPosOld().z() << std::endl; 
            

}




void
MPMBoxPhases::redefineBackgroundGrid()
{

  int extra = 2;
  
  std::vector<double> pointsXComponent;
  std::vector<double> pointsYComponent;
  std::vector<double> pointsZComponent;
//  BoxSP box = d_mpm_box->getBox();
//  std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
  std::vector<MaterialPointP> pointsArray = getGeometryPoints();
  int pointSize = pointsArray.size();
  for (int k = 0; k < pointSize; k++)
  {
    MaterialPointP cur_point = pointsArray[k];
    pointsXComponent.emplace_back(cur_point->getPosOld().x());
    pointsYComponent.emplace_back(cur_point->getPosOld().y());   
    pointsZComponent.emplace_back(cur_point->getPosOld().z()); 

  }

  Vector3D min(findMin(pointsXComponent), findMin(pointsYComponent), findMin(pointsZComponent));
  Vector3D max(findMax(pointsXComponent), findMax(pointsYComponent), findMax(pointsZComponent));
//  std::cout << d_mpm_box->getBoxFlag() << std::endl;
  BoxSP box = getGridBox();

  double dx = (max.x() - min.x())/(box->getNumElementsX());
  double dy = (max.y() - min.y())/(box->getNumElementsY());
  double dz = (max.z() - min.z())/(box->getNumElementsZ());

  min.x(min.x() - extra*dx);
  min.y(min.y() - extra*dy);
  min.z(min.z() - extra*dz);

  max.x(max.x() + extra*dx);
  max.y(max.y() + extra*dy);
  max.z(max.z() + extra*dz);


//  std::cout << "Num of Elements= " << box->getNumElementsY() << std::endl;
//  std::cout << "box " << std::endl;
//  std::cout << "min= " << min << " max= " << max << std::endl;
//  std::cout << box->getNumElementsX() << std::endl;  
  box->setMinPoint(min);    
  box->setMaxPoint(max);
//  std::cout << box->getMinPoint() << std::endl;
//  std::cout << box->getMaxPoint() << std::endl;
  box->clearStiff();
  box->clearMass();
  box->clearLumpedMass(); 
  box->clearElementsArray();     
  box->clearNodesArray();

  box->createNodesArray();

//  std::vector <NodeP> nodes = box->getNodes();
//  std::cout << "Size of nodes= " << nodes.size() << std::endl;

//  std::vector<NodeP> nodes = box->getNodes();
//  int size = nodes.size();
//  for (int ii = 0; ii < size; ii++)
//  {
//    NodeP cur_node = nodes[ii];
//    std::cout << "ID= " << cur_node->getID() << " Velocity= " << cur_node->velocity() << std::endl;
//  }

  box->createElementsArray();
//  std::cout << "Ok " << std::endl;
  if (!(d_mpm_box->getBoxFlag()))
  {
    findBoundaryNodesOfPoints();
  }

//  std::cout << "Ok " << std::endl;
//  box->makeGlobalMatrices();

}




double
MPMBoxPhases::findMax(const std::vector<double>&  vector)
{
  double max = -999999.0;
  int size = vector.size();
  for (int ii = 0; ii < size; ii++)
  {
    if (vector[ii] >= max) max = vector[ii];
  }
  return max;

}



double
MPMBoxPhases::findMin(const std::vector<double>&  vector)
{
  double min = 999999.0;
  int size = vector.size();
  for (int ii = 0; ii < size; ii++)
  {
    if (vector[ii] <= min) min = vector[ii];
  }
  return min;

}

//***********************************************************************************************************


void
MPMBoxPhases::mpmApplyPointsForce()
{
  std::vector<double> pointsXComponent;
  std::vector<double> pointsYComponent;
  std::vector<double> pointsZComponent;

  Vector3D zero(0.0);

  BoxSP box = d_mpm_box->getBox();
  Vector3D pointForce = box->getPointForce();
  Box::ForceSurface surface = box->getPointSurface();
  std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
  int pointSize = pointsArray.size();

  for (int k = 0; k < pointSize; k++)
  {
    MaterialPointP cur_point = pointsArray[k];
    pointsXComponent.emplace_back(cur_point->getPosOld().x());
    pointsYComponent.emplace_back(cur_point->getPosOld().y());   
    pointsZComponent.emplace_back(cur_point->getPosOld().z()); 

  }

  for (int j = 0; j < pointSize; j++)
  {
    MaterialPointP cur_point = pointsArray[j];
    Point3D position = cur_point->getPosOld();
    switch (surface)
    {
      case Box::ForceSurface::right:
        if (position.x() == findMax(pointsXComponent))
        {
           cur_point->setPointForce(pointForce);
           cur_point->setNoneZeroExtForce(true);
        }
      break;
      case Box::ForceSurface::left:
        if (position.x() == findMin(pointsXComponent))
        {
           cur_point->setPointForce(pointForce);
           cur_point->setNoneZeroExtForce(true);
        } 
      break;
      case Box::ForceSurface::up:
        if (position.y() == findMax(pointsYComponent))
        {
           cur_point->setPointForce(pointForce);
           cur_point->setNoneZeroExtForce(true);
        } 
      break;
      case Box::ForceSurface::down:
        if (position.y() == findMin(pointsYComponent))
        {
           cur_point->setPointForce(pointForce);
           cur_point->setNoneZeroExtForce(true);
        } 
      break;
      case Box::ForceSurface::front:
        if (position.z() == findMax(pointsZComponent))
        {
           cur_point->setPointForce(pointForce);
           cur_point->setNoneZeroExtForce(true);
        } 
      break;
      case Box::ForceSurface::back:
        if (position.z() == findMin(pointsZComponent))
        {
           cur_point->setPointForce(pointForce);
           cur_point->setNoneZeroExtForce(true);
        } 
      break;
      case Box::ForceSurface::none:
           cur_point->setPointForce(zero);
           cur_point->setNoneZeroExtForce(false);
      break;
      default:
           std::cout << "A box has just six surfaces!" << std::endl;
      break;
    }
   
  }

}


void
MPMBoxPhases::mpmApplyPointsForceGeometry()
{
  ComplicatedGeometrySP geometry = d_mpm_box->getGeometry();
  Vector3D pointForce = geometry->getPointForce();
  std::vector <MaterialPointP> points = geometry->getMaterialPointsArray();
  std::vector <int> pointsID = geometry->getPointForcePoints();
  int size = points.size();
  int sizeID = pointsID.size();
  for (int i = 0; i < size; i++)
  {
    MaterialPointP point = points[i];
    for (int j = 0; j < sizeID; j++)
    {
      int id = pointsID[j];
      if (point->getID() == id)
      {
        point->setPointForce(pointForce);
        point->setNoneZeroExtForce(true);
        continue;
      }
    }  
  }

}


void
MPMBoxPhases::mpmApplyBoundaryCondition()
{
  std::vector<double> pointsXComponent;
  std::vector<double> pointsYComponent;
  std::vector<double> pointsZComponent;

  Vector3D zero(0.0);

  BoxSP box = d_mpm_box->getBox();
  Vector3D pos = box->getBoundaryDisplacement();
  Point3D boundaryPosition(pos.x(), pos.y(), pos.z());
  Box::FixedDirichletBoundary surface = box->getFixedDirichletBoundary();
  std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
  int pointSize = pointsArray.size();

  for (int k = 0; k < pointSize; k++)
  {
    MaterialPointP cur_point = pointsArray[k];
    pointsXComponent.emplace_back(cur_point->getPosOld().x());
    pointsYComponent.emplace_back(cur_point->getPosOld().y());   
    pointsZComponent.emplace_back(cur_point->getPosOld().z()); 

  }

  for (int j = 0; j < pointSize; j++)
  {
    MaterialPointP cur_point = pointsArray[j];
    Point3D position = cur_point->getPosOld();
    switch (surface)
    {
      case Box::FixedDirichletBoundary::xMax:
        if (position.x() == findMax(pointsXComponent))
        {
//           cur_point->setPosOld(boundaryPosition);
//           cur_point->setPosNew(boundaryPosition);
           cur_point->setPosNew(position);
           cur_point->setHasBoundaryCondition(true);
        }
      break;
      case Box::FixedDirichletBoundary::xMin:
        if (position.x() == findMin(pointsXComponent))
        {
//           cur_point->setPosOld(boundaryPosition);
//           cur_point->setPosNew(boundaryPosition);
           cur_point->setPosNew(position);
           cur_point->setHasBoundaryCondition(true);
        }
      break;
      case Box::FixedDirichletBoundary::yMax:
        if (position.y() == findMax(pointsYComponent))
        {
//           cur_point->setPosOld(boundaryPosition);
//           cur_point->setPosNew(boundaryPosition);
           cur_point->setPosNew(position);
           cur_point->setHasBoundaryCondition(true);
        }
      break;
      case Box::FixedDirichletBoundary::yMin:
        if (position.y() == findMin(pointsYComponent))
        {
//           cur_point->setPosOld(boundaryPosition);
//           cur_point->setPosNew(boundaryPosition);
           cur_point->setPosNew(position);
           cur_point->setHasBoundaryCondition(true);
        }
      break;
      case Box::FixedDirichletBoundary::zMax:
        if (position.z() == findMax(pointsZComponent))
        {
//           cur_point->setPosOld(boundaryPosition);
//           cur_point->setPosNew(boundaryPosition);
           cur_point->setPosNew(position);
           cur_point->setHasBoundaryCondition(true);
        }
      break;
      case Box::FixedDirichletBoundary::zMin:
        if (position.z() == findMin(pointsZComponent))
        {
//           cur_point->setPosOld(boundaryPosition);
//           cur_point->setPosNew(boundaryPosition);
           cur_point->setPosNew(position);
           cur_point->setHasBoundaryCondition(true);
        }
      break;
      case Box::FixedDirichletBoundary::noFixedBoundary:
//           cur_point->setPointForce(zero);
           cur_point->setHasBoundaryCondition(false);
      break;
      default:
           std::cout << "A box has just six surfaces!" << std::endl;
      break;
    }
   
  }

}


void
MPMBoxPhases::mpmApplyBoundaryConditionGeometry()
{
  ComplicatedGeometrySP geometry = d_mpm_box->getGeometry();
  Vector3D boundaryDisplacement = geometry->getBoundaryDisplacement();
  std::vector <MaterialPointP> points = geometry->getMaterialPointsArray();
  std::vector <int> pointsID = geometry->getDirichletBoundaryNodes();
  int size = points.size();
  int sizeID = pointsID.size();
  for (int i = 0; i < size; i++)
  {
    MaterialPointP point = points[i];
    for (int j = 0; j < sizeID; j++)
    {
      int id = pointsID[j];
      if (point->getID() == id)
      {
        point->setPosNew(point->getPosOld() + boundaryDisplacement);
        point->setHasBoundaryCondition(true);
        continue;
      }
    }  
  }

}


void
MPMBoxPhases::findBoundaryNodesOfPoint(const MaterialPointP& point)
{
  double TOL = 0.000001;
  std::vector <NodeP> closestNodes;
  BoxSP box = getGridBox();
  std::vector<SymmetricMaterialTrilinearElementSP> elementsArray = box->getElements();
  double distance = 0.0;
  double min = d_mpm_box->getGeometry()->dist(elementsArray[0]->triElement()->getArray()[0], point);
  int size = elementsArray.size();
  for (int i = 0; i < size; i++)
  {
    SymmetricMaterialTrilinearElementSP cur_elem = elementsArray[i];
    // if (point->getHasBoundaryCondition())
    if (d_mpm_box->isPointInElement(point, cur_elem))
    {
      std::array <NodeP, 8>  nodes = cur_elem->triElement()->getArray();
//      double min = d_mpm_box->getGeometry()->dist(nodes[0], point);
      for (int k = 0; k < 8; k++)
      {
        if (min >= d_mpm_box->getGeometry()->dist(nodes[k], point));
        {
          min = d_mpm_box->getGeometry()->dist(nodes[k], point);
        }
      }   
    }
  }

  std::vector <NodeP> nodesArray = box->getNodes();
  int nodeSize = nodesArray.size();
  for (int j = 0; j < nodeSize; j++)
  {
    NodeP cur_node = nodesArray[j];
    distance = d_mpm_box->getGeometry()->dist(cur_node, point);
    if (std::abs(distance - min) < TOL)
    {
      closestNodes.emplace_back(cur_node);
      if (point->getHasBoundaryCondition())
      {
        cur_node->isBoundary(true);
      }
    }
  }
  
  point->setClosestNodes(closestNodes);  

}



void
MPMBoxPhases::findBoundaryNodesOfPoints()
{
  ComplicatedGeometrySP geometry = d_mpm_box->getGeometry();
  std::vector <MaterialPointP> points = geometry->getMaterialPointsArray();
  int size = points.size();
  for (int i = 0; i < size; i++)
  {
    MaterialPointP point = points[i];
    if (point->getHasBoundaryCondition())
    {
      findBoundaryNodesOfPoint(point);
    }
  
  }
    
}



bool
MPMBoxPhases::nodesBoundaryCondition(NodeP node)
{
  bool condition = false;
  if (d_mpm_box->getBoxFlag())
  { 
    BoxSP box = d_mpm_box->getBox();
    Box::FixedDirichletBoundary surface = box->getFixedDirichletBoundary();
      switch (surface)
      {
        case Box::FixedDirichletBoundary::xMax:
          if (node->x() == box->getMaxPoint().x())
          {
            condition = true;
          }
        break;
        case Box::FixedDirichletBoundary::xMin:
          if (node->x() == box->getMinPoint().x())
          {
            condition = true;
          }        
        break;
        case Box::FixedDirichletBoundary::yMax:
          if (node->y() == box->getMaxPoint().y())
          {
            condition = true;
          } 
        break;
        case Box::FixedDirichletBoundary::yMin:
          if (node->y() == box->getMinPoint().y())
          {
            condition = true;
          } 
        break;
        case Box::FixedDirichletBoundary::zMax:
          if (node->z() == box->getMaxPoint().z())
          {
            condition = true;
          } 
        break;
        case Box::FixedDirichletBoundary::zMin:
          if (node->z() == box->getMinPoint().z())
          {
            condition = true;
          } 
        break;
        case Box::FixedDirichletBoundary::noFixedBoundary:
           condition = false;
        break;
        default:
             std::cout << "A box has just six surfaces!" << std::endl;
        break;
      }

  }
  
  else
  {
    condition =  node->isBoundary();
  }

  return condition;
}



void
MPMBoxPhases::mpmSolver()
{
  if (d_mpm_box->getBoxFlag())
  {
    mpmApplyPointsForce();
    //////  mpmApplyBoundaryCondition();
  }
  else
  {
    redefineBackgroundGrid();

//    std::cout << "redefineBackgroundGrid();" << std::endl;

    mpmApplyPointsForceGeometry();
    mpmApplyBoundaryConditionGeometry();
//    findBoundaryNodesOfPoints();
  }


  int iter_num = 0;

  int status;
  int ret;
  std::string direcName = "outputFolder";
  makeOutputFolder(status, ret, direcName);

  char buffer[2000];
  char *str = getcwd(buffer, 2000);
  std::string currentDirec = std::string(buffer);
  std::cout << currentDirec << std::endl;
  
  while (iter_num*d_del_t < d_time)
  {  

    std::cout << "Num of iteration: " << iter_num << std::endl;

    if (iter_num % 50 == 1)
    {
//      std::string name = "initialPoints.txt";
      pointsFile(iter_num, currentDirec);

    }

//    nodesFile(iter_num, currentDirec);
//    elementsFile();

//    std::cout << "Before MPM " << std::endl;
    switch (d_solver_type)
    {
      case Lumped:
        initializationPhaseLumped();
//        std::cout << "After initializationPhaseLumped" << std::endl;
        lagrangianPhaseLumped();
//        std::cout << "After lagrangianPhaseLumped()" << std::endl;    
        convectivePhase();        
      break;

      case Consistent:
        initializationPhaseConsistent();
        lagrangianPhaseConsistent();    
        convectivePhase();        
      break;

      case LumpedUsingMomentum:
        initializationPhaseLumpedUsingMomentum();
        lagrangianPhaseLumpedUsingMomentum();    
        convectivePhase();        
      break;

      case ConsistentUsingMomentum:
        initializationPhaseConsistentUsingMomentum();
        lagrangianPhaseConsistentUsingMomentum();    
        convectivePhase();        
      break;

      default:
           std::cout << "Possible solvers: Lumped, Consistent, LumpedUsingMomentum, ConsistentUsingMomentum" << std::endl;
      break;
     }
 
//    if(d_solver_type == Lumped)
//    {
//        initializationPhaseLumped();
//        lagrangianPhaseLumped();    
//        convectivePhase();        
//    }

    updateBoundaryConditionAndPointsForce();

//      mpmPrintOld(iter_num);
      mpmPrintNew();
   
//    if((iter_num+1)*d_del_t >= d_time)
//    {
//      std::string name = "finalPoints.txt";
//      pointsFile(name);
//    }      
     iter_num += 1;    
  }
//  nodesFile();
//  elementsFile();

}



void
MPMBoxPhases::makeOutputFolder(int& status, int& ret, const std::string& direcName)
{
//  std::string direcName = "outputFolder";
  char buffer[2000];
  char *str = getcwd(buffer, 2000);
  std::string currentDirec = std::string(buffer);
//  std::cout << currentDirec << std::endl;
//  std::string direcLocation = "/home/hzar193/" + direcName;
  std::string direcLocation = currentDirec + "/" + direcName;
//  status = mkdir("/home/hzar193/outpoutFolder", S_IRWXU);
  status = mkdir(direcLocation.c_str(), S_IRWXU);
  ret = chdir(direcLocation.c_str());
  
}

void
MPMBoxPhases::updateBoundaryConditionAndPointsForce()
{

  Vector3D zero(0.0); 
  Vector3D pointForce(0.0);
  Vector3D pos(0.0);
 

  if (d_mpm_box->getBoxFlag())
  {
    BoxSP box = d_mpm_box->getBox();

    pointForce = box->getPointForce();

    pos = box->getBoundaryDisplacement();
  
//    std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
  }
  else
  {
    ComplicatedGeometrySP geometry = d_mpm_box->getGeometry();

    pointForce = geometry->getPointForce();

    pos = geometry->getBoundaryDisplacement();
  
  }


  Point3D boundaryDisplacement(pos.x(), pos.y(), pos.z());

  std::vector<MaterialPointP> pointsArray = getGeometryPoints();

  int number = 0;

  int pointSize = pointsArray.size();
  for (int k = 0; k < pointSize; k++)
  {
    MaterialPointP cur_point = pointsArray[k];
    if (cur_point->getNoneZeroExtForce())
    {
      number++;
      cur_point->setPointForce(pointForce);
      cur_point->setNoneZeroExtForce(true);
    }   

    if (cur_point->getHasBoundaryCondition())
    {
//    cur_point->setPosOld(boundaryPosition);
//    cur_point->setPosNew(boundaryPosition);
      cur_point->setPosNew(cur_point->getPosOld());
      cur_point->setHasBoundaryCondition(true);
    }

//    if (!(cur_point->getPointForce() == zero))
//    {
//      number ++;
//      std::cout << "ID= " << cur_point->getID() << " pointForce= " << cur_point->getPointForce() << std::endl;
//    } 
        

//    if (cur_point->getHasBoundaryCondition())
//    {
//      number ++;
//      std::cout << "ID= " << cur_point->getID() << " position= " << cur_point->getPosNew() << std::endl;
//    } 

  }
  std::cout << "The number of points exerted by point force= " << number << std::endl;

}


void
MPMBoxPhases::mpmPrintOld(int iter)
{
//  BoxSP box = d_mpm_box->getBox();
//  std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
  std::vector<MaterialPointP> pointsArray = getGeometryPoints();
  int pointSize = pointsArray.size();
  std::cout << "Old positions of points" << std::endl;
//  for (int k = 0; k < pointSize; k++)
//  {
//    MaterialPointP cur_point = pointsArray[k];
//  int i = 11;
/*  std::cout << "Point " << k << ": "; // << std::endl;
  std::cout << pointsArray[k]->getPosOld().x() << ",  "
            << pointsArray[k]->getPosOld().y() << ",  "
            << pointsArray[k]->getPosOld().z()  
            << "    disp= " << cur_point->getDisp()
                            << std::endl;*/


  int i = 250;
//  if (k == i)
//  {
//  std::cout << "Point " << i << ": "; // << std::endl;
  std::cout << /*iter << "  " <<*/ pointsArray[i]->getPosOld().x() << ", "
            << pointsArray[i]->getPosOld().y() << ", "
            << pointsArray[i]->getPosOld().z() 
            << " disp= " << pointsArray[i]->getDisp()
                            << std::endl;
//  }
//  }

}



void
MPMBoxPhases::mpmPrintNew()
{
//  BoxSP box = d_mpm_box->getBox();
//  std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
  std::vector<MaterialPointP> pointsArray = getGeometryPoints();
  int pointSize = pointsArray.size();
  std::cout << "New positions of points" << std::endl;
/*  for (int k = 0; k < pointSize; k++)
  {
    MaterialPointP cur_point = pointsArray[k];
//    int i = 950;
  std::cout << "Point " << k << ": "; // << std::endl;
  std::cout << pointsArray[k]->getPosNew().x() << ",  "
            << pointsArray[k]->getPosNew().y() << ",  "
            << pointsArray[k]->getPosNew().z() << std::endl; 

  }*/


  int i = 150;
//  std::cout << "Point " << i << ": "; // << std::endl;
  std::cout << pointsArray[i]->getPosNew().x() << ",  "
            << pointsArray[i]->getPosNew().y() << ",  "
            << pointsArray[i]->getPosNew().z()
            << " disp= " << pointsArray[i]->getDisp() << std::endl;

}



void
MPMBoxPhases::MUMPS(const std::vector<std::vector<double> >& mat, 
              const std::vector<double>& rhs,
              std::vector<double>& newDisp)
{
//  int size = d_box->getGlobalStiffness().size();
  int size = rhs.size();
//  std::cout << "Size= " << size << std::endl;
  VectorXd rhsVector(size);
  MatrixXd lhsMatrix(size,size);
  for (int i = 0; i < size; i++)
  {
     rhsVector[i] = rhs[i];
     for (int j = 0; j < size; j++)
     {
       lhsMatrix(i,j) = mat[i][j];
     }

  }

//    std::cout << lhsMatrix.determinant() << std::endl;	
    VectorXd x = lhsMatrix.colPivHouseholderQr().solve(rhsVector);
//  std::cout << "The solution is:\n" << x << std::endl;
    for (int k = 0; k < size; k++)
      newDisp.push_back(x[k]);  
}



  void
  MPMBoxPhases::initialMatrix(std::vector <std::vector <double> >& Matrix, int size)
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
  MPMBoxPhases::initialVector(std::vector <double>& Vector, int size)
  {
    for (int i = 0; i < size; i++)
    {
      Vector.push_back(0.0);
    }
  }

  void
  MPMBoxPhases::inverse(const std::vector <std::vector <double> >& mat,
                        std::vector <std::vector<double> >& invMatrix, 
                        int size)
  {
    MatrixXd lhsMatrix(size,size);
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        lhsMatrix(i,j) = mat[i][j];
      }
    }

    initialMatrix(invMatrix, size);

    MatrixXd inv(size, size);
    inv = lhsMatrix.inverse();
    for (int jj = 0; jj < size; jj++)
      for (int kk = 0; kk < size; kk++)
        invMatrix[jj][kk] = inv(jj,kk);
      
  }
    


  void
  MPMBoxPhases::nodesFile(int& iter_num, const std::string& currentFolder)
  {
//    int status1;
    int ret1;
//    std::string nodeFolderName = "nodesFolder";
//    if (iter_num == 0)
//    {
//      makeOutputFolder(status1, ret1, nodeFolderName);
//    }
//    else
//    {
//      std::string direcLocation = currentFolder + "/" + nodeFolderName;
//      ret1 = chdir(direcLocation.c_str());
//    }

    std::string fileName = "box_" + std::to_string(iter_num) + ".exnode";
    std::ofstream myfile(fileName);
    myfile << " Group name: box" // << "_" << std::to_string(iter_num) 
                                  << std::endl;
    myfile << " #Fields=1" << std::endl;
    myfile << " 1) coordinates, coordinate, rectangular cartesian, #Components=3" << std::endl;
    myfile << "   x.  Value index= 1, #Derivatives= 0" << std::endl;
    myfile << "   y.  Value index= 2, #Derivatives= 0" << std::endl;
    myfile << "   z.  Value index= 3, #Derivatives= 0" << std::endl;
//    BoxSP box = d_mpm_box->getBox();
    BoxSP box = getGridBox();
    std::vector<NodeP> nodesArray = box->getNodes();
    int size = nodesArray.size();
    for (int i = 0; i < size; i++)
    {
      NodeP cur_node = nodesArray[i];
      myfile << " Node:        " << cur_node->getID() << std::endl;
      myfile << "    " << cur_node->x() << std::endl;
      myfile << "    " << cur_node->y() << std::endl;
      myfile << "    " << cur_node->z() << std::endl;
    }

    myfile.close();
    ret1 = chdir(currentFolder.c_str());

  }
    
        
  void
  MPMBoxPhases::elementsFile()
  {
//    std::ofstream myfile("/people/hzar193/Desktop/elementsFile.txt");
    std::ofstream myfile("elementsFile.txt");
//    BoxSP box = d_mpm_box->getBox();
    BoxSP box = getGridBox();
    std::vector<SymmetricMaterialTrilinearElementSP> elementsArray = box->getElements();
    int size = elementsArray.size();
    for (int i = 0; i < size; i++)
    {
      SymmetricMaterialTrilinearElementSP cur_elem = elementsArray[i];
      myfile //<< "Element " 
             << cur_elem->id();
      std::array <NodeP, 8>  nodes = cur_elem->triElement()->getArray();
      myfile << " ";
      myfile << " " << nodes[0]->getID() << " "// << nodes[0]->getID() << " "  suitable for cmgui
                    << nodes[1]->getID() << " "// << nodes[1]->getID() << " "  suitable for cmgui 
                    << nodes[3]->getID() << " "// << nodes[3]->getID() << " "  suitable for cmgui
                    << nodes[2]->getID() << " "// << nodes[2]->getID() << " "  suitable for cmgui
                    << nodes[5]->getID() << " "// << nodes[5]->getID() << " "  suitable for cmgui
                    << nodes[6]->getID() << " "// << nodes[6]->getID() << " "  suitable for cmgui
                    << nodes[4]->getID() << " "// << nodes[4]->getID() << " "  suitable for cmgui
                    << nodes[7]->getID() << " ";//<< nodes[7]->getID() << " "  suitable for cmgui

//      for (int j = 0; j < 8; j++)
//      {
//        NodeP cur_node = nodes[j];
//        myfile << " " << cur_node->getID() << " ";
//      }
      myfile << std::endl;
    }
//      myfile << "Element " << cur_node->getID() << " " << cur_node->x() << " " << cur_node->y() 
//                        << " " << cur_node->z() << std::endl;
    

    myfile.close();

  }


  void
  MPMBoxPhases::pointsFile(int& iter_num, const std::string& currentFolder)
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
    myfile << " #Fields=1" << std::endl;
    myfile << " 1) coordinates, coordinate, rectangular cartesian, #Components=3" << std::endl;
    myfile << "   x.  Value index= 1, #Derivatives= 0" << std::endl;
    myfile << "   y.  Value index= 2, #Derivatives= 0" << std::endl;
    myfile << "   z.  Value index= 3, #Derivatives= 0" << std::endl;
//  BoxSP box = d_mpm_box->getBox();
//  std::vector<MaterialPointP> pointsArray = box->getMaterialPoints();
    std::vector<MaterialPointP> pointsArray = getGeometryPoints();
    int pointSize = pointsArray.size();
    std::cout << "size= " << pointSize << std::endl;
    for (int k = 0; k < pointSize; k++)
    {
      MaterialPointP cur_point = pointsArray[k];
      myfile << " Node:         " << cur_point->getID() << std::endl;
      myfile << "  " << pointsArray[k]->getPosNew().x() << "  "
             << pointsArray[k]->getPosNew().y() << "  "
             << pointsArray[k]->getPosNew().z()  
             << std::endl;
    }

    
    
    myfile.close();
    ret2 = chdir(currentFolder.c_str());

  }


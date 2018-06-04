#include <fstream>
#include <iostream>
#include <string>

#include <Solver.h>
#include <Box.h>
#include <eigen/Eigen/Dense>

using namespace FiniteElement;

Solver::Solver()
    :d_time(1.0), d_del_t(0.001), d_box(new Box()), d_geometry(new ComplicatedGeometry())
{
  setSolverType(Solver::solverType::Explicit);
  EIBoxSolvers(); 
}

Solver::Solver(BoxSP box)
    :d_time(1.0), d_del_t(0.001), d_box(box), d_geometry(new ComplicatedGeometry())
{
  setSolverType(Solver::solverType::Static);
  EIBoxSolvers(); 
}


Solver::Solver(solverType type, double time, double delT, BoxSP box)
    :d_time(time), d_del_t(delT), d_box(box), d_geometry(new ComplicatedGeometry())
{
//setBox(box);
  setSolverType(type);
  EIBoxSolvers();
}


/*Solver::Solver(ComplicatedGeometrySP geometry)
    :d_time(1.0), d_del_t(0.001), d_box(new Box()), d_geometry(geometry)
{
  setSolverType(Solver::solverType::Static);
  EIGeometrySolvers(); 
}*/


Solver::Solver(ComplicatedGeometrySP geometry)
    :d_time(1.0), d_del_t(0.001), d_box(new Box()), d_geometry(geometry)
{
  setSolverType(Solver::solverType::Static);
  EIGeometrySolvers();
  std::string name = "pointNew.exnode";
  d_geometry->nodesFile(name);  
}


Solver::Solver(solverType type, double time, double delT, ComplicatedGeometrySP geometry)
    :d_time(time), d_del_t(delT), d_box(new Box()), d_geometry(geometry)
{
  setSolverType(type);
  EIGeometrySolvers(); 
}



Solver::~Solver()
{
}


void
Solver::EIBoxSolvers()
{
  d_box->makeGlobalMatrices();
  std::vector<double> extForce = d_box->getExternalForce(); 
  std::vector<double> oldDisp = d_box->getOldGlobDisp();
  std::vector<double> disp = d_box->getGlobDisp();
  std::vector<double> newDisp; // = d_box->getNewGlobDisp();


  std::vector<std::vector<double> > stiffMatrix = d_box->getGlobalStiffness();
  std::vector<std::vector<double> > massLumpedMatrix = d_box->getGlobalLumpedMass();
  std::vector<std::vector<double> > massMatrix = d_box->getGlobalMass();
  int iter_num = 1;

  std::ofstream myfile ("nodeDisplacement.txt");



  if(d_solver_type == Static)
  {
    applyBoundaryCondition(stiffMatrix);
//    std::cout << "applyBoundaryCondition" << std::endl;
    d_box->setGlobalStiffness(stiffMatrix);
//    std::cout << "setGlobalStiffness" << std::endl;
    newDisplacementStatic(stiffMatrix, extForce, newDisp);
//    std::cout << "newDisplacementStatic" << std::endl;
    d_box->setOldGlobDisp(disp);
    d_box->setGlobDisp(newDisp);

    applyBoundaryCondition(stiffMatrix);
    d_box->setGlobalStiffness(stiffMatrix);

    oldDisp = d_box->getOldGlobDisp();
    disp = d_box->getGlobDisp();
    for (int ll = 0; ll < disp.size(); ll++)
    {
      if ((ll%3 == 0) /*&& (ll != 0)*/) std::cout << std::endl << ll/3+1 << "   ";
      std::cout << "(" << ll << ")= " << disp[ll] << "  ";
    }
    std::cout << std::endl;

    if (myfile.is_open())
    {
      myfile << "  " << iter_num*d_del_t << "  " << disp[3] << std::endl;     
    }
  }
  else
  {   
    while (iter_num*d_del_t < d_time)
    {
      newDisp.clear();

      if (d_solver_type == Explicit)
      {
        newDisplacementExplicit(stiffMatrix, massLumpedMatrix, extForce, oldDisp, disp, newDisp);
      }
      else
      {
        newDisplacementImplicit(stiffMatrix, massMatrix, extForce, oldDisp, disp, newDisp);
      }
      std::cout << "iter_num= " << iter_num;
//      for (int ll = 0; ll < newDisp.size(); ll++)
//      {
//        if ((ll%3 == 0) /*&& (ll != 0)*/) std::cout << std::endl << ll/3+1 << "   ";
//        std::cout << "(" << ll << ")= " << newDisp[ll] << "  ";
//      }
//      std::cout << std::endl;
////      d_box->getOldGlobDisp().clear();
////      d_box->initialVector(d_box->getOldGlobDisp(), size);
////      d_box->getGlobDisp().clear();
////      d_box->initialVector(d_box->getGlobDisp(), size);
      d_box->setOldGlobDisp(disp);
      d_box->setGlobDisp(newDisp);

      applyBoundaryCondition(stiffMatrix);
      d_box->setGlobalStiffness(stiffMatrix);

      oldDisp = d_box->getOldGlobDisp();
      disp = d_box->getGlobDisp();

      for (int ll = 0; ll < disp.size(); ll++)
      {
        if ((ll%3 == 0) /*&& (ll != 0)*/) std::cout << std::endl << ll/3+1 << "   ";
        std::cout << "(" << ll << ")= " << disp[ll] << "  ";
      }
      std::cout << std::endl;

      if (myfile.is_open())
      {
        myfile << "  " << iter_num*d_del_t << "  " << disp[3] << std::endl;     
      }
      else std::cout << "Unable to open file";

      iter_num++;
    }
 }
 myfile.close();
}


void
Solver::EIGeometrySolvers()
{
  d_geometry->makeGlobalMatrices();
//  std::cout << "Global matrices are made successfully." << std::endl;
  std::vector<double> extForce = d_geometry->getExternalForce(); 
  std::vector<double> oldDisp = d_geometry->getOldGlobDisp();
  std::vector<double> disp = d_geometry->getGlobDisp();
  std::vector<double> newDisp; // = d_geometry->getNewGlobDisp();


  std::vector<std::vector<double> > stiffMatrix = d_geometry->getGlobalStiffness();
  std::vector<std::vector<double> > massLumpedMatrix = d_geometry->getGlobalLumpedMass();
  std::vector<std::vector<double> > massMatrix = d_geometry->getGlobalMass();
  int iter_num = 1;

  char buffer[2000];
  char *str = getcwd(buffer, 2000);
  std::string currentDirec = std::string(buffer);

  std::ofstream myfile ("nodeDisplacement.txt");

  applyBoundaryConditionGeometry(stiffMatrix);

  if(d_solver_type == Static)
  {
//    applyBoundaryConditionGeometry(stiffMatrix);
//    std::cout << "applyBoundaryCondition" << std::endl;
    d_geometry->setGlobalStiffness(stiffMatrix);
//    std::cout << "setGlobalStiffness" << std::endl;
    newDisplacementStatic(stiffMatrix, extForce, newDisp);
//    std::cout << "newDisplacementStatic" << std::endl;
    d_geometry->setOldGlobDisp(disp);
    d_geometry->setGlobDisp(newDisp);

    applyBoundaryConditionGeometry(stiffMatrix);
    d_geometry->setGlobalStiffness(stiffMatrix);

    oldDisp = d_geometry->getOldGlobDisp();
    disp = d_geometry->getGlobDisp();
    for (int ll = 0; ll < disp.size(); ll++)
    {
      if ((ll%3 == 0) /*&& (ll != 0)*/) std::cout << std::endl << ll/3+1 << "   ";
      std::cout << "(" << ll << ")= " << disp[ll] << "  ";
//    }
//    std::cout << std::endl;

      if (myfile.is_open()) 
      {
        if ((ll%3 == 0) /*&& (ll != 0)*/) myfile << std::endl << ll/3+1 << "   ";
        myfile << "(" << ll << ")= " << disp[ll] << "  ";    
      }
    }
    myfile << std::endl;
    std::cout << std::endl;
  }
  else
  {   
    while (iter_num*d_del_t < d_time)
    {
      newDisp.clear();

      if (d_solver_type == Explicit)
      {
        newDisplacementExplicit(stiffMatrix, massLumpedMatrix, extForce, oldDisp, disp, newDisp);
      }
      else
      {
        newDisplacementImplicit(stiffMatrix, massMatrix, extForce, oldDisp, disp, newDisp);
      }
      std::cout << "iter_num= " << iter_num;
//      for (int ll = 0; ll < newDisp.size(); ll++)
//      {
//        if ((ll%3 == 0) /*&& (ll != 0)*/) std::cout << std::endl << ll/3+1 << "   ";
//        std::cout << "(" << ll << ")= " << newDisp[ll] << "  ";
//      }
//      std::cout << std::endl;
////      d_geometry->getOldGlobDisp().clear();
////      d_geometry->initialVector(d_geometry->getOldGlobDisp(), size);
////      d_geometry->getGlobDisp().clear();
////      d_geometry->initialVector(d_geometry->getGlobDisp(), size);
      d_geometry->setOldGlobDisp(disp);
      d_geometry->setGlobDisp(newDisp);

      applyBoundaryConditionGeometry(stiffMatrix);
      d_geometry->setGlobalStiffness(stiffMatrix);

      oldDisp.clear();
      disp.clear();
      oldDisp = d_geometry->getOldGlobDisp();
      disp = d_geometry->getGlobDisp();

//      for (int ll = 0; ll < disp.size(); ll++)
//      {
//        if ((ll%3 == 0) /*&& (ll != 0)*/) std::cout << std::endl << ll/3+1 << "   ";
//        std::cout << "(" << ll << ")= " << disp[ll] << "  ";
//      }

      std::cout << " (64) " << disp[189] << " " << disp[190] << " " << disp[191];
      std::cout << std::endl;
      
      updateNodesPos(disp);    


      if (iter_num % 5000 == 1)
      {
  //      std::string name = "initialPoints.txt";
        nodesFile(iter_num, currentDirec);

      }


      if (myfile.is_open())
      {
//        myfile << "  " << iter_num*d_del_t << "  " << disp[3] << std::endl;
        myfile << "  " << iter_num*d_del_t << "  " << disp[189] << "  "
                       << disp[190] << " " << disp[191] << std::endl;        
      }
      else std::cout << "Unable to open file";

      iter_num++;
    }
 }
 myfile.close();
}




void
Solver::newDisplacementExplicit(const std::vector<std::vector<double> >& stiffMatrix,
                                const std::vector<std::vector<double> >& massLumpedMatrix,
                                const std::vector<double>& extForce, const std::vector<double>& oldDisp,
                                const std::vector<double>& disp, std::vector<double>& newDisp)
{
//  int size = d_box->getGlobalStiffness().size();
  int size = disp.size();
  double sum = 0.0;
  for (int i = 0; i < size; i++)
  {
    sum = 0.0;
    for (int j = 0; j < size; j++)
    {
      sum += stiffMatrix[i][j]*disp[j];
    }
    newDisp.push_back((d_del_t*d_del_t)*(1/massLumpedMatrix[i][i])*(extForce[i]-sum)+2*disp[i]-oldDisp[i]);
  }
}



void 
Solver::implicitSolver()
{
}



std::vector<double> 
Solver::massDotVec(const std::vector<std::vector<double> >& matrix, 
                   const std::vector<double>& vec)
{
  double sum = 0.0;
  std::vector<double> result;
  int size = matrix.size();
  for (int i = 0; i < size; i++)
  {
    sum = 0.0;
    for (int j = 0; j < size; j++) 
    {
      sum += matrix[i][j]*vec[j];
    }
    result.push_back(sum);
  }
return result;
}



void
Solver::newDisplacementImplicit(const std::vector<std::vector<double> >& stiffMatrix,
                                const std::vector<std::vector<double> >& massMatrix,
                                const std::vector<double>& extForce, const std::vector<double>& oldDisp,
                                const std::vector<double>& disp, std::vector<double>& newDisp)
{
//  int size = d_box->getGlobalStiffness().size();
  int size = disp.size();
//  std::cout << "Size= " << size << std::endl;
  std::vector<std::vector<double> > mat;
  std::vector<double> rhs;
  for (int ii = 0; ii < size; ii++)
  {
//    rhs[ii] = d_del_t*d_del_t*extForce[ii] + 2*disp[ii] - oldDisp[ii];
    double r = d_del_t*d_del_t*extForce[ii] + 2*massDotVec(massMatrix,disp)[ii] - massDotVec(massMatrix,oldDisp)[ii];
    rhs.push_back(r);
//    std::cout << rhs[ii] << std::endl;
    std::vector<double> vec;
    for (int jj = 0; jj < size; jj++)
    {
//        mat[ii][jj] = massMatrix[ii][jj] + d_del_t*d_del_t*stiffMatrix[ii][jj];
      double m = massMatrix[ii][jj] + d_del_t*d_del_t*stiffMatrix[ii][jj];
      vec.push_back(m);
    }
    mat.push_back(vec);   
  }
  MUMPS(mat, rhs, newDisp); 
}


void
Solver::newDisplacementStatic(const std::vector<std::vector<double> >& stiffMatrix,
                                const std::vector<double>& extForce,
                                std::vector<double>& newDisp)
{
  MUMPS(stiffMatrix, extForce, newDisp);
//  std::cout << "MUMPS" << std::endl;
}
  

void
Solver::MUMPS(const std::vector<std::vector<double> >& mat, 
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

//  std::cout << lhsMatrix.determinant() << std::endl;	
  VectorXd x = lhsMatrix.colPivHouseholderQr().solve(rhsVector);
//  std::cout << "The solution is:\n" << x << std::endl;
  for (int k = 0; k < size; k++)
    newDisp.push_back(x[k]);  
}


void 
Solver::applyBoundaryCondition(std::vector<std::vector<double> >& stiffMatrix)
  {
    std::vector<double> displacement = d_box->getGlobDisp();
    std::vector<NodeP> nodesArray = d_box->getNodes();
    
    Box::FixedDirichletBoundary fixedBoundary = d_box->getFixedDirichletBoundary();
    Vector3D boundaryPosition = d_box->getBoundaryDisplacement();

    Vector3D minPoint = d_box->getMinPoint();
    Vector3D maxPoint = d_box->getMaxPoint();

    int size = nodesArray.size();
    for (int ii = 0; ii < size; ii++)
    {
      NodeP cur_node = nodesArray[ii];
      if (condition(cur_node, fixedBoundary, minPoint, maxPoint))
        {
          displacement[3*ii] = boundaryPosition.x();
          displacement[3*ii+1] = boundaryPosition.y();
          displacement[3*ii+2] = boundaryPosition.z();
          for (int jj = 0; jj < stiffMatrix.size(); jj++)
          {
            stiffMatrix[3*ii][jj] = 0.0;
            stiffMatrix[3*ii+1][jj] = 0.0;
            stiffMatrix[3*ii+2][jj] = 0.0;
          }
          stiffMatrix[3*ii][3*ii] = 1.0;
          stiffMatrix[3*ii+1][3*ii+1] = 1.0;
          stiffMatrix[3*ii+2][3*ii+2] = 1.0;
        }
      double newX = cur_node->x() + displacement[3*ii];
      double newY = cur_node->y() + displacement[3*ii+1];
      double newZ = cur_node->z() + displacement[3*ii+2];

      Point3D pos(newX, newY, newZ);
      cur_node->position(pos);
//      cur_node->position().x(newX);
//      cur_node->position().y(cur_node->y() + displacement[3*ii+1]);
//      cur_node->position().z(cur_node->z() + displacement[3*ii+2]);     
    }
    
    d_box->setGlobDisp(displacement);

  }


void 
Solver::applyBoundaryConditionGeometry(std::vector<std::vector<double> >& stiffMatrix)
  {
    std::vector<double> displacement = d_geometry->getGlobDisp();
//    std::vector <NodeP> nodes;
//    nodes.clear();
    std::vector<NodeP> nodesArray = d_geometry->getNodes();
    
//    Box::FixedDirichletBoundary fixedBoundary = d_box->getFixedDirichletBoundary();
    Vector3D boundaryPosition = d_geometry->getBoundaryDisplacement();

//    Vector3D minPoint = d_box->getMinPoint();
//    Vector3D maxPoint = d_box->getMaxPoint();

    int size = nodesArray.size();
    for (int ii = 0; ii < size; ii++)
    {
      NodeP cur_node = nodesArray[ii];
//      if (condition(cur_node, fixedBoundary, minPoint, maxPoint))
      if (d_geometry->isNodeNumberInArray(cur_node, d_geometry->getDirichletBoundaryNodes()))
        {
//          std::cout << cur_node->getID() << " ";
          displacement[3*ii] = boundaryPosition.x();
          displacement[3*ii+1] = boundaryPosition.y();
          displacement[3*ii+2] = boundaryPosition.z();
          for (int jj = 0; jj < stiffMatrix.size(); jj++)
          {
            stiffMatrix[3*ii][jj] = 0.0;
            stiffMatrix[3*ii+1][jj] = 0.0;
            stiffMatrix[3*ii+2][jj] = 0.0;
          }
          stiffMatrix[3*ii][3*ii] = 1.0;
          stiffMatrix[3*ii+1][3*ii+1] = 1.0;
          stiffMatrix[3*ii+2][3*ii+2] = 1.0;
        }
////      double newX = cur_node->x() + displacement[3*ii];
////      double newY = cur_node->y() + displacement[3*ii+1];
////      double newZ = cur_node->z() + displacement[3*ii+2];

////      Point3D pos(newX, newY, newZ);
////      cur_node->position(pos);
////      nodes.emplace_back(cur_node);
      
//      cur_node->position().x(newX);
//      cur_node->position().y(cur_node->y() + displacement[3*ii+1]);
//      cur_node->position().z(cur_node->z() + displacement[3*ii+2]);     
    }
    
    d_geometry->setGlobDisp(displacement);
    d_geometry->setNodes(nodesArray);

  }



bool
Solver::condition(const NodeP& node, const Box::FixedDirichletBoundary& boundary,
                  const Vector3D& min, const Vector3D& max)
  {
//    Box::FixedDirichletBoundary fixedBoundary = d_box->getFixedDirichletBoundary();
//    double boundaryPosition = d_box->getBoundaryPosition();
    bool cond;
    switch (boundary)
    {
      case Box::FixedDirichletBoundary::xMax:
        cond = (node->x() == max.x());
      break;
      case Box::FixedDirichletBoundary::xMin:
        cond = (node->x() == min.x());
      break;
      case Box::FixedDirichletBoundary::yMax:
        cond = (node->y() == max.y());
      break;
      case Box::FixedDirichletBoundary::yMin:
        cond = (node->y() == min.y());
      break;
      case Box::FixedDirichletBoundary::zMax:
        cond = (node->z() == max.z());
      break;
      case Box::FixedDirichletBoundary::zMin:
        cond = (node->z() == min.z());
      break;
      case Box::FixedDirichletBoundary::noFixedBoundary:
        cond = false;
      break;
      default:
        std::cout << "A box has just six surfaces!" << std::endl;
        cond = false;
      break;
    }
     return cond;
  }


  void
  Solver::nodesFile(const int& iter_num, const std::string& currentFolder)
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
//    BoxSP box = getGridBox();
    std::vector<NodeP> nodesArray = d_geometry->getNodes();
    int size = nodesArray.size();
    for (int i = 0; i < size; i++)
    {
      NodeP cur_node = nodesArray[i];
      myfile << " Node:        " << cur_node->getID() << std::endl;
//      myfile << "    " << cur_node->x() << std::endl;
//      myfile << "    " << cur_node->y() << std::endl;
//      myfile << "    " << cur_node->z() << std::endl;

      myfile << "    " << cur_node->positionNew().x() << std::endl;
      myfile << "    " << cur_node->positionNew().y() << std::endl;
      myfile << "    " << cur_node->positionNew().z() << std::endl;
    }

    myfile.close();
    ret1 = chdir(currentFolder.c_str());

  }


void
Solver::updateNodesPos(const std::vector<double>& disp)
{
  std::vector<NodeP> nodesArray = d_geometry->getNodes();
  int size = nodesArray.size();
  int index = 0;
  Point3D pos(0.0, 0.0, 0.0);
  for (int i = 0; i < size; i++)
  {  
    NodeP cur_node = nodesArray[i];
    index = cur_node->getID() - 1;
    pos.x(cur_node->x() + disp[index]);
    pos.y(cur_node->y() + disp[index + 1]);    
    pos.z(cur_node->z() + disp[index + 2]);    
    cur_node->positionNew(pos);
  }
}




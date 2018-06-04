#include <iostream>
#include <array>
#include <vector>
#include <Tensor.h>
#include <GeometryMath/Vector3D.h>
#include <GeometryMath/Point3D.h>
#include <GeometryMath/Matrix3D.h>
#include <GeometryMath/Types.h>
#include <NodeP.h>
#include <Node.h>
#include <TrilinearVolumeElement.h>
#include <TrilinearVolumeElementSP.h>
#include <SymmetricMaterial.h>
#include <SymmetricMaterialSP.h>
#include <SymmetricMaterialTrilinearElement.h>
#include <SymmetricMaterialTrilinearElementSP.h>
#include <Body.h>
#include <Box.h>
#include <BoxSP.h>
#include <Solver.h>

#include <PeriBoxP.h>
#include <PeriBox.h>
#include <PeriPMBMaterialP.h>
#include <PeriPMBMaterial.h>
#include <PeridynamicsP.h>
#include <Peridynamics.h>
#include <PeriCrackP.h>
#include <PeriCrack.h>


#include <ComplicatedGeometry.h>
#include <ComplicatedGeometrySP.h>
#include <string>

//#include <armadillo>

#include <ctime>
#include <random>
#include <math.h>


#include <MaterialPointP.h>  // MPM
#include <MaterialPoint.h>   // MPM

//#define  ARMA_DONT_USE_WRAPPER
//#define  ARMA_USE_LAPACK
//#define ARMA_USE_BLAS
//#define ARMA_USE_ATLAS
//#include <armadillo-4.500.0/include/armadillo>
//#include <armadillo-4.450.4/include/armadillo>

//#include <clapack.h>
//#include<armadillo-4.450.4/include/armadillo>
//#include <home/hzar193/libraries/armadillo-4.450.4/include/armadillo>
//#include <lapacke.h>
//#include <liblas/capi/liblas.h>

//typedef std::array <double, 3> Array3;
//typedef std::array <Array3, 3> Matrix3;


//using namespace arma;
using namespace FiniteElement;

int main () {


//Vector3D min_point(0.0, 0.0, 0.0);
//Vector3D max_point(1.5, 1.5, 1.5);
//Vector3D max_point(0.01, 0.007, 0.007);


//int nx = 1;
//int ny = 1;
//int nz = 1;

//double delT = 0.000025;
//double delT = 0.0001;
//double delT = 0.000048;
//double delT = 5*pow(10, -8);;


//double totalT = 2.0;
//double totalT = 0.02;
//double totalT = 0.1;
//double totalT = 5*pow(10, -7);
//double totalT = 0.0001;

double density = 2440.0;
//double density = 2700.0;

double poission = 0.25;
 
//double young = 72000000000;
double young = 720000;
//double young = 70*pow(10, 9); // E = 70 GPa

//double bodyForce = 0.0;


//double surfaceForce = -30000;
//double surfaceForce = -35*pow(10, 6);
//double surfaceForce = 0.0;
//double surfaceForce = -50000;
//double surfaceForce = -3000;


Vector3D boundaryDisplacement(0.0);


//Vector3D unitVector(-0.1581, 0.1645, -0.9736);
Vector3D unitVector(0.180291402, -0.123449913, 0.975838965);  // unitVector of (NodeID 16.pos - NodeID 291.pos)
double magnitude = 1000;

//int numPointsX = 10;
//int numPointsY = 10;
//int numPointsZ = 10;

//int numPointsX = 20;
//int numPointsY = 20;
//int numPointsZ = 20;

//int numPointsX = 21;
//int numPointsY = 15;
//int numPointsZ = 15;

//int numPointsX = 11;
//int numPointsY = 8;
//int numPointsZ = 8;



//int numPointsX = 28;
//int numPointsY = 28;
//int numPointsZ = 28;

 
Vector3D tractionForce = unitVector*magnitude;
//Box::ForceSurface surface = Box::ForceSurface::right;

//PeriBox::ForceSurface surface = PeriBox::ForceSurface::right;
//PeriBox::ForceSurface surface = PeriBox::ForceSurface::front;
//PeriBox::ForceSurface surface = PeriBox::ForceSurface::down;
//PeriBox::ForceSurface surface = PeriBox::ForceSurface::upDown;
//PeriBox::ForceSurface surface = PeriBox::ForceSurface::none;  
  
  

//Box::FixedDirichletBoundary boundary = Box::FixedDirichletBoundary::noFixedBoundary;
//Box::FixedDirichletBoundary boundary = Box::FixedDirichletBoundary::xMin;


//BoxSP test_box(new Box(min_point, max_point, nx, ny, nz, density, poission, young,
//               bodyForce, surfaceForce, surface, boundary, boundaryDisplacement));

//Solver::solverType type1 = Solver::solverType::Explicit;
//Solver solver(type1, 100, 0.0002, test_box);

//Solver::solverType type2 = Solver::solverType::Implicit;
//Solver solver(type2, 150, 0.0002, test_box);

//Solver solver(test_box);


std::string nodeFile = "BoneDataNodes.txt";
std::string elementFile = "BoneDataElements.txt";

std::string tractionNodesFile = "downRadius.exnode";
std::string dirichletBoundaryNodesFile = "upRadius.exnode";


ComplicatedGeometrySP test_bone(new ComplicatedGeometry(nodeFile, elementFile, 
                                tractionNodesFile, dirichletBoundaryNodesFile,
                                density, young, poission, tractionForce, boundaryDisplacement));

//Solver solver(test_bone);
Solver::solverType type1 = Solver::solverType::Explicit;
Solver solver(type1, 20, 0.000005, test_bone);
//Solver solver(type1, 0.000004, 0.000002, test_bone);

//Solver::solverType type2 = Solver::solverType::Implicit;
//Solver solver(type2, 100, 0.002, test_bone);

// ComplicatedGeometry test_bone(nodeFile, elementFile, density, young, poission, 5, 5, 5);


//ComplicatedGeometrySP test_bone(new ComplicatedGeometry(nodeFile, elementFile, 
//                                tractionNodesFile, dirichletBoundaryNodesFile,
//                                density, young, poission, tractionForce, boundaryDisplacement, 2, 2, 2));



/********************************* Peridynamics **************************/

//Vector3D min_point(0.0, 0.0, 0.0);
//Vector3D max_point(1.5, 1.5, 1.5);
//Vector3D max_point(0.01, 0.007, 0.007);


//int numPointsX = 21;
//int numPointsY = 15;
//int numPointsZ = 15;

//int numPointsX = 11;
//int numPointsY = 8;
//int numPointsZ = 8;


//double delT = 0.000025;
//double delT = 0.0001;
//double delT = 0.000048;
//double delT = 5*pow(10, -8);

//double totalT = 2.0;
//double totalT = 0.02;
//double totalT = 0.1;
//double totalT = 0.0001;
//double totalT = 15*pow(10, -8);
//double totalT = 1.8*pow(10, -3);
//double totalT = 1.8*pow(10, -2);


//double poission = 0.25;

//double density = 2440.0;
//double density = 2700.0;
 
//double young = 72000000000;
//double young = 70*pow(10, 9); // E = 70 GPa

//double fractureEnergy = 135;


//double horizon = 0.22;
//double horizon = 0.3;
//double horizon = 0.002;
//double horizon = 0.001;



//double bodyForce = 0.0;


//Box::FixedDirichletBoundary boundary = Box::FixedDirichletBoundary::noFixedBoundary;
////Box::FixedDirichletBoundary boundary = Box::FixedDirichletBoundary::xMin;

//Vector3D boundaryDisplacement(0.0);


////PeriBox::ForceSurface surface = PeriBox::ForceSurface::upDown;
//PeriBox::ForceSurface surface = PeriBox::ForceSurface::none;

////double surfaceForce = -30000;
////double surfaceForce = -35*pow(10, 6);
//double surfaceForce = 0.0;
////double surfaceForce = -50000;
////double surfaceForce = -3000;


//PeriBox::ForceSurface velocitySurface = PeriBox::ForceSurface::upDownX;
////PeriBox::ForceSurface velocitySurface = PeriBox::ForceSurface::none;
//Vector3D velocity(0.01, 0.0, 0.0);


////Vector3D min(0.0, 0.7805, 0.0);
////Vector3D Max(1.5, 0.7905, 0.3);
////double PreCrackLengthRatio = 0.5;


//PeriPMBMaterialP material(new PeriPMBMaterial(density, horizon, fractureEnergy, young));
//PeriBoxP peri_box(new PeriBox(min_point, max_point, numPointsX, numPointsY, numPointsZ,
//                              density, poission, young, bodyForce, surfaceForce, surface,
//                              boundary, boundaryDisplacement, horizon, material,
//                              velocitySurface, velocity));
////PeriCrackP peri_crack(new PeriCrack(min, Max));

//PeridynamicsP test_prei(new Peridynamics(totalT, delT, peri_box));                  
////PeridynamicsP test_prei(new Peridynamics(totalT, delT, peri_box, peri_crack));




/************************** Material Point Method *************************/

/*BoxSP test_box(new Box(min_point, max_point, nx, ny, nz, 10, 10, 10, density, poission, young,
               bodyForce, surfaceForce, surface, boundary, boundaryPosition));

std::vector<MaterialPointP> pointsArray = test_box->getMaterialPoints();
int size = pointsArray.size();
std::cout  << size << std::endl;
for (int kk = 0; kk < size; kk++)
{
  MaterialPointP point = pointsArray[kk];
  std::cout << "Point(" << point->getID() << "): " << std::endl 
            << "x= " << point->getPosOld().x() << "  y= " << point->getPosOld().y()
            << "  z= " << point->getPosOld().z() << std::endl <<std::endl;

}*/

/****************************************************************************/






/************Tensor tt(1);
Tensor ss(2);
//   //std::array <std::array<double,3>, 3> mat;
Matrix3D mat;
for (int ii=0; ii<3; ii++){
  for(int jj=0; jj<3; jj++){
    mat(ii,jj)= 4;
//  //  mat.set(ii, jj, 4);
  }
}
std::cout << "Tensor[1][0][2][2]= " << tt(1,0,2,2) << std::endl;
std::cout << "Contract= " << tt.Contract(ss) << std::endl;
std::cout << "doubleContract= " << std::endl;
for (int ii=0; ii<3; ii++){
  for(int jj=0; jj<3; jj++){
  //  std::cout << ss.doubleContract(mat).get(ii,jj) << " ";
    std::cout << ss.doubleContract(mat)(ii,jj) << " ";
  }
  std::cout << std::endl;
}

NodeP node1(new Node(1, 0.0, 0.0, 0.0, false));
NodeP node2(new Node(2, 10.0, 0.0, 0.0, false));
NodeP node3(new Node(3, 10.0, 10.0, 0.0, false));
NodeP node4(new Node(4, 0.0, 10.0, 0.0, false));
NodeP node5(new Node(5, 0.0, 10.0, 10.0, false));
NodeP node6(new Node(6, 0.0, 0.0, 10.0, false));
NodeP node7(new Node(7, 10.0, 0.0, 10.0, false));
NodeP node8(new Node(8, 10.0, 10.0, 10.0, false));

node1->densityNode(1000);
node2->densityNode(1200);
node3->densityNode(1400);
node4->densityNode(1600);
node5->densityNode(1600);
node6->densityNode(1000);
node7->densityNode(1200);
node1->densityNode(1400);

Vector3D v1(100);
node1->surfaceForce(v1);

// //Vector3D v2(120);
// //node2->surfaceForce(v2);

// //Vector3D v3(140);
// //node3->surfaceForce(v3);

Vector3D v4(160);
node4->surfaceForce(v4);

Vector3D v5(160);
node5->surfaceForce(v5);

Vector3D v6(100);
node6->surfaceForce(v6);

Vector3D v7(120);
node7->surfaceForce(v7);

Vector3D v8(140);
node8->surfaceForce(v8);

Vector3D body_v1(10);
node1->bodyForce(body_v1);

Vector3D body_v2(12);
node2->bodyForce(body_v2);

Vector3D body_v3(14);
node3->bodyForce(body_v3);

Vector3D body_v4(16);
node4->bodyForce(body_v4);

Vector3D body_v5(16);
node5->bodyForce(body_v5);

Vector3D body_v6(10);
node6->bodyForce(body_v6);

Vector3D body_v7(12);
node7->bodyForce(body_v7);

Vector3D body_v8(14);
node8->bodyForce(body_v8);




TrilinearVolumeElementSP element(new TrilinearVolumeElement(node1, node2, node3,
                                                            node4, node5, node6,
                                                            node7, node8));

std::cout << "testFunction(7, 0.25, 0.5, 0.75)= " << element->testFunction(7, 0.25, 0.5, 0.75) << std::endl;
std::cout << "derivativeTestFunctionToXi(3, 2, 0.25, 0.5, 0.75)= " << element->derivativeTestFunctionToXi(3, 2, 0.25, 0.5, 0.75) << std::endl;
   
Matrix3D jacob = element->jacobianMatrix(0.25, 0.5, 0.75);
std::cout << "jacobianMatrix(0.25, 0.5, 0.75)= " << std::endl;
for (int ii=0; ii<3; ii++){
  for(int jj=0; jj<3; jj++){
//  //  std::cout << jacob.get(ii,jj) << " ";
    std::cout << jacob(ii,jj) << " ";
  }
  std::cout << std::endl;
}

std::cout << "derivativeTestFunctionToPosition(3, 2, 0.25, 0.5, 0.75)= " << element->derivativeTestFunctionToPosition(3, 2, 0.25, 0.5, 0.75) << std::endl;


// //etricMaterialSP material(new SymmetricMaterial(1, 1000, 72000000000, 0.25));
// //  for (int i=0;i<3;i++){
// //    for (int j=0;j<3;j++){
// //      for (int k=0;k<3;k++){
// //        for (int l=0;l<3;l++){
// //          std::cout << "C[" << i << "][" << j << "][" << k << "][" << l << "]= "
// //                    << material->tensor()(i, j, k, l) << "  ";
// //        }
// //      std::cout << std::endl;
// //      }
// //    }
// //  }

SymmetricMaterialTrilinearElementSP materialTrilinear(new SymmetricMaterialTrilinearElement(1, 1000, 72000000000,
                                                                                           0.25, element));
TrilinearVolumeElementSP trilinearElement= materialTrilinear->triElement();
std::cout << "symTestFunction(7, 0.25, 0.5, 0.75)= " << trilinearElement->symTestFunction(7, 0.25, 0.5, 0.75) << std::endl;
std::cout << "symDerivativeTestFunctionToXi(3, 2, 0.25, 0.5, 0.75)= " << 
              trilinearElement->symDerivativeTestFunctionToXi(3, 2, 0.25, 0.5, 0.75) << std::endl;
   
Matrix3D jacobian = trilinearElement->symJacobianMatrix(0.25, 0.5, 0.75);
std::cout << "jacobianMatrix(0.25, 0.5, 0.75)= " << std::endl;
for (int ii=0; ii<3; ii++){
  for(int jj=0; jj<3; jj++){
//  //  std::cout << jacob.get(ii,jj) << " ";
    std::cout << jacobian(ii,jj) << " ";
  }
  std::cout << std::endl;
}

std::cout << "symDerivativeTestFunctionToPosition(3, 2, 0.25, 0.5, 0.75)= " << 
             trilinearElement->symDerivativeTestFunctionToPosition(3, 2, 0.25, 0.5, 0.75) << std::endl;

for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          std::cout << "C[" << i << "][" << j << "][" << k << "][" << l << "]= "
                    << materialTrilinear->tensor()(i, j, k, l) << "  ";
        }
      std::cout << std::endl;
      }
    }
}

unsigned int kk = 6;
materialTrilinear->stiffnessMatrix(kk);
std::cout << std::endl << "Stiffness Matrix: " << std::endl;
for (unsigned int ii = 0; ii < 24; ii++) {
  for (unsigned int jj = 0; jj < 24; jj++) {
    std::cout << "SM(" << ii << ", " << jj << ")= " << materialTrilinear->stiffMat()(ii, jj) << "  ";
    if (jj%3 == 2) std::cout << std::endl;
  }
}
std::cout << std::endl;

// //Matrix24 invMat;
// //arma::inv(invMat, materialTrilinear->stiffMat());
// //invMat = arma::inv(materialTrilinear->stiffMat());
// //Matrix24 invMat = arma::auxlib::inv(materialTrilinear->stiffMat());
// //for (unsigned int ii = 0; ii < 24; ii++) {
// //  for (unsigned int jj = 0; jj < 24; jj++) {
// //    std::cout << "INVSM(" << ii << ", " << jj << ")= " << invMat(ii, jj) << "  ";
// //    if (jj%3 == 2) std::cout << std::endl;
// //  }
// //}
// //std::cout << std::endl;



materialTrilinear->massMatrix(kk);
std::cout << std::endl << "Mass Matrix: " << std::endl;
for (unsigned int ii = 0; ii < 24; ii++) {
  for (unsigned int jj = 0; jj < 24; jj++) {
    std::cout << "MM(" << ii << ", " << jj << ")= " << materialTrilinear->massMat()(ii, jj) << "  ";
    if (jj%3 == 2) std::cout << std::endl;
  }
}
std::cout << std::endl;

std::cout << "Norm Vector of plane xi3 = 1: "  << trilinearElement->symNormVecXi3Max(-0.25, 0.75) << std::endl;
std::cout << "Norm Vector of plane xi3 = -1: " << trilinearElement->symNormVecXi3Min(-0.25, 0.75) << std::endl;
std::cout << "Norm Vector of plane xi2 = 1: "  << trilinearElement->symNormVecXi2Max(-0.25, 0.75) << std::endl;
std::cout << "Norm Vector of plane xi2 = -1: " << trilinearElement->symNormVecXi2Min(-0.25, 0.75) << std::endl;
std::cout << "Norm Vector of plane xi1 = 1: "  << trilinearElement->symNormVecXi1Max(-0.25, 0.75) << std::endl;
std::cout << "Norm Vector of plane xi1 = -1: " << trilinearElement->symNormVecXi1Min(-0.25, 0.75) << std::endl;




materialTrilinear->surfaceForceVector(kk);
std::cout << std::endl << "Surface Force Vector: " << std::endl;
for (unsigned int ii = 0; ii < 24; ii++) {
//  //for (unsigned int jj = 0; jj < 24; jj++) {
    std::cout << "SFV(" << ii <<  ")= " << materialTrilinear->surfaceVec()(ii) << "  ";
    if (ii%3 == 2) std::cout << std::endl;
// //  }
}
std::cout << std::endl;


materialTrilinear->bodyForceVector(kk);
std::cout << std::endl << "Body Force Vector: " << std::endl;
for (unsigned int ii = 0; ii < 24; ii++) {
    std::cout << "BFV(" << ii <<  ")= " << materialTrilinear->bodyVec()(ii) << "  ";
    if (ii%3 == 2) std::cout << std::endl;
//  }
}
std::cout << std::endl; 


Body body;
std::cout << body.globalIndex(materialTrilinear, 23) << std::endl;
// //SymmetricMaterialTrilinearElement a;
// //std::cout << "twoPointsIntegral of examFunc is: " << 
// //           a.twoPointsIntegral(&SymmetricMaterialTrilinearElement::examFunc, a) 
// //          << std::endl;

// //std::cout << "threePointsIntegral of examFunc is: " << 
// //           a.threePointsIntegral(&SymmetricMaterialTrilinearElement::examFunc, a) 
// //          << std::endl;

// //a.printExamFunc();




// //  arma::mat A = "-1 4 8 2;5 6 7 8;-10 -9 11 12; 14 18 15 20";
// //  int numberRows = A.n_rows;
// // mat A = randu<mat>(5,5);
// //  std::cout << "detA= " << det(A) << std::endl;
// // mat B = A.i();
// //  arma::mat B = arma::inv(A);
// //  for (int i=0; i<numberRows; i++) {
// //    for (int j=0; j<numberRows; j++) {
// //      std::cout << "A[" << i << "][" << j << "]= " << A(i, j) << "  ";
// //    }
// //  std::cout << std::endl;
// //  }

// //  for (int i=0; i<4; i++) {
// //    for (int j=0; j<4; j++) {
// //      std::cout << "  InvA[" << i << "][" << j << "]= " << B(i, j) << "  ";
// //    }
// //  std::cout << std::endl;
// //  }
********/




//std::cout << surface << std::endl;


//Box test_box(min_point, max_point, nx, ny, nz, density, poission, young);

//Vector3D vec1(1.0);
//Vector3D vec2(2.0);
//Vector3D vec3(3.0);

//for (auto node_iter = test_box.getNodes().begin(); node_iter != test_box.getNodes().end(); node_iter++)
//for (int node_iter = 0; node_iter < test_box->getNodes().size(); node_iter++)
//{
//  NodeP cur_node = *node_iter;
//  NodeP cur_node = test_box->getNodes()[node_iter];
/*  std::cout << "Node number= " << cur_node->getID() << "   (" << cur_node->x()
                                                    << ", " << cur_node->y()
                                                    << ", " << cur_node->z() 
                                                    << ")" << " surface node= (" << cur_node->surfaceForce().x() 
                                                    << ", " << cur_node->surfaceForce().y() << ", " 
                                                    << cur_node->surfaceForce().z() << ") " << std::endl; */
/*  if (cur_node->getID()%3 == 1)
      cur_node->newDisplacement(vec1);
  if (cur_node->getID()%3 == 2)
      cur_node->newDisplacement(vec2);
  if (cur_node->getID()%3 == 0)
      cur_node->newDisplacement(vec3); */
//}


//Vector3D vec_zero(0.0);
//for (auto element_iter = test_box.getElements().begin(); element_iter != test_box.getElements().end(); element_iter++)
/*for (int element_iter = 0; element_iter < test_box->getElements().size(); element_iter++)
{
//  SymmetricMaterialTrilinearElementSP cur_element = *element_iter;
  SymmetricMaterialTrilinearElementSP cur_element = test_box->getElements()[element_iter];
  std::cout << "Element Number= " << cur_element->id() << std::endl << "(";
  for (int i=0; i < 8; i++)
  {
    NodeP cur_node = cur_element->triElement()->getArray()[i];
    if (!(cur_node->surfaceForce() == vec_zero))
       std::cout << "surface force of local node(" << i+1 << ")= (" << 
       cur_node->surfaceForce().x() << "," << cur_node->surfaceForce().y() << "," << cur_node->surfaceForce().z() 
        << ") " << std::endl; 
    if (i < 7)
    std::cout << cur_node->getID() << ", ";
    else 
    std::cout << cur_node->getID() << ")" << std::endl;
    
  }
}*/


//Solver::solverType type1 = Solver::solverType::Explicit;
//Solver solver(type1, 100, 0.0002, test_box);
//Solver::solverType type2 = Solver::solverType::Implicit;
//Solver solver(type2, 100, 0.0002, test_box);

//test_box->makeGlobalMatrices();

//int num = 8;


/*for (int jj = 0; jj < test_box->getGlobalStiffness().size(); jj++)
{
  for (int kk = 0; kk < test_box->getGlobalStiffness().size(); kk++)
  {
    std::cout << "(" << jj << "," << kk << ")="
              << test_box->getGlobalStiffness()[jj][kk] << " ";
//    if (kk%num == num-1)
//      std::cout << std::endl;
  }
}
std::cout << std::endl;*/ 



//std::cout << "Global Mass:" << std::endl;
//for (int jj = 0; jj < test_box->getGlobalMass().size(); jj++)
//{
//  for (int kk = 0; kk < test_box->getGlobalMass().size(); kk++)
//  {
//    if (test_box->getGlobalMass()[jj][kk] != 0)
//    std::cout << "(" << jj << "," << kk << ")="
//              << test_box->getGlobalMass()[jj][kk] << " ";
//    if (kk%num == num-1)
//      std::cout << std::endl;
//  }
//} 
//std::cout << std::endl << std::endl;

/*for (int jj = 0; jj < test_box->getGlobalLumpedMass().size(); jj++)
{
  for (int kk = 0; kk < test_box->getGlobalLumpedMass().size(); kk++)
  {
    if (test_box->getGlobalLumpedMass()[jj][kk] != 0.0)
    std::cout << "(" << jj << "," << kk << ")="
              << test_box->getGlobalLumpedMass()[jj][kk] << " ";
//    if (kk%num == num-1)
//      std::cout << std::endl;
  }
} 
std::cout << std::endl;*/

//std::cout << "External Force:" << std::endl;
//for (int ll = 0; ll < test_box->getExternalForce().size(); ll++)
//{
//  if ((ll%3 == 0) /*&& (ll != 0)*/) std::cout << std::endl << ll/3+1 << "   ";
//  std::cout << "(" << ll << ")= " << test_box->getExternalForce()[ll] << "  ";
//}
//std::cout << std::endl;

//for (int lll = 0; lll < test_box->getExternalForce().size(); lll++)
//{
//  if ((lll%3 == 0) /*&& (ll != 0)*/) std::cout << std::endl << lll/3+1 << "   ";
//  std::cout << "(" << lll << ")= " << test_box->getNewGlobDisp()[lll] << "  ";

//}
//std::cout << std::endl;


//std::vector<SymmetricMaterialTrilinearElementSP> elements = test_box->getElements();
//int elem_size = elements.size();

//for (int nn = 0; nn < 8; nn++)
//{
//  std::cout << test_box->getElements()[0]->triElement()->densityArray()[nn] << std::endl;
//}

/*std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> dis;*/


/*int M = RAND_MAX;
std::srand(std::time(NULL));
double xi1 = (2*std::rand()/M) - 1;
std::srand(std::time(NULL));
double xi2 = (2*std::rand()/M) - 1;
std::srand(std::time(NULL));
double xi3 = (2*std::rand()/M) - 1;*/


//double xi1 = -1.00;
//double xi2 = -0.22;
//double xi3 = 0.95;
//std::cout << "Det= " << test_box->getElements()[0]->triElement()->symJacobianMatrix(xi1, xi2, xi3).Determinant();
//std::cout << std::endl << "Norm Vector=(" << test_box->getElements()[1]->triElement()->symNormVecXi1Max(xi1, xi2).x()
//                       << "," << test_box->getElements()[1]->triElement()->symNormVecXi1Max(xi1, xi2).y()
//                       << "," << test_box->getElements()[1]->triElement()->symNormVecXi1Max(xi1, xi2).z()
//                       << ")" << std::endl; 
//std::cout << "Surface Jacob= " << test_box->getElements()[1]->triElement()->symNormVecXi1Max(xi1, xi2).length();
//std::cout << "  (" << xi1 << "," << xi2 << "," << xi3 << ")" << std::endl;
//std::cout << test_box->getElements()[0]->triElement()->symInterpolatedValue
//                       (test_box->getElements()[0]->triElement()->densityArray(), xi1, xi2, xi3) << std::endl;


//for (int mm = 0; mm < 24; mm++)
//  for (int pp = 0; pp < 24; pp++)
//     std::cout << "(" << mm << "," << pp << ")=" << test_box->getElements()[0]->massMat()(mm,pp) << "  ";
//std::cout << std::endl;


//for (int qq = 0; qq < 24; qq++)
//  std::cout << "(" << qq << ")=" << test_box->getElements()[1]->surfaceVec()(qq) << "  ";  
//std::cout << std::endl;

//std::cout << "dN1/dz(xi1=-0.25, xi2=0.75)=" 
//          << test_box->getElements()[0]->triElement()->symDerivativeTestFunctionToPosition(1,3,-0.25,0.75,0.25)
//          << std::endl;


//for (int rr = 0; rr < 23; rr++)
//  for (int ss = 0; ss < 23; ss++) {
//    if (test_box->getElements()[0]->stiffMat()(rr,ss) > 1)
//      std::cout << "(" << rr << "," << ss << ")="
//                << test_box->getElements()[0]->stiffMat()(rr,ss) << "  "; 
//  }
  


/*std::cout << "Lambda= " << test_box->getElements()[1]->firstLame() << "  Mui= " 
                        << test_box->getElements()[1]->shearModulus() << std::endl;

  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          std::cout << "C[" << i << "][" << j << "][" << k << "][" << l << "]= "
                    << test_box->getElements()[1]->tensor()(i, j, k, l) << "  ";
        }
      std::cout << std::endl;
      }
    }
  }*/



return 0;
}

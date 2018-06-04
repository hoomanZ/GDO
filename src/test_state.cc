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

//using namespace arma;
using namespace FiniteElement;

int main () {


Vector3D unitVector(-0.1581, 0.1645, -0.9736);
double magnitude = -0.01;



std::string nodeFile = "BoneDataNodes.txt";
std::string elementFile = "BoneDataElements.txt";

std::string tractionNodesFile = "downRadius.exnode";
std::string dirichletBoundaryNodesFile = "upRadius.exnode";



/********************************* Peridynamics **************************/

Vector3D min_point(0.0, 0.0, 0.0);
//Vector3D max_point(1.5, 1.5, 1.5);
//Vector3D max_point(0.01, 0.007, 0.007);

////Vector3D max_point(0.254, 0.0308, 0.0308);

Vector3D max_point(0.0616, 0.0308, 0.0308);


//Vector3D min(0.12059, 0.01248, 0.0154); // 0.00641
//Vector3D Max(0.13341, 0.01796, 0.0154);


////Vector3D min(0.12059, 0.01248, 0.03);
////Vector3D Max(0.13341, 0.01796, 0.03);


Vector3D min(0.03721, 0.01248, 0.03);
Vector3D Max(0.02439, 0.01796, 0.03);


//Vector3D min(0.02534, 0.01248, 0.0154);  // 0.254/8 - (0.01282/2)
//Vector3D Max(0.03816, 0.01796, 0.0154);  // 0.254/8 + (0.01282/2)



//int numPointsX = 21;
//int numPointsY = 15;
//int numPointsZ = 15;

//int numPointsX = 11;
//int numPointsY = 8;
//int numPointsZ = 8;


////int numPointsX = 100;
int numPointsX = 25;
int numPointsY = 13;
int numPointsZ = 13;



//double delT = 0.000025;
//double delT = 0.0001;
//double delT = 0.000048;
//double delT = 5*pow(10, -8);
double delT = pow(10, -8);
////double delT = pow(10, -9);
//////double delT = pow(10, -9);

////////double delT = pow(10, -10);


//double delT = pow(10, -8);
////double delT = pow(10, -7);


//double totalT = 2.0;
//double totalT = 0.02;
//double totalT = 0.1;
//double totalT = 0.0001;
//double totalT = 10*pow(10, -8);
////double totalT = 15*pow(10, -9);
//double totalT = 5*pow(10, -9);
//double totalT = 1.8*pow(10, -3);
//double totalT = 1.8*pow(10, -2);
//double totalT = 4.5*pow(10, -3);
//double totalT = 1.75*pow(10, -5);
//double totalT = 100*pow(10, -9);
//////double totalT = 1.75*pow(10, -4);
////double totalT = 1.75*pow(10, -6);
////////double totalT = 2.5*pow(10, -5)
///////////double totalT = 2.5*pow(10, -4);
//////double totalT = 3.29*pow(10, -6);
//double totalT = 4.5*pow(10, -6);

int k = 490;
//int k = 4;
//int k = 1;
//double totalT = k*pow(10, -7);

////double totalT = k*pow(10, -7);

double totalT = delT*k;

//double poission = 0.25;
double poission = 0.33;

//double density = 2440.0;
double density = 2700.0;
 
//double young = 72000000000;
double young = 70*pow(10, 9); // E = 70 GPa

//double fractureEnergy = 135;
//double yield = 255*pow(10, 6);

double criticalEquivalentStrain = 0.0042;


//double horizon = 0.22;
//double horizon = 0.3;
//double horizon = 0.002;
//double horizon = 0.002;
double horizon = 0.010262;



double bodyForce = 0.0;


Box::FixedDirichletBoundary boundary = Box::FixedDirichletBoundary::noFixedBoundary;
//Box::FixedDirichletBoundary boundary = Box::FixedDirichletBoundary::xMin;

Vector3D boundaryDisplacement(0.0);


//PeriBox::ForceSurface surface = PeriBox::ForceSurface::upDownX;
PeriBox::ForceSurface surface = PeriBox::ForceSurface::none;

//double surfaceForce = -30000;
//double surfaceForce = -35*pow(10, 6);
double surfaceForce = 0.0;
//double surfaceForce = -50000;
//double surfaceForce = -3000;


////PeriBox::ForceSurface velocitySurface = PeriBox::ForceSurface::upDownX;
PeriBox::ForceSurface velocitySurface = PeriBox::ForceSurface::tear;
////PeriBox::ForceSurface velocitySurface = PeriBox::ForceSurface::none;
////Vector3D velocity(0.01, 0.0, 0.0);
//////Vector3D velocity(100, 0.0, 0.0);
//Vector3D velocity(0.0, 100, 0.0);
Vector3D velocity(0.0, 20.0, 0.0);

////Vector3D velocity(0.0, 0.0, 0.0);


//Vector3D min(0.0, 0.7805, 0.0);
//Vector3D Max(1.5, 0.7905, 0.3);

////Vector3D min(0.12059, 0.01248, 0.0154);
////Vector3D Max(0.13341, 0.01796, 0.0154);

//double PreCrackLengthRatio = 0.5;


PeriBox::PeriType periType = PeriBox::PeriType::State;


//Matrix3D initialVelocityGradient(7874.0, 0.0, 0.0,
//                                 0.0, -2622.4, 0.0,
//                                 0.0, 0.0, -2622.4);


//Matrix3D initialVelocityGradient(78740.0, 0.0, 0.0,
//                                 0.0, -26220.4, 0.0,
//                                 0.0, 0.0, -26220.4);

////Matrix3D initialVelocityGradient(7874.0, 0.0, 0.0,
////                                 0.0, -2622.4, 0.0,
////                                 0.0, 0.0, -2622.4);

////double k = 100.0;
////initialVelocityGradient = initialVelocityGradient*k;

//PeriPMBMaterialP material(new PeriPMBMaterial(density, horizon, fractureEnergy, young));
PeriBoxP peri_box(new PeriBox(min_point, max_point, numPointsX, numPointsY, numPointsZ,
                              density, poission, young, criticalEquivalentStrain,
                              bodyForce, surfaceForce, surface,
                              boundary, boundaryDisplacement, horizon, periType,
                              velocitySurface, velocity));

//PeriBoxP peri_box(new PeriBox(min_point, max_point, numPointsX, numPointsY, numPointsZ,
//                              density, poission, young, criticalEquivalentStrain,
//                              bodyForce, surfaceForce, surface,
//                              boundary, boundaryDisplacement, horizon, periType,
//                              velocitySurface, velocity, initialVelocityGradient));


PeriCrackP peri_crack(new PeriCrack(min, Max));

//Peridynamics::BasisPolynomialOrders order = Peridynamics::BasisPolynomialOrders::One;
Peridynamics::BasisPolynomialOrders order = Peridynamics::BasisPolynomialOrders::Two;

////Peridynamics::GoverningEquationType type = Peridynamics::GoverningEquationType::StateBasedPeridynamics;
Peridynamics::GoverningEquationType governingType = Peridynamics::GoverningEquationType::Divergance;


////Peridynamics::TimeIntegration timeIntegration = Peridynamics::TimeIntegration::Explicit;
Peridynamics::TimeIntegration timeIntegration = Peridynamics::TimeIntegration::Implicit;



//PeridynamicsP test_prei(new Peridynamics(totalT, delT, peri_box));
//PeridynamicsP test_state(new Peridynamics(totalT, delT, peri_box));
//PeridynamicsP test_state(new Peridynamics(totalT, delT, peri_box, initialVelocityGradient));
////PeridynamicsP test_state(new Peridynamics(totalT, delT, peri_box, initialVelocityGradient, peri_crack)); 
//////PeridynamicsP test_state(new Peridynamics(totalT, delT, peri_box, initialVelocityGradient, peri_crack, order));

////PeridynamicsP test_state(new Peridynamics(totalT, delT, peri_box, initialVelocityGradient, peri_crack, order, type));
                  
//PeridynamicsP test_prei(new Peridynamics(totalT, delT, peri_box, peri_crack));




// ************************************************* Implicit Time Integration ***************************************

double alpha = 0.5;
double beta = 0.25;
double tolerance = 1.5*pow(10, -13);
//double tolerance = 1.5*pow(10, -20);
//double tolerance = 1.5*pow(10, -7);



//PeridynamicsP test_state(new Peridynamics(totalT, delT, peri_box, peri_crack, order, governingType));  
PeridynamicsP test_state(new Peridynamics(totalT, delT, peri_box, peri_crack, order, governingType, timeIntegration, alpha, beta, tolerance));  


// ********************************************************************************************************************


return 0;
}

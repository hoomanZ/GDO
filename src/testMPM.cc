#include <GeometryMath/Vector3D.h>
#include <MPMBoxPhases.h>
#include <MPMBoxP.h>
#include <BoxSP.h>

#include <MPMBox.h>
#include <Box.h>


using namespace FiniteElement;

int main () {


Vector3D min_point(0.0, 0.0, 0.0);
Vector3D max_point(1.5, 1.5, 1.5);

int xNumElem = 3;
int yNumElem = 3;
int zNumElem = 3;


int xNumPoints = 10;
int yNumPoints = 10;
int zNumPoints = 10;


double density = 2440.0;
double poission = 0.25;
//double young = 720000000.0;
double young = 720000.0;

//double implicitParameter = 1.0;
//double implicitParameter = 0.6;
double implicitParameter = 0.0;


////Vector3D pointForce(2250.0, 0.0, 0.0);
////Box::ForceSurface pointSurface = Box::ForceSurface::right;

//Box::FixedDirichletBoundary boundary = Box::FixedDirichletBoundary::noFixedBoundary;
Box::FixedDirichletBoundary boundary = Box::FixedDirichletBoundary::xMin;
Vector3D boundaryPosition(0.0, 0.0, 0.0);

//MPMBox::basisFunction basis = MPMBox::basisFunction::LinearHexahedral;
//MPMBox::basisFunction basis = MPMBox::basisFunction::PiecewiseLinear;
//MPMBox::basisFunction basis = MPMBox::basisFunction::BSpline;
MPMBox::basisFunction basis = MPMBox::basisFunction::GIMPfunction;


MPMBoxPhases::solverType type1 = MPMBoxPhases::solverType::Lumped;
//MPMBoxPhases::solverType type1 = MPMBoxPhases::solverType::Consistent;
//MPMBoxPhases::solverType type1 = MPMBoxPhases::solverType::LumpedUsingMomentum;
//MPMBoxPhases::solverType type1 = MPMBoxPhases::solverType::ConsistentUsingMomentum;


MPMBoxPhases::updateStressType update = MPMBoxPhases::updateStressType::USL;
//MPMBoxPhases::updateStressType update = MPMBoxPhases::updateStressType::MUSL;
//MPMBoxPhases::updateStressType update = MPMBoxPhases::updateStressType::USF;


int timeSteps = 10000;//10000;//0000;
double deltaT = 0.0005;
//double deltaT = 0.000002;
//double totalTime = 0.04;
double totalTime = timeSteps*deltaT;

double bodyForce = 0.0;
double surForce = 0.0;

//Box::ForceSurface tracSurface = Box::ForceSurface::none;



//BoxSP test_box(new Box(min_point, max_point,
//                       xNumElem, yNumElem, zNumElem,
//                       xNumPoints, yNumPoints, zNumPoints,
//                       density, poission, young,
//                       bodyForce, pointForce, surForce,
//                       tracSurface, pointSurface,
//                       boundary, boundaryPosition));


Vector3D unitVector(-0.1581, 0.1645, -0.9736);
//double magnitude = -0.01;
double magnitude = 50000;
Vector3D tractionForce = unitVector*magnitude;
Vector3D boundaryDisplacement(0.0);

std::string nodeFile = "BoneDataNodes.txt";
std::string elementFile = "BoneDataElements.txt";

std::string tractionNodesFile = "downRadius.exnode";
std::string dirichletBoundaryNodesFile = "upRadius.exnode";


//ComplicatedGeometrySP test_bone(new ComplicatedGeometry(nodeFile, elementFile, 
//                                tractionNodesFile, dirichletBoundaryNodesFile,
//                                density, young, poission, tractionForce, boundaryDisplacement, 2, 2, 2));

ComplicatedGeometrySP test_bone(new ComplicatedGeometry(nodeFile, elementFile, 
                                tractionNodesFile, dirichletBoundaryNodesFile,
                                density, young, poission, tractionForce, boundaryDisplacement,
                                5, 5, 5, xNumElem, yNumElem, zNumElem));



//MPMBoxP mpmBox(new MPMBox(test_box));
//MPMBoxP mpmBox(new MPMBox(test_box, basis));

MPMBoxP mpmBox(new MPMBox(test_bone, basis));


//MPMBoxPhases test_mpm(mpmBox, type1, deltaT, totalTime);
//MPMBoxPhases test_mpm(mpmBox, implicitParameter, deltaT, totalTime);
//MPMBoxPhases test_mpm(mpmBox, type1, implicitParameter, deltaT, totalTime);
MPMBoxPhases test_mpm(mpmBox, type1, update, implicitParameter, deltaT, totalTime);
//MPMBoxPhases test_mpm();




      
return 0;

}

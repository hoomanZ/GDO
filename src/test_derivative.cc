#include <derivative_operator.h>
#include <derivative_operatorP.h>
#include <cmath>

#include <eigen/Eigen/Eigenvalues>
#include <GeometryMath/Types.h>

using namespace FiniteElement;


double dx = 0.0001;

double point = M_PI;
double horizonSize = 0.0034;
int numRHSPoint = 54; 
int numLHSPoint = 54;

bool useTenTerms(true);
//bool useTenTerms(false);  


double
Function(double x)
{
  return 2*x*x*x + 4*x*x - 10*x + cos(3*x) - 15;

}


//******************************************************* Derivatives of Multi_variable functions ****************************************************




Vector3D pointVector(1.0, 2.0, 3.0);

//Vector3D pointVector(3.0, 0.0, 4.0);

//double lPercent = 0.0;
double lPercent = 100.0;

//double rPercent = 0.0;
double rPercent = 100.0;

double gridSize = 0.0002;
////double gridSize = 0.02; // for poster

double horizon(0.001);
////double horizon = 0.1;  // for poster


Vector3D vec1(1.0, 1.0, 0.0);
Vector3D vec2(0.0, 1.0, 1.0);
Vector3D vec3(1.0, 0.0, 1.0);
double gridSizeVec1Direction(gridSize);
double gridSizeVec2Direction(gridSize);
double gridSizeVec3Direction(gridSize);
double rhsPercentV1(rPercent);
double rhsPercentV2(rPercent);
double rhsPercentV3(rPercent);
double lhsPercentV1(lPercent);
double lhsPercentV2(lPercent);
double lhsPercentV3(lPercent);

//double lhsPercentV1(0.0);
//double lhsPercentV2(0.0);
//double lhsPercentV3(0.0);


//bool useBondDerivative(false);
bool useBondDerivative(true);

bool stateForce(false);
//bool stateForce(true);




Vector3D
example(Vector3D vec)
{
  Vector3D output(0.0);
  double x = vec.x();
  double y = vec.y();
  double z = vec.z();

//  output.x(2*x*x - 5*x*y + 7*y*z*z);
  output.x(x);

//  output.y(3*x*y*z + 6*pow(x,2)*pow(y,2)*pow(z,2) - 7*pow(z,4)*pow(y,3));
  output.y(y);

//  output.z(5*z*z + 12*x*x + 15*y*y);
  output.z(z);

  return output;

}



Matrix3
exampleDivergence(Vector3D vec)
{
  Matrix3 matrix = Matrix3::Zero();
  double x = vec.x();
  double y = vec.y();
  double z = vec.z();

  double r = pow(x*x + y*y +z*z, 0.5);

  matrix(0,0) = x/r;

  matrix(0,1) = y/r;

  matrix(0,2) = z/r;

  matrix(1,0) = x*x*y;

  matrix(1,1) = x*y*z;

  matrix(1,2) = (-1)*x*x*y*y;

  matrix(2,0) = 1/r;

  matrix(2,1) = 1/r;

  matrix(2,2) = 1/r;

  return matrix;

}


void
calculateDerivativeFunction()
{

  derivative_operatorP diff(new derivative_operator(&Function, point, horizonSize, numRHSPoint, numLHSPoint, useBondDerivative, useTenTerms)); 

}


void
calculateDerivativeVector()
{

  derivative_operatorP diff(new derivative_operator(&example, horizon, pointVector, vec1, vec2, vec3, 
                                                    gridSizeVec1Direction, gridSizeVec2Direction, gridSizeVec3Direction,
                                                    rhsPercentV1, rhsPercentV2, rhsPercentV3,
                                                    lhsPercentV1, lhsPercentV2, lhsPercentV3, useBondDerivative));

}


void
calculateDivergenceMatrix()
{

  derivative_operatorP diff(new derivative_operator(&exampleDivergence, horizon, pointVector, vec1, vec2, vec3, 
                                                    gridSizeVec1Direction, gridSizeVec2Direction, gridSizeVec3Direction,
                                                    rhsPercentV1, rhsPercentV2, rhsPercentV3,
                                                    lhsPercentV1, lhsPercentV2, lhsPercentV3, 
                                                    useBondDerivative, stateForce));

}



int main () {

//  calculateDerivativeFunction();

  calculateDerivativeVector();

//  calculateDivergenceMatrix();
   

////  derivative_operatorP diff(new derivative_operator()); 


}





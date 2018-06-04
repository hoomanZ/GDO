#include <derivative_operator.h>
#include <cmath>



using namespace FiniteElement;






derivative_operator::derivative_operator()
 : d_horizon_size(1.0), d_point(1.0, 2.0, 3.0),
   d_vec_1(1.0, 0.0, 0.0), 
   d_vec_2(0.0, 1.0, 0.0),
   d_vec_3(0.0, 0.0, 1.0),
   d_grid_size_vec1_direction(0.01),
   d_grid_size_vec2_direction(0.01),
   d_grid_size_vec3_direction(0.01),
   d_RHS_percent_vec1_direction(100),
   d_RHS_percent_vec2_direction(100),
   d_RHS_percent_vec3_direction(100),
   d_LHS_percent_vec1_direction(100),
   d_LHS_percent_vec2_direction(100),
   d_LHS_percent_vec3_direction(100)
{
  createNeighbourPoints();
  multiVariableFunction function = &derivative_operator::example;
  Matrix3 answer = derivative(function);
  printMatrix3(answer);  

}


derivative_operator::derivative_operator(double (*f)(double x), double point, double horizonSize, 
                                         int numRHSPoints, int numLHSPoints,
                                         bool useBondDerivative, bool useTenTerms)
 : d_horizon_size(1.0), d_point(1.0, 2.0, 3.0),
   d_vec_1(1.0, 0.0, 0.0), 
   d_vec_2(0.0, 1.0, 0.0),
   d_vec_3(0.0, 0.0, 1.0),
   d_grid_size_vec1_direction(0.01),
   d_grid_size_vec2_direction(0.01),
   d_grid_size_vec3_direction(0.01),
   d_RHS_percent_vec1_direction(100),
   d_RHS_percent_vec2_direction(100),
   d_RHS_percent_vec3_direction(100),
   d_LHS_percent_vec1_direction(100),
   d_LHS_percent_vec2_direction(100),
   d_LHS_percent_vec3_direction(100)
{

  if (useBondDerivative)
  {
    derivativeBondOfFunction(f, point, horizonSize, numRHSPoints, numLHSPoints);
  }
  else
  {
    if (useTenTerms)
    {
      tenDerivativeOfFunction(f, point, horizonSize, numRHSPoints, numLHSPoints);
    }
    else
    {
      derivativeOfFunction(f, point, horizonSize, numRHSPoints, numLHSPoints);
    }
  }
  
}



derivative_operator::derivative_operator(Vector3D (*f)(Vector3D vec),
                                         double horizon, Vector3D point,
                                         Vector3D vec1, Vector3D vec2, Vector3D vec3,
                                         double size1, double size2, double size3,
                                         double rhsPercentV1,  double rhsPercentV2,  double rhsPercentV3,
                                         double lhsPercentV1,  double lhsPercentV2,  double lhsPercentV3, bool useBondDerivative)
{
  setHorizonSize(horizon);
  setPoint(point);
  setVec1(vec1.normalized());
  setVec2(vec2.normalized());
  setVec3(vec3.normalized());
  setGridSizeVec1Direction(size1);
  setGridSizeVec2Direction(size2);
  setGridSizeVec3Direction(size3);
  setRHSPercentVec1Direction(rhsPercentV1);
  setRHSPercentVec2Direction(rhsPercentV2);
  setRHSPercentVec3Direction(rhsPercentV3);
  setLHSPercentVec1Direction(lhsPercentV1);
  setLHSPercentVec2Direction(lhsPercentV2);
  setLHSPercentVec3Direction(lhsPercentV3);


  createNeighbourPoints();

  Matrix3 answer = Matrix3::Zero();

  if (useBondDerivative)
  {
    answer = bondDerivative(f);
  }
  else
  {
    answer = derivative(f);
  }


  printMatrix3(answer);
}    



derivative_operator::derivative_operator(Matrix3 (*f)(Vector3D vec),
                                         double horizon, Vector3D point,
                                         Vector3D vec1, Vector3D vec2, Vector3D vec3,
                                         double size1, double size2, double size3,
                                         double rhsPercentV1,  double rhsPercentV2,  double rhsPercentV3,
                                         double lhsPercentV1,  double lhsPercentV2,  double lhsPercentV3, 
                                         bool useBondDerivative, bool state)
{
  setHorizonSize(horizon);
  setPoint(point);
  setVec1(vec1.normalized());
  setVec2(vec2.normalized());
  setVec3(vec3.normalized());
  setGridSizeVec1Direction(size1);
  setGridSizeVec2Direction(size2);
  setGridSizeVec3Direction(size3);
  setRHSPercentVec1Direction(rhsPercentV1);
  setRHSPercentVec2Direction(rhsPercentV2);
  setRHSPercentVec3Direction(rhsPercentV3);
  setLHSPercentVec1Direction(lhsPercentV1);
  setLHSPercentVec2Direction(lhsPercentV2);
  setLHSPercentVec3Direction(lhsPercentV3);


  createNeighbourPoints();

  Matrix3_1 answer = Matrix3_1::Zero();

  if (useBondDerivative)
  {
    if (!state)
    {
      answer = bondDivergenceMatrix3_3(f);
    }
    else
    {
      answer = forceStateMatrix3_3(f);
    }
  }
  else
  {
    answer = divergenceMatrix3_3(f);
  }


  printDivergence(answer);
}    




derivative_operator::~derivative_operator()
{

}




std::array <double, 4>
derivative_operator::Polybasis(double point, double neigbour)
{
  std::array <double, 4> poly{{0.0, 0.0, 0.0, 0.0}};
  double delX = neigbour - point;
  std::fill_n(poly.begin(), 1, 1);
  std::fill_n(poly.begin() + 1, 1, delX);
  std::fill_n(poly.begin() + 2, 1, delX*delX);
  std::fill_n(poly.begin() + 3, 1, delX*delX*delX);

  return poly;

}



Matrix4_4
derivative_operator::tensorProduct(std::array <double, 4> poly)
{
  Matrix4_4 tensor = Matrix4_4::Zero(); 
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      tensor(i, j) = poly[i]*poly[j];
    }
  }

  return tensor;

}



Matrix4_4
derivative_operator::shape(double point, double horizonSize, int numRHSPoints, int numLHSPoints)
{
  Matrix4_4 shapeTensor = Matrix4_4::Zero();  
  int n = numRHSPoints;
  int m = numLHSPoints;
  double dx = (2*horizonSize)/(m + n);
  for (int i = (-1)*m; i <= n; i++)
  { 
    if (i != 0)
    { 
      shapeTensor = shapeTensor + tensorProduct(Polybasis(point, point + i*dx))*dx;
    }

  }
  return shapeTensor;

}



double
derivative_operator::derivative(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints)
{
  int n = numRHSPoints;
  int m = numLHSPoints;
  double dx = (2*horizonSize)/(m + n);
  std::array <double, 4> poly{{0.0, 0.0, 0.0, 0.0}};
  double result = 0.0;  
  Matrix1_4 d(1, 4);
  d << 0, 1, 0, 0;
  Matrix4_1 P = Matrix4_1::Zero();
  Matrix4_1 integrand = Matrix4_1::Zero(); 
  Matrix4_4 invShape = shape(point, horizonSize, n, m).inverse();
  for (int i = (-1)*m; i <= n; i++)
  { 
    poly = Polybasis(point, point + i*dx);
    for (int j = 0; j < 4; j++)
    {    
      P(j, 0) = poly[j];      

    }
    if (i != 0)
    { 
      integrand = integrand + P*(*f)(point + i*dx)*dx;
    }

  }

  result = d*invShape*integrand;
  return result;

}



double
derivative_operator::secondDerivative(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints)
{
  int n = numRHSPoints;
  int m = numLHSPoints;
  double dx = (2*horizonSize)/(m + n);
  std::array <double, 4> poly{{0.0, 0.0, 0.0, 0.0}};
  double result = 0.0;  
  Matrix1_4 d(1, 4);
  d << 0, 0, 1, 0;
  Matrix4_1 P = Matrix4_1::Zero();
  Matrix4_1 integrand = Matrix4_1::Zero(); 
  Matrix4_4 invShape = shape(point, horizonSize, n, m).inverse();
  for (int i = (-1)*m; i <= n; i++)
  { 
    poly = Polybasis(point, point + i*dx);
    for (int j = 0; j < 4; j++)
    {    
      P(j, 0) = poly[j];      

    }
    if (i != 0)
    { 
      integrand = integrand + P*(*f)(point + i*dx)*dx;
    }

  }

  result = 2*d*invShape*integrand;
  return result;

}



double
derivative_operator::thirdDerivative(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints)
{
  int n = numRHSPoints;
  int m = numLHSPoints;
  double dx = (2*horizonSize)/(m + n);
  std::array <double, 4> poly{{0.0, 0.0, 0.0, 0.0}};
  double result = 0.0;  
  Matrix1_4 d(1, 4);
  d << 0, 0, 0, 1;
  Matrix4_1 P = Matrix4_1::Zero();
  Matrix4_1 integrand = Matrix4_1::Zero(); 
  Matrix4_4 invShape = shape(point, horizonSize, n, m).inverse();
  for (int i = (-1)*m; i <= n; i++)
  { 
    poly = Polybasis(point, point + i*dx);
    for (int j = 0; j < 4; j++)
    {    
      P(j, 0) = poly[j];      

    }
    if (i != 0)
    { 
      integrand = integrand + P*(*f)(point + i*dx)*dx;
    }

  }

  result = 6*d*invShape*integrand;
  return result;

}


void
derivative_operator::derivativeOfFunction(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints)
{
  double dx = horizonSize/(numRHSPoints + numLHSPoints);
  std::cout <<"horizonSize= " << horizonSize <<" dx= " << dx << ", numRHSPoints= " << numRHSPoints << ", numLHSPoints= " << numLHSPoints << std::endl;
  std::cout << "Derivative= " << derivative(f, point, horizonSize, numRHSPoints, numLHSPoints) << std::endl;
  std::cout << "Second derivative= " << secondDerivative(f, point, horizonSize, numRHSPoints, numLHSPoints) << std::endl;  
  std::cout << "Third derivative= " << thirdDerivative(f, point, horizonSize, numRHSPoints, numLHSPoints) << std::endl;  

}


//***************************************************** Derivatives of One_variable functions using ten terms of Taylor series ************************************



std::array <double, 10>
derivative_operator::tenPolybasis(double point, double neigbour)
{
  std::array <double, 10> poly{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  double delX = neigbour - point;
  std::fill_n(poly.begin(), 1, 1.0);
  std::fill_n(poly.begin() + 1, 1, delX);
  std::fill_n(poly.begin() + 2, 1, delX*delX);
  std::fill_n(poly.begin() + 3, 1, delX*delX*delX);
  std::fill_n(poly.begin() + 4, 1, pow(delX, 4));
  std::fill_n(poly.begin() + 5, 1, pow(delX, 5));
  std::fill_n(poly.begin() + 6, 1, pow(delX, 6));
  std::fill_n(poly.begin() + 7, 1, pow(delX, 7));
  std::fill_n(poly.begin() + 8, 1, pow(delX, 8));
  std::fill_n(poly.begin() + 9, 1, pow(delX, 9));

  return poly;

}



Matrix10
derivative_operator::tenTensorProduct(std::array <double, 10> poly)
{
  Matrix10 tensor = Matrix10::Zero(); 
  for (int i = 0; i < 10; i++)
  {
    for (int j = 0; j < 10; j++)
    {
      tensor(i, j) = poly[i]*poly[j];
    }
  }

  return tensor;

}



Matrix10
derivative_operator::tenShape(double point, double horizonSize, int numRHSPoints, int numLHSPoints)
{
  Matrix10 shapeTensor = Matrix10::Zero();  
  int n = numRHSPoints;
  int m = numLHSPoints;
  double dx = (2*horizonSize)/(m + n);
  for (int i = (-1)*m; i <= n; i++)
  { 
    if (i != 0)
    { 
      shapeTensor = shapeTensor + tenTensorProduct(tenPolybasis(point, point + i*dx))*dx;
    }

  }
  return shapeTensor;

}



double
derivative_operator::tenDerivative(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints)
{
  int n = numRHSPoints;
  int m = numLHSPoints;
  double dx = (2*horizonSize)/(m + n);
  std::array <double, 10> poly{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  double result = 0.0;  
  Matrix1_10 d(1, 10);
  d << 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
  Matrix10_1 P = Matrix10_1::Zero();
  Matrix10_1 integrand = Matrix10_1::Zero(); 
  Matrix10 invShape = tenShape(point, horizonSize, n, m).inverse();
  for (int i = (-1)*m; i <= n; i++)
  { 
    poly = tenPolybasis(point, point + i*dx);
    for (int j = 0; j < 10; j++)
    {    
      P(j, 0) = poly[j];      

    }
    if (i != 0)
    { 
      integrand = integrand + P*(*f)(point + i*dx)*dx;
    }

  }

  result = d*invShape*integrand;
  return result;

}



double
derivative_operator::tenSecondDerivative(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints)
{
  int n = numRHSPoints;
  int m = numLHSPoints;
  double dx = (2*horizonSize)/(m + n);
  std::array <double, 10> poly{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  double result = 0.0;  
  Matrix1_10 d(1, 10);
  d << 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
  Matrix10_1 P = Matrix10_1::Zero();
  Matrix10_1 integrand = Matrix10_1::Zero(); 
  Matrix10 invShape = tenShape(point, horizonSize, n, m).inverse();
  for (int i = (-1)*m; i <= n; i++)
  { 
    poly = tenPolybasis(point, point + i*dx);
    for (int j = 0; j < 10; j++)
    {    
      P(j, 0) = poly[j];      

    }
    if (i != 0)
    { 
      integrand = integrand + P*(*f)(point + i*dx)*dx;
    }

  }

  result = 2*d*invShape*integrand;
  return result;

}



double
derivative_operator::tenThirdDerivative(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints)
{
  int n = numRHSPoints;
  int m = numLHSPoints;
  double dx = (2*horizonSize)/(m + n);
  std::array <double, 10> poly{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  double result = 0.0;  
  Matrix1_10 d(1, 10);
  d << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
  Matrix10_1 P = Matrix10_1::Zero();
  Matrix10_1 integrand = Matrix10_1::Zero(); 
  Matrix10 invShape = tenShape(point, horizonSize, n, m).inverse();
  for (int i = (-1)*m; i <= n; i++)
  { 
    poly = tenPolybasis(point, point + i*dx);
    for (int j = 0; j < 10; j++)
    {    
      P(j, 0) = poly[j];      

    }
    if (i != 0)
    { 
      integrand = integrand + P*(*f)(point + i*dx)*dx;
    }

  }

  result = 6*d*invShape*integrand;
  return result;

}


void
derivative_operator::tenDerivativeOfFunction(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints)
{
  double dx = horizonSize/(numRHSPoints + numLHSPoints);
  std::cout << "horizonSize= " << horizonSize <<" dx= " << dx << ", numRHSPoints= " << numRHSPoints << ", numLHSPoints= " << numLHSPoints << std::endl;
  std::cout << "Derivative= " << tenDerivative(f, point, horizonSize, numRHSPoints, numLHSPoints) << std::endl;
  std::cout << "Second derivative= " << tenSecondDerivative(f, point, horizonSize, numRHSPoints, numLHSPoints) << std::endl;  
  std::cout << "Third derivative= " << tenThirdDerivative(f, point, horizonSize, numRHSPoints, numLHSPoints) << std::endl;  

}





// ********************************************************************** 1D Bond Derivatives ********************************************************
  

std::array <double, 3>
derivative_operator::PolybasisBond(double point, double neigbour)
{
  std::array <double, 3> poly{{0.0, 0.0, 0.0}};
  double delX = neigbour - point;
  std::fill_n(poly.begin(), 1, delX);
  std::fill_n(poly.begin() + 1, 1, delX*delX);
  std::fill_n(poly.begin() + 2, 1, delX*delX*delX);
  return poly;

}


Matrix3
derivative_operator::tensorProductBond(std::array <double, 3> poly)
{
  Matrix3 tensor = Matrix3::Zero(); 
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      tensor(i, j) = poly[i]*poly[j];
    }
  }

  return tensor;

}



Matrix3
derivative_operator::shapeBond(double point, double dx, int numRHSPoints, int numLHSPoints)
{
  Matrix3 shapeTensor = Matrix3::Zero();  
  int n = numRHSPoints;
  int m = numLHSPoints;
  for (int i = (-1)*m; i <= n; i++)
  { 
    if (i != 0)
    { 
      shapeTensor = shapeTensor + tensorProductBond(PolybasisBond(point, point + i*dx))*dx;
    }

  }
  return shapeTensor;

}



double
derivative_operator::derivativeBond(double (*f)(double x), double point, double dx, int numRHSPoints, int numLHSPoints)
{
  int n = numRHSPoints;
  int m = numLHSPoints;
  std::array <double, 3> poly{{0.0, 0.0, 0.0}};
  double result = 0.0;  
  Matrix1_3 d(1, 3);
  d << 1, 0, 0;
  Matrix3_1 P = Matrix3_1::Zero();
  Matrix3_1 integrand = Matrix3_1::Zero(); 
  Matrix3 invShape = shapeBond(point, dx, n, m).inverse();
  for (int i = (-1)*m; i <= n; i++)
  { 
    poly = PolybasisBond(point, point + i*dx);
    for (int j = 0; j < 3; j++)
    {    
      P(j, 0) = poly[j];      

    }
    if (i != 0)
    { 
      integrand = integrand + P*((*f)(point + i*dx) - (*f)(point))*dx;
    }

  }

  result = d*invShape*integrand;
  return result;

}


double
derivative_operator::secondDerivativeBond(double (*f)(double x), double point, double dx, int numRHSPoints, int numLHSPoints)
{
  int n = numRHSPoints;
  int m = numLHSPoints;
  std::array <double, 3> poly{{0.0, 0.0, 0.0}};
  double result = 0.0;  
  Matrix1_3 d(1, 3);
  d << 0, 1, 0;
  Matrix3_1 P = Matrix3_1::Zero();
  Matrix3_1 integrand = Matrix3_1::Zero(); 
  Matrix3 invShape = shapeBond(point, dx, n, m).inverse();
  for (int i = (-1)*m; i <= n; i++)
  { 
    poly = PolybasisBond(point, point + i*dx);
    for (int j = 0; j < 3; j++)
    {    
      P(j, 0) = poly[j];      

    }
    if (i != 0)
    { 
      integrand = integrand + P*((*f)(point + i*dx) - (*f)(point))*dx;
    }

  }

  result = 2*d*invShape*integrand;
  return result;

}



double
derivative_operator::thirdDerivativeBond(double (*f)(double x), double point, double dx, int numRHSPoints, int numLHSPoints)
{
  int n = numRHSPoints;
  int m = numLHSPoints;
  std::array <double, 3> poly{{0.0, 0.0, 0.0}};
  double result = 0.0;  
  Matrix1_3 d(1, 3);
  d << 0, 0, 1;
  Matrix3_1 P = Matrix3_1::Zero();
  Matrix3_1 integrand = Matrix3_1::Zero(); 
  Matrix3 invShape = shapeBond(point, dx, n, m).inverse();
  for (int i = (-1)*m; i <= n; i++)
  { 
    poly = PolybasisBond(point, point + i*dx);
    for (int j = 0; j < 3; j++)
    {    
      P(j, 0) = poly[j];      

    }
    if (i != 0)
    { 
      integrand = integrand + P*((*f)(point + i*dx) - (*f)(point))*dx;
    }

  }

  result = 6*d*invShape*integrand;
  return result;

}



void
derivative_operator::derivativeBondOfFunction(double (*f)(double x), double point, double dx, int numRHSPoints, int numLHSPoints)
{

  std::cout << "dx= " << dx << ", numRHSPoints= " << numRHSPoints << ", numLHSPoints= " << numLHSPoints << std::endl;
  std::cout << "Derivative= " << derivativeBond(f, point, dx, numRHSPoints, numLHSPoints) << std::endl;
  std::cout << "Second derivative= " << secondDerivativeBond(f, point, dx, numRHSPoints, numLHSPoints) << std::endl;  
  std::cout << "Third derivative= " << thirdDerivativeBond(f, point, dx, numRHSPoints, numLHSPoints) << std::endl;  

}




//**************************************************************************** Derivatives of Multi_variable functions ***********************************************************


double
derivative_operator::makeBetweenZeroAndOne(double number)
{
  if (number < 0 || number > 100)
  {
    return 1.0;
  }
  else
  {
    return number/100;
  }
}


int
derivative_operator::extVal(double number, double percent)
{
  return number*makeBetweenZeroAndOne(percent);

}




void
derivative_operator::createNeighbourPoints()
{
  std::vector <Vector3D> neigbours;
  int maxN1 = getHorizonSize()/getGridSizeVec1Direction();
  maxN1 = extVal(1*maxN1, getRHSPercentVec1Direction());   
  int maxN2 = getHorizonSize()/getGridSizeVec2Direction();
  maxN2 = extVal(1*maxN2, getRHSPercentVec1Direction());
  int maxN3 = getHorizonSize()/getGridSizeVec3Direction();
  maxN3 = extVal(1*maxN3, getRHSPercentVec1Direction());

  int minN1 = extVal((-1)*maxN1, getLHSPercentVec1Direction());
  int minN2 = extVal((-1)*maxN2, getLHSPercentVec2Direction());
  int minN3 = extVal((-1)*maxN3, getLHSPercentVec1Direction());

  for (int i = minN1; i <= maxN1; i++)
  {
    for (int j = minN2; j <= maxN2; j++)
    {
      for (int k = minN3; k <= maxN3; k++)
      {
         Vector3D vector = getVec1()*(i*getGridSizeVec1Direction())+
                           getVec2()*(j*getGridSizeVec2Direction())+
                           getVec3()*(k*getGridSizeVec3Direction());
         if (vector.length() < getHorizonSize())
         {
           neigbours.push_back(getPoint() + vector);
         
         }            

       }

    }


  }

  setNeighbourPoints(neigbours);

}


Vector3D
derivative_operator::example(Vector3D vec)
{
  Vector3D output(0.0);
  output.x(vec.x()*vec.x());
  output.y(vec.y()*vec.y());
  output.z(vec.z()*vec.z());
  return output;

}




std::array <double, 10>  
derivative_operator::multiPolybasis(Vector3D neigbour)
{

  std::array <double, 10> poly;
  Vector3D sai = neigbour - getPoint();
  
  std::fill_n(poly.begin(), 1, 1);
  std::fill_n(poly.begin() + 1, 1, sai.x());
  std::fill_n(poly.begin() + 2, 1, sai.y());
  std::fill_n(poly.begin() + 3, 1, sai.z());
  std::fill_n(poly.begin() + 4, 1, sai.x()*sai.x());
  std::fill_n(poly.begin() + 5, 1, sai.y()*sai.y());
  std::fill_n(poly.begin() + 6, 1, sai.z()*sai.z());
  std::fill_n(poly.begin() + 7, 1, sai.x()*sai.y());
  std::fill_n(poly.begin() + 8, 1, sai.x()*sai.z());
  std::fill_n(poly.begin() + 9, 1, sai.y()*sai.z());

  return poly;

}


Matrix10
derivative_operator::multiTensorProduct(std::array <double, 10> poly)
{
  Matrix10 tensor = Matrix10::Zero(); 
  for (int i = 0; i < 10; i++)
  {
    for (int j = 0; j < 10; j++)
    {
      tensor(i, j) = poly[i]*poly[j];
    }
  }

  return tensor;
  
}


Matrix10 
derivative_operator::multiShape()
{
  double volume = getGridSizeVec1Direction()*getGridSizeVec2Direction()*getGridSizeVec3Direction();
  Matrix10 shapeTensor = Matrix10::Zero();  
  std::vector <Vector3D> neighbours = getNeighbourPoints();
  int size = neighbours.size();

  for (int i = 0; i < size; i++)
  { 
    Vector3D cur_neigh = neighbours[i];     
    shapeTensor = shapeTensor + multiTensorProduct(multiPolybasis(cur_neigh))*volume;
   
  }
  return shapeTensor;

}


Matrix3
derivative_operator::derivative(Vector3D (*f)(Vector3D vec))
{
  double volume = getGridSizeVec1Direction()*getGridSizeVec2Direction()*getGridSizeVec3Direction();

  Matrix3 result = Matrix3::Zero();
  std::array <double, 10> poly{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};  
  Matrix10_3 d(10, 3);
  d << 0, 0, 0,
       1, 0, 0,
       0, 1, 0,
       0, 0, 1,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0;

  Matrix3_10 integrand = Matrix3_10::Zero(); 
  Matrix10 invShape = multiShape().inverse();
  std::vector <Vector3D> neighbours = getNeighbourPoints();
  int size = neighbours.size();

  for (int i = 0; i < size; i++)
  { 
    Vector3D cur_neigh = neighbours[i];     
    poly = multiPolybasis(cur_neigh);

    for (int i1 = 0; i1 < 3; i1++)
    { 
      for (int i2 = 0; i2 < 10; i2++)
      {
        integrand(i1, i2) = integrand(i1, i2) + (*f)(cur_neigh)[i1]*poly[i2]*volume;
      }
    }    
    
  }

  result = integrand*invShape*d;
  return result;

}


Matrix3
derivative_operator::derivative(multiVariableFunction f)
{
  double volume = getGridSizeVec1Direction()*getGridSizeVec2Direction()*getGridSizeVec3Direction();

  Matrix3 result = Matrix3::Zero();
  std::array <double, 10> poly{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};  
  Matrix10_3 d(10, 3);
  d << 0, 0, 0, 
       1, 0, 0,
       0, 1, 0,
       0, 0, 1,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0;

  Matrix3_10 integrand = Matrix3_10::Zero(); 
  Matrix10 invShape = multiShape().inverse();
  std::vector <Vector3D> neighbours = getNeighbourPoints();
  int size = neighbours.size();

  for (int i = 0; i < size; i++)
  { 
    Vector3D cur_neigh = neighbours[i];     
    poly = multiPolybasis(cur_neigh);

    for (int i1 = 0; i1 < 3; i1++)
    { 
      for (int i2 = 0; i2 < 10; i2++)
      {
        integrand(i1, i2) = integrand(i1, i2) + (this->*f)(cur_neigh)[i1]*poly[i2]*volume;
      }
    }    
    
  }

  result = integrand*invShape*d;
  return result;

}


Matrix3_1
derivative_operator::divergenceMatrix3_3(Matrix3 (*f)(Vector3D vec))
{
  double volume = getGridSizeVec1Direction()*getGridSizeVec2Direction()*getGridSizeVec3Direction();
  Matrix3_1 integrand = Matrix3_1::Zero();
  Matrix3_10 d(3, 10);

  d << 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 1, 0, 0, 0, 0, 0, 0;

  std::array <double, 10> polyBasis{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  Matrix10 invShape = multiShape().inverse();
  std::vector <Vector3D> neighbours = getNeighbourPoints();
  int size = neighbours.size();

  Matrix10_1 poly = Matrix10_1::Zero();

  for (int i = 0; i < size; i++)
  { 
    Vector3D cur_neigh = neighbours[i];
    polyBasis = multiPolybasis(cur_neigh);     

    for (int j = 0; j < 10; j++)
    {
       poly(j, 0) = polyBasis[j]; 
     
    }

    integrand = integrand + ((*f)(cur_neigh)*d*invShape*poly)*volume;;    
  }

  return integrand;

}



void
derivative_operator::printMatrix3(Matrix3 matrix)
{
  std::cout << " Partial derivative of f_x WRT x = " << matrix(0,0) << std::endl; 
  std::cout << " Partial derivative of f_x WRT y = " << matrix(0,1) << std::endl; 
  std::cout << " Partial derivative of f_x WRT z = " << matrix(0,2) << std::endl; 
  std::cout << " Partial derivative of f_y WRT x = " << matrix(1,0) << std::endl; 
  std::cout << " Partial derivative of f_y WRT y = " << matrix(1,1) << std::endl; 
  std::cout << " Partial derivative of f_y WRT z = " << matrix(1,2) << std::endl; 
  std::cout << " Partial derivative of f_z WRT x = " << matrix(2,0) << std::endl; 
  std::cout << " Partial derivative of f_z WRT y = " << matrix(2,1) << std::endl; 
  std::cout << " Partial derivative of f_z WRT z = " << matrix(2,2) << std::endl; 
}


void
derivative_operator::printDivergence(Matrix3_1 divergence)
{
  std::cout << " Divergance of the  matrix is:" << std::endl;
  std::cout << "(" << divergence(0,0) << ", " << divergence(1,0) << ", " << divergence(2,0) << ")" << std::endl;

}




// ********************************************************************** 3D Bond Derivatives ********************************************************


std::array <double, 9>  
derivative_operator::bondMultiPolybasis(Vector3D neigbour)
{

  std::array <double, 9> poly;
  Vector3D sai = neigbour - getPoint();
  
  std::fill_n(poly.begin(), 1, sai.x());
  std::fill_n(poly.begin() + 1, 1, sai.y());
  std::fill_n(poly.begin() + 2, 1, sai.z());
  std::fill_n(poly.begin() + 3, 1, sai.x()*sai.x());
  std::fill_n(poly.begin() + 4, 1, sai.y()*sai.y());
  std::fill_n(poly.begin() + 5, 1, sai.z()*sai.z());
  std::fill_n(poly.begin() + 6, 1, sai.x()*sai.y());
  std::fill_n(poly.begin() + 7, 1, sai.x()*sai.z());
  std::fill_n(poly.begin() + 8, 1, sai.y()*sai.z());

  return poly;

}


Matrix9
derivative_operator::bondMultiTensorProduct(std::array <double, 9> poly)
{
  Matrix9 tensor = Matrix9::Zero(); 
  for (int i = 0; i < 9; i++)
  {
    for (int j = 0; j < 9; j++)
    {
      tensor(i, j) = poly[i]*poly[j];
    }
  }

  return tensor;
  
}


Matrix9 
derivative_operator::bondMultiShape()
{
  double volume = getGridSizeVec1Direction()*getGridSizeVec2Direction()*getGridSizeVec3Direction();
  Matrix9 shapeTensor = Matrix9::Zero();  
  std::vector <Vector3D> neighbours = getNeighbourPoints();
  int size = neighbours.size();

  for (int i = 0; i < size; i++)
  { 
    Vector3D cur_neigh = neighbours[i];     
    shapeTensor = shapeTensor + bondMultiTensorProduct(bondMultiPolybasis(cur_neigh))*volume;
   
  }
  return shapeTensor;

}


Matrix3
derivative_operator::bondDerivative(Vector3D (*f)(Vector3D vec))
{
  double volume = getGridSizeVec1Direction()*getGridSizeVec2Direction()*getGridSizeVec3Direction();

  Matrix3 result = Matrix3::Zero();
  std::array <double, 9> poly{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};  
  Matrix9_3 d(9, 3);

  d << 1, 0, 0,
       0, 1, 0,
       0, 0, 1,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0;

  Matrix3_9 integrand = Matrix3_9::Zero(); 
  Matrix9 invShape = bondMultiShape().inverse();
  std::vector <Vector3D> neighbours = getNeighbourPoints();
  int size = neighbours.size();

  for (int i = 0; i < size; i++)
  { 
    Vector3D cur_neigh = neighbours[i];     
    poly = bondMultiPolybasis(cur_neigh);

    for (int i1 = 0; i1 < 3; i1++)
    { 
      for (int i2 = 0; i2 < 9; i2++)
      {
        integrand(i1, i2) = integrand(i1, i2) + ((*f)(cur_neigh) - (*f)(getPoint()))[i1]*poly[i2]*volume;
      }
    }    
    
  }

  result = integrand*invShape*d;
  return result;

}


Matrix3
derivative_operator::bondDerivative(multiVariableFunction f)
{
  double volume = getGridSizeVec1Direction()*getGridSizeVec2Direction()*getGridSizeVec3Direction();

  Matrix3 result = Matrix3::Zero();
  std::array <double, 9> poly{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};  
  Matrix9_3 d(9, 3);
  d << 1, 0, 0,
       0, 1, 0,
       0, 0, 1,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0;

  Matrix3_9 integrand = Matrix3_9::Zero(); 
  Matrix9 invShape = bondMultiShape().inverse();
  std::vector <Vector3D> neighbours = getNeighbourPoints();
  int size = neighbours.size();

  for (int i = 0; i < size; i++)
  { 
    Vector3D cur_neigh = neighbours[i];     
    poly = bondMultiPolybasis(cur_neigh);

    for (int i1 = 0; i1 < 3; i1++)
    { 
      for (int i2 = 0; i2 < 9; i2++)
      {
        integrand(i1, i2) = integrand(i1, i2) + ((this->*f)(cur_neigh) - (this->*f)(getPoint()))[i1]*poly[i2]*volume;
      }
    }    
    
  }

  result = integrand*invShape*d;
  return result;

}



Matrix3_1
derivative_operator::bondDivergenceMatrix3_3(Matrix3 (*f)(Vector3D vec))
{
  double volume = getGridSizeVec1Direction()*getGridSizeVec2Direction()*getGridSizeVec3Direction();
  Matrix3_1 integrand = Matrix3_1::Zero();
  Matrix3_9 d(3, 9);

  d << 1, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 1, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 1, 0, 0, 0, 0, 0, 0;

  std::array <double, 9> polyBasis{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  Matrix9 invShape = bondMultiShape().inverse();
  std::vector <Vector3D> neighbours = getNeighbourPoints();
  int size = neighbours.size();

  Matrix9_1 poly = Matrix9_1::Zero();

  for (int i = 0; i < size; i++)
  { 
    Vector3D cur_neigh = neighbours[i];
    polyBasis = bondMultiPolybasis(cur_neigh);     

    for (int j = 0; j < 9; j++)
    {
       poly(j, 0) = polyBasis[j]; 
     
    }

    integrand = integrand + (((*f)(cur_neigh) - (*f)(getPoint()))*d*invShape*poly)*volume;    
  }

  return integrand;

}



Matrix3_1
derivative_operator::forceStateMatrix3_3(Matrix3 (*f)(Vector3D vec))
{
  double volume = getGridSizeVec1Direction()*getGridSizeVec2Direction()*getGridSizeVec3Direction();
  Matrix3_1 integrand = Matrix3_1::Zero();
  Matrix3_9 d(3, 9);

  d << 1, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 1, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 1, 0, 0, 0, 0, 0, 0;

  std::array <double, 9> polyBasis{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  Matrix9 invShape = bondMultiShape().inverse();
  std::vector <Vector3D> neighbours = getNeighbourPoints();
  int size = neighbours.size();

  Matrix9_1 poly = Matrix9_1::Zero();
  Matrix9_1 neighbourPoly = Matrix9_1::Zero();

  for (int i = 0; i < size; i++)
  { 
    Vector3D cur_neigh = neighbours[i];
    polyBasis = bondMultiPolybasis(cur_neigh);     

    for (int j = 0; j < 9; j++)
    {
       poly(j, 0) = polyBasis[j];
       if (j < 3)
       { 
         neighbourPoly(j, 0) = (-1)*polyBasis[j];
       }
       else
       {
         neighbourPoly(j, 0) = polyBasis[j];
       }         
               
    }

    integrand = integrand + ((*f)(getPoint())*d*invShape*poly)*volume - ((*f)(cur_neigh)*d*invShape*neighbourPoly)*volume;    
  }

  return integrand;

}






/*********************************************************************************************************************************************************************************
double
derivative_operator::derivative(scalar f, double point, double dx, int numOfPoints)
{
  int n = numOfPoints;
  std::array <double, 4> poly{{0.0, 0.0, 0.0, 0.0}};
  double result = 0.0;  
  Matrix1_4 d(1, 4);
  d << 0, 1, 0, 0;
  Matrix4_1 P = Matrix4_1::Zero();
  Matrix4_1 integrand = Matrix4_1::Zero(); 
  Matrix4_4 invShape = shape(point, dx, n).inverse();
  for (int i = (-1)*n; i <= n; i++)
  { 
    poly = Polybasis(point, point + i*dx);
    for (int j = 0; j < 4; j++)
    {    
      P(j, 0) = poly[j];      

    }
    if (i != 0)
    { 
      integrand = integrand + P*(this->*f)(point + i*dx)*dx;
    }

  }

  result = d*invShape*integrand;
  return result;

}


double
derivative_operator::secondDerivative(scalar f, double point, double dx, int numOfPoints)
{
  int n = numOfPoints;
  std::array <double, 4> poly{{0.0, 0.0, 0.0, 0.0}};
  double result = 0.0;  
  Matrix1_4 d(1, 4);
  d << 0, 0, 1, 0;
  Matrix4_1 P = Matrix4_1::Zero();
  Matrix4_1 integrand = Matrix4_1::Zero(); 
  Matrix4_4 invShape = shape(point, dx, n).inverse();
  for (int i = (-1)*n; i <= n; i++)
  { 
    poly = Polybasis(point, point + i*dx);
    for (int j = 0; j < 4; j++)
    {    
      P(j, 0) = poly[j];      

    }
    if (i != 0)
    { 
      integrand = integrand + P*(this->*f)(point + i*dx)*dx;
    }

  }

  result = 2*d*invShape*integrand;
  return result;

}


double
derivative_operator::thirdDerivative(scalar f, double point, double dx, int numOfPoints)
{
  int n = numOfPoints;
  std::array <double, 4> poly{{0.0, 0.0, 0.0, 0.0}};
  double result = 0.0;  
  Matrix1_4 d(1, 4);
  d << 0, 0, 0, 1;
  Matrix4_1 P = Matrix4_1::Zero();
  Matrix4_1 integrand = Matrix4_1::Zero(); 
  Matrix4_4 invShape = shape(point, dx, n).inverse();
  for (int i = (-1)*n; i <= n; i++)
  { 
    poly = Polybasis(point, point + i*dx);
    for (int j = 0; j < 4; j++)
    {    
      P(j, 0) = poly[j];      

    }
    if (i != 0)
    { 
      integrand = integrand + P*(this->*f)(point + i*dx)*dx;
    }

  }

  result = 6*d*invShape*integrand;
  return result;

}



double 
derivative_operator::Function(double x)
{
  return 2*x*x*x + 4*x*x - 10*x + cos(3*x) - 15;

}


void
derivative_operator::derivativeOfFunction(double point, double dx, int numOfPoints)
{
  scalar function = &derivative_operator::Function;
//  double point = 2.0;
//  double dx = 0.0001;
//  int numOfPoint = 20;
  std::cout << "Derivative= " << derivative(function, point, dx, numOfPoints) << std::endl;
  std::cout << "Second derivative= " << secondDerivative(function, point, dx, numOfPoints) << std::endl;  
  std::cout << "Third derivative= " << thirdDerivative(function, point, dx, numOfPoints) << std::endl;  


}

*/



#ifndef FINITEELEMENT_DERIVATIVE_OPERATOR_H__
#define FINITEELEMENT_DERIVATIVE_OPERATOR_H__



#include <GeometryMath/Vector3D.h>
#include <GeometryMath/Point3D.h>
#include <GeometryMath/Matrix3D.h>
#include <GeometryMath/Types.h>

#include <eigen/Eigen/Eigenvalues>
#include <string>


namespace FiniteElement {
  
  class derivative_operator {

    public:

      derivative_operator();
      
      derivative_operator(double (*f)(double x), double point, double horizonSize,
                          int numRHSPoints, int numLHSPoints, 
                          bool useBondDerivative, bool useTenTerms);

      derivative_operator(Vector3D (*f)(Vector3D vec), 
                          double horizon, Vector3D point,
                          Vector3D vec1, Vector3D vec2, Vector3D vec3,
                          double size1, double size2, double size3,
                          double rhsPercentV1,  double rhsPercentV2,  double rhsPercentV3,
                          double lhsPercentV1,  double lhsPercentV2,  double lhsPercentV3, bool useBondDerivative);    

      derivative_operator(Matrix3 (*f)(Vector3D vec), 
                          double horizon, Vector3D point,
                          Vector3D vec1, Vector3D vec2, Vector3D vec3,
                          double size1, double size2, double size3,
                          double rhsPercentV1,  double rhsPercentV2,  double rhsPercentV3,
                          double lhsPercentV1,  double lhsPercentV2,  double lhsPercentV3,
                          bool useBondDerivative, bool state);    


      ~derivative_operator();

      typedef Vector3D (derivative_operator::*multiVariableFunction)(Vector3D vec);
      typedef double (derivative_operator::*scalar)(double x);

//      inline void setFunction(const scalar& func) {d_function = func;}
//      inline scalar getPeriBox() const {return d_function;}
      
//***************************************************** Derivatives of One_variable functions using four terms of Taylor series ********************************

      std::array <double, 4>  Polybasis(double point, double neigbour);
      Matrix4_4 tensorProduct(std::array <double, 4> poly);
      Matrix4_4 shape(double point, double horizonSize, int numRHSPoints, int numLHSPoints);

      double derivative(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints);
      double secondDerivative(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints);
      double thirdDerivative(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints);

      void derivativeOfFunction(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints);


//***************************************************** Derivatives of One_variable functions using ten terms of Taylor series ************************************

      std::array <double, 10>  tenPolybasis(double point, double neigbour);
      Matrix10 tenTensorProduct(std::array <double, 10> poly);
      Matrix10 tenShape(double point, double horizonSize, int numRHSPoints, int numLHSPoints);

      double tenDerivative(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints);
      double tenSecondDerivative(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints);
      double tenThirdDerivative(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints);

      void tenDerivativeOfFunction(double (*f)(double x), double point, double horizonSize, int numRHSPoints, int numLHSPoints);


//****************************************** 1D Bond Derivatives **************************


      std::array <double, 3>  PolybasisBond(double point, double neigbour);
      Matrix3 tensorProductBond(std::array <double, 3> poly);
      Matrix3 shapeBond(double point, double dx, int numRHSPoints, int numLHSPoints);

      double derivativeBond(double (*f)(double x), double point, double dx, int numRHSPoints, int numLHSPoints);
      double secondDerivativeBond(double (*f)(double x), double point, double dx, int numRHSPoints, int numLHSPoints);
      double thirdDerivativeBond(double (*f)(double x), double point, double dx, int numRHSPoints, int numLHSPoints);

      void derivativeBondOfFunction(double (*f)(double x), double point, double dx, int numRHSPoints, int numLHSPoints);


//********************************************************************* Derivatives of Multi_variable functions *********************************************



      inline void setHorizonSize(const double& horizon) {d_horizon_size = horizon;}
      inline double getHorizonSize() const {return d_horizon_size;}

      inline void setPoint(const Vector3D& point) {d_point = point;}
      inline Vector3D getPoint() const {return d_point;}

      inline void setVec1(const Vector3D& vec1) {d_vec_1 = vec1;}
      inline Vector3D getVec1() const {return d_vec_1;}

      inline void setVec2(const Vector3D& vec2) {d_vec_2 = vec2;}
      inline Vector3D getVec2() const {return d_vec_2;}

      inline void setVec3(const Vector3D& vec3) {d_vec_3 = vec3;}
      inline Vector3D getVec3() const {return d_vec_3;}

      inline void setGridSizeVec1Direction(const double& size) {d_grid_size_vec1_direction = size;}
      inline double getGridSizeVec1Direction() const {return d_grid_size_vec1_direction;}

      inline void setGridSizeVec2Direction(const double& size) {d_grid_size_vec2_direction = size;}
      inline double getGridSizeVec2Direction() const {return d_grid_size_vec2_direction;}

      inline void setGridSizeVec3Direction(const double& size) {d_grid_size_vec3_direction = size;}
      inline double getGridSizeVec3Direction() const {return d_grid_size_vec3_direction;}

      inline void setRHSPercentVec1Direction(const double& percent) {d_RHS_percent_vec1_direction = percent;}
      inline double getRHSPercentVec1Direction() const {return d_RHS_percent_vec1_direction;}

      inline void setRHSPercentVec2Direction(const double& percent) {d_RHS_percent_vec2_direction = percent;}
      inline double getRHSPercentVec2Direction() const {return d_RHS_percent_vec2_direction;}

      inline void setRHSPercentVec3Direction(const double& percent) {d_RHS_percent_vec3_direction = percent;}
      inline double getRHSPercentVec3Direction() const {return d_RHS_percent_vec3_direction;}

      inline void setLHSPercentVec1Direction(const double& percent) {d_LHS_percent_vec1_direction = percent;}
      inline double getLHSPercentVec1Direction() const {return d_LHS_percent_vec1_direction;}

      inline void setLHSPercentVec2Direction(const double& percent) {d_LHS_percent_vec2_direction = percent;}
      inline double getLHSPercentVec2Direction() const {return d_LHS_percent_vec2_direction;}

      inline void setLHSPercentVec3Direction(const double& percent) {d_LHS_percent_vec3_direction = percent;}
      inline double getLHSPercentVec3Direction() const {return d_LHS_percent_vec3_direction;}


      inline void setNeighbourPoints(const std::vector <Vector3D>& vectors) {d_neighbour_points = vectors;}
      inline std::vector <Vector3D> getNeighbourPoints() const {return d_neighbour_points;}


      void createNeighbourPoints();
      int extVal(double number, double percent);
      double makeBetweenZeroAndOne(double number);      

      std::array <double, 10>  multiPolybasis(Vector3D neigbour);
      Matrix10 multiTensorProduct(std::array <double, 10> poly);
      Matrix10 multiShape();

      Matrix3 derivative(Vector3D (*f)(Vector3D vec));
      Matrix3 derivative(multiVariableFunction f);

      Matrix3_1 divergenceMatrix3_3(Matrix3 (*f)(Vector3D vec));

      void printMatrix3(Matrix3 matrix);
      void printDivergence(Matrix3_1 divergence);

      Vector3D example(Vector3D vec);      


//****************************************** 3D Bond Derivatives **************************


      std::array <double, 9>  bondMultiPolybasis(Vector3D neigbour);
      Matrix9 bondMultiTensorProduct(std::array <double, 9> poly);
      Matrix9 bondMultiShape();

      Matrix3 bondDerivative(Vector3D (*f)(Vector3D vec));
      Matrix3 bondDerivative(multiVariableFunction f);

      Matrix3_1 bondDivergenceMatrix3_3(Matrix3 (*f)(Vector3D vec));
      Matrix3_1 forceStateMatrix3_3(Matrix3 (*f)(Vector3D vec));


//**********************************************************************************************************************************************************           
//      double derivative(scalar f, double point, double dx, int numOfPoints);
//      double secondDerivative(scalar f, double point, double dx, int numOfPoints);
//      double thirdDerivative(scalar f, double point, double dx, int numOfPoints);  


      double Function(double x);

//      void derivativeOfFunction(double point, double dx, int numOfPoints);  
      
    


    private:
//      scalar d_function;

      double d_horizon_size;
       
      Vector3D d_point;

      Vector3D d_vec_1;
      Vector3D d_vec_2;
      Vector3D d_vec_3;

      double d_grid_size_vec1_direction;
      double d_grid_size_vec2_direction;
      double d_grid_size_vec3_direction;
      
      double d_RHS_percent_vec1_direction;
      double d_RHS_percent_vec2_direction;
      double d_RHS_percent_vec3_direction;

      double d_LHS_percent_vec1_direction;
      double d_LHS_percent_vec2_direction;
      double d_LHS_percent_vec3_direction;


      std::vector <Vector3D>  d_neighbour_points;

  };

}

#endif

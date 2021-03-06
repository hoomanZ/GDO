#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <vector>
#include <utility>
//#include <armadillo-4.500.0/include/armadillo>
//#include<armadillo-4.450.4/include/armadillo>
#include <armadillo>
//#include<home/hzar193/libraries/armadillo-4.450.4/include/armadillo>

namespace FiniteElement {

  template<class T, int rows, int cols>
  class Matrix
  {
    public:
      Matrix();
      Matrix(const T& initialValue);
          
      const T& get(int row, int column) const;
      T& get(int row, int column);
          
      void set(int row, int column, const T& value);
          
      inline int numRows() const {return d_num_rows;}
      inline int numColumns() const {return d_num_columns;}

      Matrix inverse(const Matrix& mat);

    private:
      int d_num_rows;
      int d_num_columns;
      std::vector<std::vector<T> > d_data;

  }; // end class

  static const int DIM = 3;  //dimension
  static const int SHAPESIZE = 27;  // Max allowed search area (GIMP)

  typedef Matrix<double, 1, DIM>  MatrixVec;
  typedef std::vector<MatrixVec> ArrayMatrixVec;
  typedef Matrix<double, DIM, DIM> Matrix3;
  typedef std::vector<Matrix3> ArrayMatrix;

  typedef Matrix<int, 1, SHAPESIZE>  IntMatrixVecShape;
  typedef std::vector<IntMatrixVecShape> ArrayIntMatrixVecShape;

  typedef Matrix<double, 1, SHAPESIZE>  MatrixVecShape;
  typedef std::vector<MatrixVecShape> ArrayMatrixVecShape;
  typedef Matrix<double, SHAPESIZE, DIM> MatrixShape;
  typedef std::vector<MatrixShape> ArrayMatrixShape;

  typedef std::pair<int, ArrayMatrixVec> MatArrayMatrixVec;
  typedef std::pair<int, ArrayMatrix> MatArrayMatrix;
  typedef std::pair<int, ArrayIntMatrixVecShape> MatArrayIntMatrixVecShape;
  typedef std::pair<int, ArrayMatrixVecShape> MatArrayMatrixVecShape;
  typedef std::pair<int, ArrayMatrixShape> MatArrayMatrixShape;

} // end namespace
              





































#endif

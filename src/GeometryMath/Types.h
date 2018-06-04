#ifndef MATITI_TYPES_H
#define MATITI_TYPES_H

#include <cstdint>
#include <array>
//#include <armadillo-4.500.0/include/armadillo>
//#include <armadillo-4.450.4/include/armadillo>
////#include <armadillo>
#include <eigen/Eigen/Dense>
//#include<armadillo-4.450.4/include/armadillo>
//#include<home/hzar193/libraries/armadillo-4.450.4/include/armadillo>

using namespace Eigen;

namespace FiniteElement {
  
  typedef uint8_t u8;
  typedef uint32_t u32;
  typedef uint64_t u64;
  typedef uint64_t u128;
  typedef int64_t long64;

  typedef std::array<int, 3> IntArray3;
  typedef std::array<double, 3> Array3;

  typedef Matrix<double, 24, 24> Matrix24;
  typedef Matrix<double, 24, 1> Vector24;

  typedef Matrix<double, 2, 2> Matrix2;

  typedef Matrix<double, 3, 1> Matrix3_1; 
  typedef Matrix<double, 3, 3> Matrix3;
  typedef Matrix<double, 1, 3> Matrix1_3;

  typedef Matrix<double, 9, 9> Matrix9;
  typedef Matrix<double, 9, 3> Matrix9_3;
  typedef Matrix<double, 9, 1> Matrix9_1;
  typedef Matrix<double, 3, 9> Matrix3_9;
  typedef Matrix<double, 1, 9> Matrix1_9;

  typedef Matrix<double, 10, 10> Matrix10;
  typedef Matrix<double, 10, 3> Matrix10_3;
  typedef Matrix<double, 10, 1> Matrix10_1;
  typedef Matrix<double, 1, 10> Matrix1_10;
  typedef Matrix<double, 3, 10> Matrix3_10;

  typedef Matrix<double, 4, 4> Matrix4_4;
  typedef Matrix<double, 1, 4> Matrix1_4;
  typedef Matrix<double, 4, 1> Matrix4_1;

//  typedef arma::mat::fixed<24, 24> Matrix24;
//  typedef arma::vec::fixed<24> Vector24;

  typedef Matrix<double, Dynamic, Dynamic> MatrixXd;
  typedef Matrix<double, Dynamic, 1> VectorXd;


  typedef union 
  {
    long long i64;
    double d64;
  } dbl_64;

}

#endif

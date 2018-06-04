#ifndef TENSOR_H_
#define TENSOR_H_

#include <GeometryMath/Vector3D.h>
#include <GeometryMath/Matrix3D.h>

#include <cstdint>
#include <array>
 
namespace FiniteElement{
class Tensor
   {

//   typedef std::array <double, 3>  Array3;
//   typedef std::array <Vector3D, 3>  Matrix3;
   public:
      Tensor();
      Tensor(double value);
      Tensor(const Tensor& t);
      Tensor(const Vector3D& v1, const Vector3D& v2,
             const Vector3D& v3, const Vector3D& v4);
      Tensor(const Matrix3D& mat00,
             const Matrix3D& mat01,
             const Matrix3D& mat02,
             const Matrix3D& mat10,
             const Matrix3D& mat11,
             const Matrix3D& mat12,
             const Matrix3D& mat20,
             const Matrix3D& mat21,
             const Matrix3D& mat22);

      ~Tensor();

      void set(double value);
      void set(int i, int j, int k, int l, double value);
      double get(int i, int j, int k, int l);
      void operator = (const Tensor& t);
      void operator * (const double value);
      void operator / (const double value);
      Tensor operator *= (const double value) const;
      Tensor operator /= (const double value) const;
      double operator () (int i, int j, int k, int l) const;
      double & operator () (int i, int j, int k, int l);
      Tensor operator + (const Tensor& t) const;
      Tensor operator - (const Tensor& t) const;
      Tensor operator * (const Tensor& t) const;
//      Tensor operator * (double c, Tensor& t) const;
//      std::array <std::array<double,3>, 3> doubleContract(const Matrix3D& m) const;
      Matrix3D doubleContract(const Matrix3D& m) const;
      double Contract(const Tensor& t1) const;
      bool operator == (const Tensor& t) const;
      bool operator != (const Tensor& t) const
               {return !(*this == t); }

   private:
      double tensor[3][3][3][3];
   };

  












}
#endif

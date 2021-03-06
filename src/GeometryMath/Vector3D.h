#ifndef __MATITI_VECTOR3D_H__
#define __MATITI_VECTOR3D_H__

#include <GeometryMath/Types.h>
#include <GeometryMath/Point3D.h>
#include <GeometryMath/IntVector3D.h>
#include <GeometryMath/Matrix3D.h>
#include <cmath>

namespace FiniteElement {

  class Vector3D {

  public:

    friend std::ostream& operator<<(std::ostream& os, const FiniteElement::Vector3D& p);

  public:

    Vector3D(): d_vec({{std::numeric_limits<double>::max(), 
                        std::numeric_limits<double>::max(), 
                        std::numeric_limits<double>::max()}}) {}
    Vector3D(double val) : d_vec({{val, val, val}}) {}
    Vector3D(double x, double y, double z): d_vec({{x, y, z}}) {}
    Vector3D(const Vector3D& vec);
    Vector3D(const Point3D& start, const Point3D& end);

    ~Vector3D() {}

    //  **WARNING** Not checking index.  Do checking outside.
    double operator[](int index) const {return d_vec[index];}
    double& operator[](int index) {return d_vec[index];}

    void set(const double& val) {d_vec[0] = val; d_vec[1] = val; d_vec[2] = val;}
    void x(const double xx) {d_vec[0] = xx;}
    double x() const {return d_vec[0];}
    void y(const double yy) {d_vec[1] = yy;}
    double y() const {return d_vec[1];}
    void z(const double zz) {d_vec[2] = zz;}
    double z() const {return d_vec[2];}
    
    double length() const;
    double lengthSq() const;
    Vector3D normalized() const;

    Vector3D invDirection() const;

    double dot(const Vector3D& vec) const;
    Vector3D cross(const Vector3D& vec) const;

    bool operator==(const Vector3D& vec) const;
    bool operator!=(const Vector3D& vec) const;

    void operator=(const Vector3D& vec);

    Vector3D operator*(const double fac) const;
    Vector3D operator/(const double fac) const;
    Vector3D operator*(const Vector3D& vec) const;
    Vector3D operator/(const Vector3D& vec) const;

    Vector3D operator+(const Vector3D& vec) const;
    Vector3D operator-(const Vector3D& vec) const;

    Vector3D& operator*=(const double);
    Vector3D& operator/=(const double);
    Vector3D& operator+=(const Vector3D& vec);
    Vector3D& operator-=(const Vector3D& vec);

//    Matrix3D tensProduct(const Vector3D& vec);

    IntVector3D floor() const;
    bool operator>(double value) const {
      return (d_vec[0] > value && d_vec[1] > value && d_vec[2] > value);
    }
    bool operator<(double value) const {
      return (d_vec[0] < value && d_vec[1] < value && d_vec[2] < value);
    }

    void reset() {d_vec[0] = 0.0; d_vec[1] = 0.0; d_vec[2] = 0.0;}
    bool isnan() const {return std::isnan(d_vec[0]) || std::isnan(d_vec[1]) || std::isnan(d_vec[2]);}

    double max() const;
    double min() const;

  private:

    //double d_x, d_y, d_z;
    Array3 d_vec;

  }; // end class

} // end namespace

namespace FiniteElement {

 Vector3D min(const Vector3D& v1, const Vector3D& v2);
 Vector3D max(const Vector3D& v1, const Vector3D& v2);

} // end namespace


#endif


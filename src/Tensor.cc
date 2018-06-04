#include "Tensor.h"



#include <iostream>
#include <iomanip>

using namespace std;
using namespace FiniteElement;

Tensor::Tensor()
{
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          tensor[i][j][k][l]=0.0;
        }
      }
    }
  }
}


Tensor::Tensor(double value)
{
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          tensor[i][j][k][l]=value;
        }
      }
    }
  }
}


Tensor::Tensor(const Tensor& t)
{ 
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          tensor[i][j][k][l]=t.tensor[i][j][k][l];
        }
      }
    }
  }
}  
  

Tensor::Tensor(const Vector3D& v1, const Vector3D& v2,
               const Vector3D& v3, const Vector3D& v4)
{
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          tensor[i][j][k][l]=v1[i]*v2[j]*v3[k]*v4[l];
        }
      }
    }
  }
}  


Tensor::Tensor(const Matrix3D& mat00,
               const Matrix3D& mat01,
               const Matrix3D& mat02,
               const Matrix3D& mat10,
               const Matrix3D& mat11,
               const Matrix3D& mat12,
               const Matrix3D& mat20,
               const Matrix3D& mat21,
               const Matrix3D& mat22)
{
  for (int k=0;k<3;k++){
    for (int l=0;l<3;l++){
      tensor[0][0][k][l]=mat00(k,l);
    }
  }
  for (int k=0;k<3;k++){
    for (int l=0;l<3;l++){
      tensor[0][1][k][l]=mat01(k,l);
    }
  }
  for (int k=0;k<3;k++){
    for (int l=0;l<3;l++){
      tensor[0][2][k][l]=mat02(k,l);
    }
  }
  for (int k=0;k<3;k++){
    for (int l=0;l<3;l++){
      tensor[1][0][k][l]=mat10(k,l);
    }
  }
  for (int k=0;k<3;k++){
    for (int l=0;l<3;l++){
      tensor[1][1][k][l]=mat11(k,l);
    }
  }
  for (int k=0;k<3;k++){
    for (int l=0;l<3;l++){
      tensor[1][2][k][l]=mat12(k,l);
    }
  }
  for (int k=0;k<3;k++){
    for (int l=0;l<3;l++){
      tensor[2][0][k][l]=mat20(k,l);
    }
  }
  for (int k=0;k<3;k++){
    for (int l=0;l<3;l++){
      tensor[2][1][k][l]=mat21(k,l);
    }
  }
  for (int k=0;k<3;k++){
    for (int l=0;l<3;l++){
      tensor[2][2][k][l]=mat22(k,l);
    }
  }
}


Tensor::~Tensor()
{
}


void Tensor::set(double value)
{
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          tensor[i][j][k][l]=value;
        }
      }
    }
  }
}


void Tensor::set(int i, int j, int k, int l, double value)
{
  tensor[i][j][k][l]=value;
}


double Tensor::get(int i, int j, int k, int l)
{
  return tensor[i][j][k][l];
}

void Tensor::operator = (const Tensor& t)
{ 
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          tensor[i][j][k][l]=t.tensor[i][j][k][l];
        }
      }
    }
  }
}

void Tensor::operator * (const double value)
{ 
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          tensor[i][j][k][l] *= value;
        }
      }
    }
  }
}

void Tensor::operator / (const double value)
{ 
  double ivalue = 1./value;
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          tensor[i][j][k][l] *= ivalue;
        }
      }
    }
  }
}

Tensor Tensor::operator *= (const double value) const
{
  Tensor t; 
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          t.tensor[i][j][k][l] = value*tensor[i][j][k][l];
        }
      }
    }
  }
  return t;
}        
  
Tensor Tensor::operator /= (const double value) const
{
  double ivalue = 1./value;
  Tensor t; 
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          t.tensor[i][j][k][l] = ivalue*tensor[i][j][k][l];
        }
      }
    }
  }
  return t;
}        

double Tensor::operator () (int i, int j, int k, int l) const
{
  return tensor[i][j][k][l];
}

double &Tensor::operator () (int i, int j, int k, int l)
{
  return tensor[i][j][k][l];
}

Tensor Tensor::operator + (const Tensor& t) const
{
  Tensor ten;
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          ten.tensor[i][j][k][l] = tensor[i][j][k][l] + t.tensor[i][j][k][l];
        }
      }
    }
  }
  return ten;
}

Tensor Tensor::operator - (const Tensor& t) const
{
  Tensor ten;
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          ten.tensor[i][j][k][l] = tensor[i][j][k][l] - t.tensor[i][j][k][l];
        }
      }
    }
  }
  return ten;
} 

Tensor Tensor::operator * (const Tensor& t) const
{
  Tensor ten;
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          ten.tensor[i][j][k][l] = tensor[i][j][k][l] * t.tensor[i][j][k][l];
        }
      }
    }
  }
  return ten;
} 

//std::array <std::array<double,3>, 3> Tensor::doubleContract(const Matrix3D& m) const
Matrix3D Tensor::doubleContract(const Matrix3D& m) const
{
  double sum = 0.0;
  Matrix3D mat3;
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      sum = 0.0;
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          sum += tensor[i][j][k][l]*m(k,l);
        }
      }
      mat3(i,j) = sum;
    }
  }
  return mat3;
}

double Tensor::Contract(const Tensor& t1) const
{
  double sum = 0.0;
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          sum += tensor[i][j][k][l]*t1.tensor[i][j][k][l];
        }
      }
    }
  }
  return sum;
}                     




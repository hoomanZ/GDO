#include <StateBond.h>
#include <StateMaterialPoint.h>

#include <math.h>

using namespace FiniteElement;


StateBond::StateBond()
  : d_center_point(new StateMaterialPoint()), 
    d_second_point(new StateMaterialPoint()),
    d_influence_function(1.0),
    d_broken_bond(false)
{
  Vector3D vec(1.0, 1.0, 1.0);
  setForceVectorState(vec);
  setFirstPiolaDivergenceIntegrand(vec);
  setSai(getSecondPoint()->getPosOld() - getCenterPoint()->getPosOld());
  setDeformationVector(getSecondPoint()->getPosNew() - getCenterPoint()->getPosNew());
//  setEta(getSecondPoint()->getDisp() - getCenterPoint()->getDisp());
  fillPolynomialBasisOrder2();
}


StateBond::StateBond(StateMaterialPointP center, StateMaterialPointP second, bool broken)
  : d_center_point(center), d_second_point(second), 
    d_influence_function(1.0), d_broken_bond(false)
{
  Vector3D vec(1.0, 1.0, 1.0);
  setForceVectorState(vec);
  setFirstPiolaDivergenceIntegrand(vec);
  setSai(getSecondPoint()->getPosOld() - getCenterPoint()->getPosOld());
  setDeformationVector(getSecondPoint()->getPosNew() - getCenterPoint()->getPosNew());
//  setEta(getSecondPoint()->getDisp() - getCenterPoint()->getDisp());
  fillPolynomialBasisOrder2();
}


StateBond::StateBond(StateMaterialPointP center, StateMaterialPointP second,
                     const double& horizon, bool broken)
  : d_center_point(center), d_second_point(second), 
    d_broken_bond(false)
{
  Vector3D vec(1.0, 1.0, 1.0);
  setForceVectorState(vec);
  setFirstPiolaDivergenceIntegrand(vec);
  setSai(getSecondPoint()->getPosOld() - getCenterPoint()->getPosOld());
  setDeformationVector(getSecondPoint()->getPosNew() - getCenterPoint()->getPosNew());
  fillPolynomialBasisOrder2();
//  setEta(getSecondPoint()->getDisp() - getCenterPoint()->getDisp());
//  findInfluenceFunction(horizon);
  normalExponentialInfluenceFunction(horizon);

 
}



void
StateBond::clearBond()
{
  StateMaterialPointP center(new StateMaterialPoint());
  StateMaterialPointP second(new StateMaterialPoint());
  setCenterPoint(center);  
  setSecondPoint(second);
  setBrokenBond(false);
  Vector3D vec(1.0, 1.0, 1.0);
  setForceVectorState(0);
  setFirstPiolaDivergenceIntegrand(vec);
  setInfluenceFunction(1.0);
  setSai(getSecondPoint()->getPosOld() - getCenterPoint()->getPosOld());
  setDeformationVector(getSecondPoint()->getPosNew() - getCenterPoint()->getPosNew());
//  setEta(getSecondPoint()->getDisp() - getCenterPoint()->getDisp());

}



double
StateBond::lengthOfSai()
{
  return (d_second_point->getPosOld() - d_center_point->getPosOld()).length();
}



void
StateBond::calculateVolume(double horizon, double rj, double delX)
{
  Vector3D sai = getSecondPoint()->getPosOld() - getCenterPoint()->getPosOld();
  double distance = sai.length();  
  double volume = delX*delX*delX;
//  double volume = getSecondPoint()->getVolume();
  double coefficient = 0.5 + (horizon - distance)/(2*rj);
  if (((horizon - rj) <= distance) && (distance <= horizon))
  {
    volume *= coefficient;
  }
  else if (distance <= (horizon - rj))
  {
    volume *= 1;
  }
  else 
  {
    volume = 0.0;
  }

  setVolume(volume);

}



void
StateBond::averageStrainTensor()
{
/*  getCenterPoint()->shapeTensor();
  getCenterPoint()->deformationGradient();
  getCenterPoint()->LagrangianStrainTensor();

  getSecondPoint()->shapeTensor();
  getSecondPoint()->deformationGradient();
  getSecondPoint()->LagrangianStrainTensor();*/

//  if (getCenterPoint()->getID() == 700)
//    printScreenMatrix3D(getCenterPoint()->getShapeTensor(), getCenterPoint()->getID(), "shapeTensor");

  
  Matrix3D secondStrain = getSecondPoint()->getLagrangianStrainTensor();
  Matrix3D centerStrain = getCenterPoint()->getLagrangianStrainTensor();
  setAverageStrainTensor((secondStrain + centerStrain)*(0.5));

}



void
StateBond::equivalentStrain()
{
  averageStrainTensor();
  Matrix3D E = getAverageStrainTensor();
  double diag = (E(0,0) - E(1,1))*(E(0,0) - E(1,1));
  diag = diag + (E(2,2) - E(1,1))*(E(2,2) - E(1,1));
  diag = diag + (E(2,2) - E(0,0))*(E(2,2) - E(0,0));
  diag = diag*(2/9);
  
  double antiDiag = E(0,1)*E(0,1);
  antiDiag = antiDiag + E(0,2)*E(0,2); 
  antiDiag = antiDiag + E(1,2)*E(1,2);
  antiDiag = antiDiag*(4/3);
  
  setEquivalentStrain(pow((diag + antiDiag), 0.5));

}



void
StateBond::volumetricStrain()
{
  Matrix3D E = getAverageStrainTensor();
  double volumeStrain = E(1,1) + E(2,2) + E(0,0)
                      + E(1,1)*E(2,2) + E(2,2)*E(0,0)
                      + E(0,0)*E(1,1) + E(1,1)*E(2,2)*E(0,0);
  setVolumetricStrain(volumeStrain);

}



void
StateBond::printScreenMatrix3D(const Matrix3D& matrix, const int& id, const std::string& name)
{
  std::cout << "The " << name << " of point number " << id << ":" << std::endl;
  std::cout << matrix(0,0) << " " << matrix(0,1) << " " << matrix(0,2) << std::endl;
  std::cout << matrix(1,0) << " " << matrix(1,1) << " " << matrix(1,2) << std::endl;
  std::cout << matrix(2,0) << " " << matrix(2,1) << " " << matrix(2,2) << std::endl << std::endl;

}



Matrix3D
StateBond::tensorProduct(const Vector3D& vec1, const Vector3D& vec2)
{
  Matrix3D tensor(0.0);
  int i = 0;
  int j = 0;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      tensor(i,j) = vec1[i]*vec2[j];   
 
    }

  }

  return tensor;

}


Matrix3D
StateBond::saiTensorSai()
{
  setSai(getSecondPoint()->getPosOld() - getCenterPoint()->getPosOld());
  return tensorProduct(getSai(), getSai());
}


Matrix3D
StateBond::deformVecTensorSai()
{
  setSai(getSecondPoint()->getPosOld() - getCenterPoint()->getPosOld());
  setDeformationVector(getSecondPoint()->getPosNew() - getCenterPoint()->getPosNew());
  return tensorProduct(getDeformationVector(), getSai());
}


Matrix3D
StateBond::velocityVecTensorSai()
{
  setSai(getSecondPoint()->getPosOld() - getCenterPoint()->getPosOld());
  setVelocityVector(getSecondPoint()->getVelNew() - getCenterPoint()->getVelNew());
  return tensorProduct(getVelocityVector(), getSai());
}



void
StateBond::findInfluenceFunction(const double& horizon)
{
  double fraction = (getSai().length())/(horizon);
  double influence = pow(fraction, 10);
  influence *= 100;
  setInfluenceFunction(influence);

}



void
StateBond::normalExponentialInfluenceFunction(const double& horizon)
{
  double power = ((-1)*getSai().length()*getSai().length())/(horizon*horizon);
  setInfluenceFunction(exp(power));

//  std::cout << getInfluenceFunction() << std::endl;


}


void
StateBond::calculateEta()
{
  setEta(getSecondPoint()->getDisp() - getCenterPoint()->getDisp());
}  



Vector3D
StateBond::hourglassForceDensity(const double& springConstant, const Matrix3D& deformation)
{
  Point3D zero(0.0, 0.0, 0.0);
  Vector3D x_direct = (getCenterPoint()->getPosNew() - zero) + deformation*getSai();
  Vector3D h = x_direct - (getSecondPoint()->getPosNew() - zero);
  Vector3D Y = getSecondPoint()->getPosNew() - getCenterPoint()->getPosNew();
  double h_project = h.dot(Y);
  double lengthSai = getSai().length();
  double lengthY = Y.length();
  Vector3D unitY = Y/(lengthY);
  return (unitY*((h_project*springConstant)/(lengthSai)));

}


void
StateBond::fillPolynomialBasisOrder1()
{
  Vector3D  polynomial(0.0, 0.0, 0.0);
  setSai(getSecondPoint()->getPosOld() - getCenterPoint()->getPosOld()); 
  polynomial.x(getSai().x());
  polynomial.y(getSai().y());
  polynomial.z(getSai().z());
  setPolynomialBasisOrder1(polynomial);

}



void
StateBond::fillPolynomialBasisOrder2()
{
  setSai(getSecondPoint()->getPosOld() - getCenterPoint()->getPosOld()); 
  std::array <double, 9>  polynomial{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  std::fill_n(polynomial.begin(), 1, getSai().x());
  std::fill_n(polynomial.begin()+1, 1, getSai().y());
  std::fill_n(polynomial.begin()+2, 1, getSai().z());  
  std::fill_n(polynomial.begin()+3, 1, getSai().x()*getSai().x());
  std::fill_n(polynomial.begin()+4, 1, getSai().y()*getSai().y());
  std::fill_n(polynomial.begin()+5, 1, getSai().z()*getSai().z());
  std::fill_n(polynomial.begin()+6, 1, getSai().x()*getSai().y());
  std::fill_n(polynomial.begin()+7, 1, getSai().x()*getSai().z());
  std::fill_n(polynomial.begin()+8, 1, getSai().y()*getSai().z());

  setPolynomialBasisOrder2(polynomial);

}



Matrix9
StateBond::polynomialsOrder2TensorPolynomialsOrder2()
{
  setSai(getSecondPoint()->getPosOld() - getCenterPoint()->getPosOld());  
  fillPolynomialBasisOrder2();
  std::array <double, 9> polynomial = getPolynomialBasisOrder2();
  Matrix9 tensor9;
  for (int i = 0; i < 9; i++)
  {
     for (int j = 0; j < 9; j++)
     {
        tensor9(i, j) = polynomial[i]*polynomial[j];

     }

  }

  return tensor9;

}


Matrix3_9
StateBond::deformVecTensorPolynomialsOrder2()
{
//  setSai(getSecondPoint()->getPosOld() - getCenterPoint()->getPosOld());
  setDeformationVector(getSecondPoint()->getPosNew() - getCenterPoint()->getPosNew());  
  fillPolynomialBasisOrder2();
  std::array <double, 9> polynomial = getPolynomialBasisOrder2();
  Matrix3_9 tensor3_9;
  for (int i = 0; i < 3; i++)
  {
     for (int j = 0; j < 9; j++)
     {
        tensor3_9(i, j) = getDeformationVector()[i]*polynomial[j];

     }

  }

  return tensor3_9;
  

}



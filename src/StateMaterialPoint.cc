#include <StateMaterialPoint.h>
#include <PeriBondP.h>
#include <PeriBond.h>

#include <fstream>

#include <StateBondPArray.h>
#include <StateBondP.h>
#include <StateBond.h>

#include <math.h>





using namespace FiniteElement;
using namespace Eigen;


StateMaterialPoint::StateMaterialPoint()
  : PeriMaterialPoint::PeriMaterialPoint(), d_shape_tensor(0.0),
    /*d_deformation_gradient(0.0),*/ d_velocity_gradient(0.0),
    d_deformation_rate_tensor(0.0), d_spin_tensor(0.0),
//    d_rotation_tensor(0.0), 
    d_accele_new(0.0),
    d_residual_vector(0.0),
    d_disp_new(0.0),
    d_angular_velocity(0.0, 0.0, 0.0),
    d_unrotated_defromation_rate(0.0), d_total_elastic_strian_increment(0.0),
    d_deviatoric_strian_increment(0.0), d_unrotated_cauchy_stress(0.0),
    d_rotated_cauchy_stress(0.0), d_first_Piola_Kirchhoff_stress(0.0), d_order_flag(0)
//    d_shape_tensor_polynomial_order_2::Zero()

{
  d_deformation_gradient.Identity();
  d_rotation_tensor.Identity();
  d_right_stretch.Identity();
  d_left_stretch.Identity();
  d_Broyden_matrix.Identity(); 
  d_state_bonds_array.reserve(6000);

////  deformationGradientUsingPolynomials();

//  infinitesimalStrainTensor();  

}



StateMaterialPoint::StateMaterialPoint(PeriMaterialPointP point)
  : /*PeriMaterialPoint::PeriMaterialPoint(),*/ d_shape_tensor(0.0),
    /*d_deformation_gradient(0.0),*/ d_velocity_gradient(0.0),
    d_deformation_rate_tensor(0.0), d_spin_tensor(0.0),
    d_accele_new(0.0),
    d_residual_vector(0.0),
    d_disp_new(0.0),
//    d_rotation_tensor(0.0),
    d_angular_velocity(0.0, 0.0, 0.0),
    d_unrotated_defromation_rate(0.0), d_total_elastic_strian_increment(0.0),
    d_deviatoric_strian_increment(0.0), d_unrotated_cauchy_stress(0.0),
    d_rotated_cauchy_stress(0.0), d_first_Piola_Kirchhoff_stress(0.0), d_order_flag(0)

{
  d_deformation_gradient.Identity();
  d_rotation_tensor.Identity();
  d_right_stretch.Identity();
  d_left_stretch.Identity();
  d_Broyden_matrix.Identity(); 
  d_state_bonds_array.reserve(6000);
  setVelMid(point->getVelMid());
  setDispOld(point->getDispOld());
  setDispNew(point->getDispNew());
  setAccele(point->getAccele());
  setHorizonBonds(point->getHorizonBonds());
  setLocalDamage(point->getLocalDamage());
  setUniaxialTensionX(point->getUniaxialTensionX());
  setUniaxialTensionY(point->getUniaxialTensionY());
  setUniaxialTensionZ(point->getUniaxialTensionZ());
  setVelocityBoundaryFlag(point->getVelocityBoundaryFlag());
  setPeriStress(point->getPeriStress());
  setPeriStrain(point->getPeriStrain());

   PeriMaterialPoint::setPosOld(point->getPosOld()); 
   PeriMaterialPoint::setPosNew(point->getPosNew());
   PeriMaterialPoint::setDisp(point->getDisp());
   PeriMaterialPoint::setVelOld(point->getVelOld());
   PeriMaterialPoint::setVelNew(point->getVelNew());
   PeriMaterialPoint::setBodyForce(point->getBodyForce());
   PeriMaterialPoint::setPointForce(point->getPointForce());
   PeriMaterialPoint::setDeformationGradientOld(point->getDeformationGradientOld());
   PeriMaterialPoint::setDeformationGradientNew(point->getDeformationGradientNew());
   PeriMaterialPoint::setStrainOld(point->getStrainOld());
   PeriMaterialPoint::setStrainNew(point->getStrainNew());
   PeriMaterialPoint::setStressOld(point->getStressOld());
   PeriMaterialPoint::setStressNew(point->getStressNew());
   PeriMaterialPoint::setStrainIncrement(point->getStrainIncrement());
   PeriMaterialPoint::setStressIncrement(point->getStressIncrement());
   PeriMaterialPoint::setID(point->getID());
   PeriMaterialPoint::setMass(point->getMass());
   PeriMaterialPoint::setVolume(point->getVolume());
   PeriMaterialPoint::setInitialVolume(point->getInitialVolume());
   PeriMaterialPoint::setNoneZeroExtForce(point->getNoneZeroExtForce());
   PeriMaterialPoint::setHasBoundaryCondition(point->getHasBoundaryCondition());

////  deformationGradientUsingPolynomials();

//  infinitesimalStrainTensor(); 

}



StateMaterialPoint::StateMaterialPoint(PeriMaterialPointP point, const Matrix3D& InitialVelocityGradient)
  : /*PeriMaterialPoint::PeriMaterialPoint(),*/ d_shape_tensor(0.0),
    /*d_deformation_gradient(0.0),*/ d_velocity_gradient(InitialVelocityGradient),
    d_deformation_rate_tensor(0.0), d_spin_tensor(0.0),
    d_accele_new(0.0),
    d_residual_vector(0.0),
    d_disp_new(0.0),
//    d_rotation_tensor(0.0),
    d_angular_velocity(0.0, 0.0, 0.0),
    d_unrotated_defromation_rate(0.0), d_total_elastic_strian_increment(0.0),
    d_deviatoric_strian_increment(0.0), d_unrotated_cauchy_stress(0.0),
    d_rotated_cauchy_stress(0.0), d_first_Piola_Kirchhoff_stress(0.0), d_order_flag(0)

{
  d_deformation_gradient.Identity();
  d_rotation_tensor.Identity();
  d_right_stretch.Identity();
  d_left_stretch.Identity();
  d_Broyden_matrix.Identity(); 
  d_state_bonds_array.reserve(6000);
  setVelMid(point->getVelMid());
  setDispOld(point->getDispOld());
  setDispNew(point->getDispNew());
  setAccele(point->getAccele());
  setHorizonBonds(point->getHorizonBonds());
  setLocalDamage(point->getLocalDamage());
  setUniaxialTensionX(point->getUniaxialTensionX());
  setUniaxialTensionY(point->getUniaxialTensionY());
  setUniaxialTensionZ(point->getUniaxialTensionZ());
  setVelocityBoundaryFlag(point->getVelocityBoundaryFlag());
  setPeriStress(point->getPeriStress());
  setPeriStrain(point->getPeriStrain());

   PeriMaterialPoint::setPosOld(point->getPosOld()); 
   PeriMaterialPoint::setPosNew(point->getPosNew());
   PeriMaterialPoint::setDisp(point->getDisp());
   PeriMaterialPoint::setVelOld(point->getVelOld());
   PeriMaterialPoint::setVelNew(point->getVelNew());
   PeriMaterialPoint::setBodyForce(point->getBodyForce());
   PeriMaterialPoint::setPointForce(point->getPointForce());
   PeriMaterialPoint::setDeformationGradientOld(point->getDeformationGradientOld());
   PeriMaterialPoint::setDeformationGradientNew(point->getDeformationGradientNew());
   PeriMaterialPoint::setStrainOld(point->getStrainOld());
   PeriMaterialPoint::setStrainNew(point->getStrainNew());
   PeriMaterialPoint::setStressOld(point->getStressOld());
   PeriMaterialPoint::setStressNew(point->getStressNew());
   PeriMaterialPoint::setStrainIncrement(point->getStrainIncrement());
   PeriMaterialPoint::setStressIncrement(point->getStressIncrement());
   PeriMaterialPoint::setID(point->getID());
   PeriMaterialPoint::setMass(point->getMass());
   PeriMaterialPoint::setVolume(point->getVolume());
   PeriMaterialPoint::setInitialVolume(point->getInitialVolume());
   PeriMaterialPoint::setNoneZeroExtForce(point->getNoneZeroExtForce());
   PeriMaterialPoint::setHasBoundaryCondition(point->getHasBoundaryCondition());


////  deformationGradientUsingPolynomials();

//  infinitesimalStrainTensor(); 

}



// ******************************************************************************************************************************************************************************
// (1) **************************************************************** CalculateDeformationGradientTensorFromCurrentPosition ***************************************************


Vector3D
StateMaterialPoint::returnCurrentPosition()
{
  Vector3D position = newPositionVector();
  return position; 

}

Vector3D
StateMaterialPoint::returnCurrentVelocity()
{
  Vector3D velocity = getVelNew();
  return velocity;

}


Vector3D
StateMaterialPoint::returnDisplacement()
{
  Vector3D disp = getDisp();
  return disp; 

}


Vector3D
StateMaterialPoint::returnDisplacementNew()
{
  Vector3D dispNew = getDispNew();
  return dispNew; 

}


Vector3D
StateMaterialPoint::returnVectorExample()
{
  Vector3D vector(0.0);
  Vector3D pos = newPositionVector();
  double x = pos.x();
  double y = pos.y();
  double z = pos.z();
  vector.x(2*x*x - 5*x*y + 7*y*z*z);
  vector.y(3*x*y*z + 6*pow(x,2)*pow(y,2)*pow(z,2) - 7*pow(z,4)*pow(y,3));  
  vector.z(5*z*z + 12*x*x + 15*y*y);
  return vector;
}


void
StateMaterialPoint::derivative(const bool& secondOrder)
{
  Matrix3D derivativeExample(0.0);
  field example = &StateMaterialPoint::returnVectorExample;
  if (secondOrder == true)
  {
    std::cout << "SecondOrder Polynomials" << std::endl;
    derivativeExample = derivativeOperatorApproximationUsingSecondOrderPolynomialBasis(example);
  }
  else
  {
    std::cout << "FirstOrder Polynomials" << std::endl;
    derivativeExample = derivativeOperatorApproximationUsingFirstOrderPolynomialBasis(example);
  }
  printScreenMatrix3D(derivativeExample, "derivative of example= ");

}




Matrix3D
StateMaterialPoint::returnMatrixExample()
{
  Matrix3D matrix(0.0);
  Vector3D pos = newPositionVector();
  double x = pos.x();
  double y = pos.y();
  double z = pos.z();
  double r = pow(x*x + y*y +z*z, 0.5);

  matrix.set(0, 0, x/r);
  matrix.set(0, 1, y/r); 
  matrix.set(0, 2, z/r);
  matrix.set(1, 0, x*x*y);
  matrix.set(1, 1, x*y*z);
  matrix.set(1, 2, (-1)*x*x*y*y);
  matrix.set(2, 0, 1/r);
  matrix.set(2, 1, 1/r);
  matrix.set(2, 2, 1/r);

  return matrix;
}



/*void
StateMaterialPoint::divergance(const bool& diverganceMethod, const bool& secondOrder)
{
  if (secondOrder == true)
  {
    setOrderFlag(2);
    std::cout << "SecondOrder Polynomials" << std::endl;
  }
  else
  {
    setOrderFlag(1);
    std::cout << "FirstOrder Polynomials" << std::endl;
  }
  Vector3D diverganceExample(0.0);
  tensor example = &StateMaterialPoint::returnMatrixExample;
  if (diverganceMethod == true)
  {
    std::cout << "Divergance Method" << std::endl;
    diverganceExample = divergenceMatrix3D(example);
  }
  else
  {
    std::cout << "StateBased Peridynamics Method" << std::endl;
    diverganceExample = nonLocalDivergenceMatrix3D(example);
  }
  std::cout << "divergance of example= " << diverganceExample << std::endl;

}*/





//

void
StateMaterialPoint::shapeTensor()
{
  Matrix3D shape(0.0);
  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  int i = 0;  
//  Vector3D vecBond(0.0);
  Matrix3D tensor(0.0);
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    if (bond->getBrokenBond() == false)
    {
//      vecBond = bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      vecBond = bond->getSai();
//      Matrix3D tensor(vecBond, vecBond);
      tensor = bond->saiTensorSai();
      shape = shape + tensor*(bond->getVolume()*bond->getInfluenceFunction());
         

//      if ((bond->getCenterPoint()->getID() == 4036)) // && (bond->getSecondPoint()->getID() == 1626))
//      {
//        std::cout << "Bond(" << bond->getCenterPoint()->getID() << ", " << bond->getSecondPoint()->getID() << "):" << std::endl;
//        std::cout << bond->getSecondPoint()->getID() << std::endl;
//        std::cout << "Sai= " << bond->getSai() << std::endl;
//        std::cout << "SaiTensorSai= " << tensor << std::endl;
//        std::cout << "BondVolume= " << bond->getVolume() << std::endl;
//        std::cout << "Bond Influenve Function= " << bond->getInfluenceFunction() << std::endl;
//        std::cout << "Bond SaiTensorSai*InfFunction*bond->getInfluenceFunction()Volume= " << tensor*(bond->getVolume()*bond->getInfluenceFunction());
//        std::cout << "Integrand= " << shape << std::endl;        
//        printScreenMatrix3D(tensor, "saiTensorSai");
//        printScreenMatrix3D(shape, "shape");
//      }       
 
    }

  }
  if (shape.Determinant() == 0)
  {
    std::cout << "Shape Function of the node number " << getID() << " has zero value of determinant." << std::endl;
//    shape.Identity();
  }  

  setShapeTensor(shape);

//  std::ofstream myfile;
//  myfile.open(fileName, std::ofstream::app);
//  myfile << std::endl << "Determinant of shape matrix is= " << getShapeTensor().Determinant() << std::endl;
//  myfile.close();
//  ret = chdir(currentDirec.c_str());

//  printScreenMatrix3D(getShapeTensor());
//  std::cout << std::endl;

}



Matrix3D
StateMaterialPoint::derivativeOperatorApproximationUsingFirstOrderPolynomialBasis(field f)
{
  Vector3D fCenter = (this->*f)();
  Vector3D polynomial(0.0, 0.0, 0.0);
  shapeTensor();
  Matrix3D tensor(0.0);
  Matrix3D deformation(0.0);
  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  int i = 0;  
  Vector3D fSecond(0.0);
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    if (bond->getBrokenBond() == false)
    {
//      vecBond = bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
        bond->fillPolynomialBasisOrder1();
        polynomial = bond->getPolynomialBasisOrder1();        
        StateMaterialPointP horizonPoint = bond->getSecondPoint();
        StateMaterialPoint point = *horizonPoint;
        fSecond = (point.*f)();

        for (int ii = 0; ii < 3; ii++)
        {
          for (int jj = 0; jj < 3; jj++)
          {
            tensor(ii,jj) = (fSecond[ii] - fCenter[ii])*polynomial[jj];
          }
        }

//      tensor = bond->deformVecTensorSai();  //for simpleDeformationTensor

        deformation = deformation + tensor*(bond->getVolume()*bond->getInfluenceFunction()); 
 
    } 

  }


  if (getShapeTensor().Determinant() != 0.0)
  {
    deformation = deformation*invMatrix3D(getShapeTensor());
  }

//  if (deformation.Determinant() == 0)
//  {
//    deformation.Identity();
//  }  

  return deformation;  

}

/***************************

void
StateMaterialPoint::shapeTensorPolynomialOrder2()
{


  Matrix9 shape = Matrix9::Zero();
  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  int i = 0;  
//  Vector3D vecBond(0.0);
  Matrix9 tensor = Matrix9::Zero();
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    if (bond->getBrokenBond() == false)
    {
//      vecBond = bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      vecBond = bond->getSai();
//      Matrix3D tensor(vecBond, vecBond);
      tensor = bond->polynomialsOrder2TensorPolynomialsOrder2();
      shape = shape + tensor*(bond->getVolume()*bond->getInfluenceFunction());
//      if ((bond->getCenterPoint()->getID() == 700)) // && (bond->getSecondPoint()->getID() == 1626))
//      {
//        std::cout << bond->getSecondPoint()->getID() << std::endl;
//        std::cout << "Sai= " << bond->getSai() << std::endl;
//        std::cout << "BondVolume= " << bond->getVolume() << std::endl;
//        printScreenMatrix3D(tensor, "saiTensorSai");
//        printScreenMatrix3D(shape, "shape");
//      }       
 
    }

  }
  if (shape.determinant() == 0)
  {
    std::cout << "Shape Function of the node number " << getID() << " has zero value of determinant." << std::endl;
//    shape = Matrix9::Identity();
  }  

  setShapeTensorPolynomialOrder2(shape);
//  printScreenMatrix3D(getShapeTensor());
//  std::cout << std::endl;

}



****************/







void
StateMaterialPoint::shapeTensorPolynomialOrder2()
{

  int id = 3909;
  int id2 = 4184;

  char buffer[2000];
  char *str = getcwd(buffer, 2000);
  std::string currentDirec = std::string(buffer);
  int ret;
  std::string fileName = "shape_tensor_detail_" + std::to_string(getID()) + "_" + std::to_string(getNumIter()) + ".com";


  Matrix9 shape = Matrix9::Zero();
  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  int i = 0;  
//  Vector3D vecBond(0.0);
  Matrix9 tensor = Matrix9::Zero();
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    if (bond->getBrokenBond() == false)
    {
//      vecBond = bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      vecBond = bond->getSai();
//      Matrix3D tensor(vecBond, vecBond);
      tensor = bond->polynomialsOrder2TensorPolynomialsOrder2();
      shape = shape + tensor*(bond->getVolume()*bond->getInfluenceFunction());
      if ((bond->getCenterPoint()->getID() == id) || (bond->getCenterPoint()->getID() == id2)) // && (bond->getSecondPoint()->getID() == 1626))
      {
//        std::string fileName = "shape_tensor_detail_" + std::to_string(getID()) + "_" + std::to_string(getNumIter()) + ".com";
        std::ofstream myfile;
        myfile.open(fileName, std::ofstream::app);
 
        myfile << i << ") Bond(" << bond->getCenterPoint()->getID() << ", " << bond->getSecondPoint()->getID() << "):" << std::endl; 
        myfile << "Bond Center Point Reference Position= " << bond->getCenterPoint()->getPosOld() << std::endl;
        myfile << "Bond Second Point Reference Position= " << bond->getSecondPoint()->getPosOld() << std::endl;
//        std::cout << bond->getSecondPoint()->getID() << std::endl;
        myfile << "Sai= " << bond->getSai() << std::endl;
        myfile << "SaiTensorSai= " << std::endl << tensor << std::endl;
        myfile << "BondVolume= " << bond->getVolume() << std::endl;
        myfile << "Bond Influenve Function= " << bond->getInfluenceFunction() << std::endl;
        myfile << "Bond SaiTensorSai*InfFunction*bond->getInfluenceFunction()Volume= " << std::endl << tensor*(bond->getVolume()*bond->getInfluenceFunction()) << std::endl;
        myfile << "Integrand= " << std::endl << shape << std::endl << std::endl; 

        myfile.close();
        ret = chdir(currentDirec.c_str());
       
//        printScreenMatrix3D(tensor, "saiTensorSai");
//        printScreenMatrix3D(shape, "shape");
      }
//      if ((bond->getCenterPoint()->getID() == 700)) // && (bond->getSecondPoint()->getID() == 1626))
//      {
//        std::cout << bond->getSecondPoint()->getID() << std::endl;
//        std::cout << "Sai= " << bond->getSai() << std::endl;
//        std::cout << "BondVolume= " << bond->getVolume() << std::endl;
//        printScreenMatrix3D(tensor, "saiTensorSai");
//        printScreenMatrix3D(shape, "shape");
//      }       
 
    }

  }
  if (shape.determinant() == 0)
  {
    std::cout << "Shape Function of the node number " << getID() << " has zero value of determinant." << std::endl;
//    shape = Matrix9::Identity();
  }  

  setShapeTensorPolynomialOrder2(shape);

//  std::string fileName = "shape_tensor_detail_" + std::to_string(getID()) + "_" + std::to_string(getNumIter()) + ".com";
  if ((getID() == id) || (getID() == id2))
  {
    std::ofstream myfile;
    myfile.open(fileName, std::ofstream::app);
    myfile << "Determinant of shape matrix of " << getID() << " is= " << getShapeTensorPolynomialOrder2().determinant() << std::endl << std::endl;
    myfile.close();
    ret = chdir(currentDirec.c_str());
  }


//  printScreenMatrix3D(getShapeTensor());
//  std::cout << std::endl;

} 



Matrix3D
StateMaterialPoint::derivativeOperatorApproximationUsingSecondOrderPolynomialBasis(field f)
{
  Matrix9_3 d(9,3);

  d << 1, 0, 0,
       0, 1, 0,
       0, 0, 1,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0;

  shapeTensorPolynomialOrder2();
  Vector3D fCenter = (this->*f)();

  Matrix3D deformation(0.0);
  Matrix3 deform = Matrix3::Zero();
  Matrix3_9 sigmaTensor = Matrix3_9::Zero();
  Matrix3_9 tensorPolynomialsOrder2 = Matrix3_9::Zero();
  std::array <double, 9> polynomial{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  int i = 0;  
  Vector3D fSecond(0.0);
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    if (bond->getBrokenBond() == false)
    {
//      vecBond = bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
        bond->fillPolynomialBasisOrder2();
        polynomial = bond->getPolynomialBasisOrder2();        
        StateMaterialPointP horizonPoint = bond->getSecondPoint();
        StateMaterialPoint point = *horizonPoint;
        fSecond = (point.*f)();
//        fSecond = (horizonPoint->*)f();

//        fSecond = std::invoke(*horizonPoint, f);


        for (int ii = 0; ii < 3; ii++)
        {
          for (int jj = 0; jj < 9; jj++)
          {
            tensorPolynomialsOrder2(ii,jj) = (fSecond[ii] - fCenter[ii])*polynomial[jj];
          }
        }

//      tensor = bond->deformVecTensorSai();  //for simpleDeformationTensor

        sigmaTensor = sigmaTensor + tensorPolynomialsOrder2*(bond->getVolume()*bond->getInfluenceFunction()); 
 
    } 

  }


  if (getShapeTensorPolynomialOrder2().determinant() != 0.0)
  {
    deform = sigmaTensor*(getShapeTensorPolynomialOrder2().inverse())*d;
  }
  else
  {
    deform = sigmaTensor*(pseudoInverseMoorePenrose9by9(getShapeTensorPolynomialOrder2()))*d;
  }


/*  if (getID() == 3909)
  {
    Matrix9 m9(9, 9);
//    m9 << 1, 1, 1, 1, 1, 1, 1, 1, 1,
//          1, 1, 1, 1, 1, 1, 1, 1, 1,
//          1, 1, 1, 1, 1, 1, 1, 1, 1,
//          1, 1, 1, 1, 1, 1, 1, 1, 1,
//          1, 1, 1, 1, 1, 1, 1, 1, 1,
//          1, 1, 1, 1, 1, 1, 1, 1, 1,
//          1, 1, 1, 1, 1, 1, 1, 1, 1,
//          1, 1, 1, 1, 1, 1, 1, 1, 1,
//          1, 1, 1, 1, 1, 1, 1, 1, 1;

    m9 << 1, 1, 1, 1, 1, 1, 1, 1, 0,
          1, 1, 1, 1, 1, 1, 1, 1, 0,
          1, 1, 1, 1, 1, 1, 1, 1, 0,
          1, 1, 1, 1, 1, 1, 1, 1, 0,
          1, 1, 1, 1, 1, 1, 1, 1, 0,
          1, 1, 1, 1, 1, 1, 1, 1, 0,
          1, 1, 1, 1, 1, 1, 1, 1, 0,
          1, 1, 1, 1, 1, 1, 1, 1, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0;


    Matrix9 minv9(9, 9);

    minv9 = pseudoInverseMoorePenrose9by9(m9);

    std::cout << "PseudoInverse of m= " << std::endl << minv9 << std::endl;
    std::cout << "A*pinvA= " << std::endl << m9*minv9 << std::endl;
    std::cout << "pinvA*A= " << std::endl << minv9*m9 << std::endl;
    std::cout << "A*pinvA*A= " << std::endl << m9*minv9*m9 << std::endl;
    std::cout << "pinvA*A*pinvA= " << std::endl << minv9*m9*minv9 << std::endl;

  }*/
  

/*  if (getID() == 3909)
  {
    int a = 10;
    int b = 10;
    int multiple = 1;
    int c = a*multiple;
    int d = b*multiple;

    Matrix2 m2(2, 2);
    m2 << a, b, 
          c, d; 

    Matrix2 minv2(2, 2);

    minv2 = pseudoInverseMoorePenrose2by2(m2);

    std::cout << "PseudoInverse of m= " << std::endl << minv2 << std::endl;
    std::cout << "A*pinvA= " << std::endl << m2*minv2 << std::endl;
    std::cout << "pinvA*A= " << std::endl << minv2*m2 << std::endl;
    std::cout << "A*pinvA*A= " << std::endl << m2*minv2*m2 << std::endl;
    std::cout << "pinvA*A*pinvA= " << std::endl << minv2*m2*minv2 << std::endl;
  }*/

/*  if (getID() == 3909)
  {
    Vector3D row1(1, 2, 0);
    Vector3D row2(-3, 4, 0);
    double coeff1 = 0.0;
    double coeff2 = 0.0;
    Vector3D row3 = row1*coeff1 + row2*coeff2;
//    Vector3D row3(2, 2, 2);

    Matrix3 m3(3, 3);
    m3 << row1.x(), row1.y(), row1.z(), 
          row2.x(), row2.y(), row2.z(),
          row3.x(), row3.y(), row3.z(); 

    Matrix3 minv3(3, 3);

    minv3 = pseudoInverseMoorePenrose3by3(m3);

    int n = 0;
    double value = pow(10, n);
    std::cout << "Epsilon Machine= " << std::endl << machine_eps(value) << std::endl;

    std::cout << "PseudoInverse of m= " << std::endl << minv3 << std::endl;
    if (m3.determinant() != 0)
    {
      std::cout << "invA= " << std::endl << m3.inverse() << std::endl;
      std::cout << "A*invA= " << std::endl << m3*m3.inverse() << std::endl;
      std::cout << "invA*A= " << std::endl << m3.inverse()*m3 << std::endl;
    }
    std::cout << "A*pinvA= " << std::endl << m3*minv3 << std::endl;
    std::cout << "pinvA*A= " << std::endl << minv3*m3 << std::endl;
    std::cout << "A*pinvA*A= " << std::endl << m3*minv3*m3 << std::endl;
    std::cout << "pinvA*A*pinvA= " << std::endl << minv3*m3*minv3 << std::endl;
  }*/


  deformation = eigenMatrixToMatrix3D(deform);

//  if (deformation.Determinant() == 0)
//  {
//    deformation.Identity();
//  }  

  return deformation;

}


Matrix9
StateMaterialPoint::pseudoInverseMoorePenrose9by9(const Matrix9& matrix9)
{
//  std::cout << "Here is the matrix m= " << std::endl << matrix9 << std::endl;
  JacobiSVD<Matrix9> svd(matrix9, ComputeThinU | ComputeThinV);
//  std::cout << "Its singular values are:" << std::endl << svd.singularValues().asDiagonal() << std::endl;
//  std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svd.matrixU() << std::endl;
//  std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svd.matrixV() << std::endl;
//  std::cout << "Matrix=" << std::endl << svd.matrixU()*svd.singularValues().asDiagonal()*svd.matrixV().transpose() << std::endl;
  Matrix9 diag = svd.singularValues().asDiagonal();
  Matrix9 leftUnitary = svd.matrixU();
  Matrix9 rightUnitary = svd.matrixV(); 

  double tolerance = 0.0;
  int matrixSize = 9;
  double maxDiagonalMatrix = diag(0, 0);
  for (int j = 0; j < matrixSize; j++)
  {
    if (diag(j, j) > maxDiagonalMatrix) maxDiagonalMatrix = diag(j, j);
  } 
//  std::cout << "maxDiagonalMatrix= " << maxDiagonalMatrix << std::endl;
//  std::cout << "matrixSize= " << matrixSize << std::endl;   
  int n = 0;
  double value = pow(10, n);
  tolerance = (double) matrixSize*machine_eps(value)*maxDiagonalMatrix;
//  std::cout << "tolerance= " << tolerance << std::endl;
  for (int k = 0; k < matrixSize; k++)
  {
    if (std::abs(diag(k, k)) < tolerance) diag(k, k) = 0.0;
    for (int s = 0; s < matrixSize; s++)
    {
      if (std::abs(leftUnitary(k, s)) < tolerance) leftUnitary(k, s) = 0.0;
      if (std::abs(rightUnitary(k, s)) < tolerance) rightUnitary(k, s) = 0.0;
    }   
  } 
//  std::cout << "Its left singular vectors are the columns of the thin U matrix AFTER MODIFICATION:" << std::endl << leftUnitary << std::endl;
//  std::cout << "Its right singular vectors are the columns of the thin V matrix AFTER MODIFICATION:" << std::endl << rightUnitary << std::endl;



  Matrix9 diagonal = Matrix9::Zero(matrixSize, matrixSize);
  for (int i = 0; i < matrixSize; i++)
  {
    if (diag(i, i) != 0)
    {
      diagonal(i, i) = 1.0/diag(i, i);
    }
  }

//  std::cout << "diagonal= " << std::endl << diag << std::endl;
//  std::cout << "invDiagonal= " << std::endl << diagonal << std::endl;

  return rightUnitary*diagonal*leftUnitary.transpose();
}



Matrix2
StateMaterialPoint::pseudoInverseMoorePenrose2by2(const Matrix2& matrix2)
{
  std::cout << "Here is the matrix m= " << std::endl << matrix2 << std::endl;
  JacobiSVD<Matrix2> svd(matrix2, ComputeThinU | ComputeThinV);
//  std::cout << "Its singular values are:" << std::endl << svd.singularValues().asDiagonal() << std::endl;
  std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svd.matrixU() << std::endl;
  std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svd.matrixV() << std::endl;
  std::cout << "Matrix=" << std::endl << svd.matrixU()*svd.singularValues().asDiagonal()*svd.matrixV().transpose() << std::endl;
  Matrix2 diag = svd.singularValues().asDiagonal();
  Matrix2 leftUnitary = svd.matrixU();
  Matrix2 rightUnitary = svd.matrixV(); 

  double tolerance = 0.0;
  int matrixSize = 2;
  double maxDiagonalMatrix = diag(0, 0);
  for (int j = 0; j < matrixSize; j++)
  {
    if (diag(j, j) > maxDiagonalMatrix) maxDiagonalMatrix = diag(j, j);
  } 
  std::cout << "maxDiagonalMatrix= " << maxDiagonalMatrix << std::endl;
//  std::cout << "matrixSize= " << matrixSize << std::endl;   
  int n = 0;
  double value = pow(10, n);
  tolerance = (double) matrixSize*machine_eps(value)*maxDiagonalMatrix;
//  std::cout << "tolerance= " << tolerance << std::endl;
  for (int k = 0; k < matrixSize; k++)
  {
    if (diag(k, k) < tolerance) diag(k, k) = 0.0;
    for (int s = 0; s < matrixSize; s++)
    {
      if (std::abs(leftUnitary(k, s)) < tolerance) leftUnitary(k, s) = 0.0;
      if (std::abs(rightUnitary(k, s)) < tolerance) rightUnitary(k, s) = 0.0;
    }   
  } 
  std::cout << "Its left singular vectors are the columns of the thin U matrix AFTER MODIFICATION:" << std::endl << leftUnitary << std::endl;
  std::cout << "Its right singular vectors are the columns of the thin V matrix AFTER MODIFICATION:" << std::endl << rightUnitary << std::endl;



  Matrix2 diagonal = Matrix2::Zero(matrixSize, matrixSize);
  for (int i = 0; i < matrixSize; i++)
  {
    if (diag(i, i) != 0)
    {
      diagonal(i, i) = 1.0/diag(i, i);
    }
  }

  std::cout << "diagonal= " << std::endl << diag << std::endl;
  std::cout << "invDiagonal= " << std::endl << diagonal << std::endl;

  return rightUnitary*diagonal*leftUnitary.transpose();
}


Matrix3
StateMaterialPoint::pseudoInverseMoorePenrose3by3(const Matrix3& matrix3)
{
  std::cout << "Here is the matrix m= " << std::endl << matrix3 << std::endl;
  JacobiSVD<Matrix3> svd(matrix3, ComputeThinU | ComputeThinV);
//  std::cout << "Its singular values are:" << std::endl << svd.singularValues().asDiagonal() << std::endl;
  std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svd.matrixU() << std::endl;
  std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svd.matrixV() << std::endl;
  std::cout << "Matrix=" << std::endl << svd.matrixU()*svd.singularValues().asDiagonal()*svd.matrixV().transpose() << std::endl;
  Matrix3 diag = svd.singularValues().asDiagonal();
  Matrix3 leftUnitary = svd.matrixU();
  Matrix3 rightUnitary = svd.matrixV(); 

  double tolerance = 0.0;
  int matrixSize = 3;
  double maxDiagonalMatrix = diag(0, 0);
  for (int j = 0; j < matrixSize; j++)
  {
    if (diag(j, j) > maxDiagonalMatrix) maxDiagonalMatrix = diag(j, j);
  } 
  std::cout << "maxDiagonalMatrix= " << maxDiagonalMatrix << std::endl;
//  std::cout << "matrixSize= " << matrixSize << std::endl;   
  int n = 0;
  double value = pow(10, n);
  tolerance = (double) matrixSize*machine_eps(value)*maxDiagonalMatrix;
//  std::cout << "tolerance= " << tolerance << std::endl;
  for (int k = 0; k < matrixSize; k++)
  {
    if (diag(k, k) < tolerance) diag(k, k) = 0.0;
    for (int s = 0; s < matrixSize; s++)
    {
      if (std::abs(leftUnitary(k, s)) < tolerance) leftUnitary(k, s) = 0.0;
      if (std::abs(rightUnitary(k, s)) < tolerance) rightUnitary(k, s) = 0.0;
    }   
  } 
  std::cout << "Its left singular vectors are the columns of the thin U matrix AFTER MODIFICATION:" << std::endl << leftUnitary << std::endl;
  std::cout << "Its right singular vectors are the columns of the thin V matrix AFTER MODIFICATION:" << std::endl << rightUnitary << std::endl;



  Matrix3 diagonal = Matrix3::Zero(matrixSize, matrixSize);
  for (int i = 0; i < matrixSize; i++)
  {
    if (diag(i, i) != 0)
    {
      diagonal(i, i) = 1.0/diag(i, i);
    }
  }

  std::cout << "diagonal= " << std::endl << diag << std::endl;
  std::cout << "invDiagonal= " << std::endl << diagonal << std::endl;

  return rightUnitary*diagonal*leftUnitary.transpose();
}


double
StateMaterialPoint::machine_eps (double value)
{
    dbl_64 s;
    s.d64 = value;
    s.i64++;
    return (s.i64 < 0 ? value - s.d64 : s.d64 - value);
}



void
StateMaterialPoint::deformationGradient()
{

  shapeTensorPolynomialOrder2();

  Matrix3D deformation(0.0);

  Matrix9_3 d(9,3);
  d << 1, 0, 0,
       0, 1, 0,
       0, 0, 1,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0;

  Matrix3 deform = Matrix3::Zero();

  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  int i = 0;  
//  Vector3D vecNewBond(0.0);
//  Vector3D vecBond(0.0);
////  Matrix3D tensor(0.0); //for simpleDeformationTensor
  Matrix3_9 sigmaTensor = Matrix3_9::Zero();
  Matrix3_9 tensorPolynomialsOrder2 = Matrix3_9::Zero();
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    if (bond->getBrokenBond() == false)
    {
//      vecBond = bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      vecNewBond = bond->getSecondPoint()->getPosNew() - bond->getCenterPoint()->getPosNew();
//      Matrix3D tensor(vecNewBond, vecBond);
////      tensor = bond->deformVecTensorSai();  //for simpleDeformationTensor
      
      tensorPolynomialsOrder2 = bond->deformVecTensorPolynomialsOrder2();

////      deformation = deformation + tensor*(bond->getVolume()*bond->getInfluenceFunction());  //for simpleDeformationTensor

      sigmaTensor = sigmaTensor + tensorPolynomialsOrder2*(bond->getVolume()*bond->getInfluenceFunction());


/*      if ((bond->getCenterPoint()->getID() == 700)) // && (bond->getSecondPoint()->getID() == 1626))
      {
        std::cout << bond->getSecondPoint()->getID() << std::endl;
        std::cout << bond->getSai() << std::endl;
        std::cout << bond->getVolume() << std::endl;
//        std::cout << bond->getDeformationVector() << std::endl;
//        printScreenMatrix3D(tensor, "yaiTensorSai");
//        printScreenMatrix3D(deformation, "deformation");
        std::cout << "sigmaTensor= " << std::endl; 
        for (int ii = 0; ii < 3; ii++)
        {
           for (int jj = 0; jj < 9; jj++)
           {
             std::cout << "(" << ii << ", " << jj << ")= " << sigmaTensor(ii,jj) << " ";
           }
        }
      }*/ 
    }

//    else
//    {
//      deformation.Identity();
//    }


  }

////  deformation = deformation*getShapeTensor().Inverse();  //for simpleDeformationTensor
////  if (getShapeTensor().Determinant() == 0)  //for simpleDeformationTensor  //for simpleDeformationTensor
////  {
////    std::cout << "Hazard: Node " << getID() << " Shape Tensor has a zero determinant." << std::endl;  //for simpleDeformationTensor
////  }
////  deformation = deformation*invMatrix3D(getShapeTensor());  //for simpleDeformationTensor


/*      if (getID() == 700) // && (bond->getSecondPoint()->getID() == 1626))
      {
//        std::cout << bond->getSecondPoint()->getID() << std::endl;
//        std::cout << bond->getSai() << std::endl;
//        std::cout << bond->getVolume() << std::endl;
//        std::cout << bond->getDeformationVector() << std::endl;
//        printScreenMatrix3D(tensor, "yaiTensorSai");
//        printScreenMatrix3D(deformation, "deformation");
        std::cout << std::endl << "ShapeTensor= " << std::endl; 
        for (int ii = 0; ii < 9; ii++)
        {
           for (int jj = 0; jj < 9; jj++)
           {
             std::cout << "(" << ii << ", " << jj << ")= " << getShapeTensorPolynomialOrder2()(ii,jj) << " ";
           }
        }
      } */
   

  if (getShapeTensorPolynomialOrder2().determinant() != 0.0)
  {
    deform = sigmaTensor*(getShapeTensorPolynomialOrder2().inverse())*d;
  }

  deformation = eigenMatrixToMatrix3D(deform);


 

//  if (deformation.Determinant() == 0)
//  {
//    deformation.Identity();
//  }  


  setDeformationGradient(deformation);
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getDeformationGradient(), "DeformationGradient");
//  }

}


void
StateMaterialPoint::deformationGradientUsingPolynomials()
{

/*   if (getID() == 1055)
   {
     printScreenMatrix3D(getDeformationGradient(), "First DeformationGradient");
//     std::cout << "Displacement= " << detDisp() << std::endl;
   }


  if (getID() == 1055)
  {
    std::cout << "Old Position= " << getPosOld() << std::endl;
    std::cout << "Disp= " << getDisp() << std::endl;
    std::cout << "New Position= " << getPosNew() << std::endl << std::endl;

  }*/


//  if ((getShapeTensorPolynomialOrder2().determinant() != 0.0) || (getShapeTensor().Determinant() != 0.0))
//  { 
    Matrix3D deform(0.0);

    Matrix3D ident(0.0);
    ident.Identity();

    field currentPosition = &StateMaterialPoint::returnCurrentPosition;

    field displacement = &StateMaterialPoint::returnDisplacement;  

    if (getOrderFlag() == 1)
    {

//      deform = derivativeOperatorApproximationUsingFirstOrderPolynomialBasis(currentPosition);
      deform = derivativeOperatorApproximationUsingFirstOrderPolynomialBasis(displacement) + ident;

    }
  
    if (getOrderFlag() == 2)
    {
 
//      deform = derivativeOperatorApproximationUsingSecondOrderPolynomialBasis(currentPosition);
      deform = derivativeOperatorApproximationUsingSecondOrderPolynomialBasis(displacement) + ident;

    }

    setDeformationGradient(deform);

//     if (getID() == 1055)
//        printScreenMatrix3D(getDeformationGradient(), "Second DeformationGradient");


//  }

}

// END OF (1) *******************************************************************************************************************************************************************



// ******************************************************************************************************************************************************************************
// (2) ***************************************************** CalculateDeformationAndStrainTensorsFromDeformationGradientTensor **************************************************



void
StateMaterialPoint::rightCauchyGreenDeformation()
{

  Matrix3D F = getDeformationGradient();
  setRightCauchyGreenDeformation(F.Transpose()*F);

}


void
StateMaterialPoint::leftCauchyGreenDeformation()
{

  Matrix3D F = getDeformationGradient();
  setLeftCauchyGreenDeformation(F*F.Transpose());

}



void
StateMaterialPoint::LagrangianStrainTensor()
{
//  shapeTensor();
//  deformationGradient();

  
//  deformationGradientUsingPolynomials();

  rightCauchyGreenDeformation();

  Matrix3D C = getRightCauchyGreenDeformation();
  Matrix3D I;
  I.Identity(); 
  setLagrangianStrainTensor((C - I)*(0.5));
    
}


void
StateMaterialPoint::infinitesimalStrainTensor()
{
//  shapeTensor();
//  deformationGradient();

//  deformationGradientUsingPolynomials();

  Matrix3D F = getDeformationGradient();

//  if (getID() == 1055)
//     printScreenMatrix3D(F, "DeformationGradientTensor");


  Matrix3D I;
  I.Identity(); 
  setInfinitesimalStrainTensor((F.Transpose() + F)*(0.5) - I);
    
}


void
StateMaterialPoint::unrotatedDeformationRateTensor(const double& delT)   //  Green-Naghdi
{
  LagrangianStrainTensor();
  spinTensor();
  updateRotationTensorMoreStable(delT);
  leftStretchRateTensor();
  updateLeftStretch(delT); 
  unrotatedDeformationRate();
  strainIncrement(delT);  

}




// END OF (2) *************************************************************************************************************************************************************




// ************************************************************************************************************************************************************************
// (3) ************************************************************** CalculateStressTensorFromStrainTensor ***************************************************************



void
StateMaterialPoint::isotropicLinearElasticityCauchyStressTensor(const double& firstLamme, const double& secondLamme)
{
  Matrix3D eps = getInfinitesimalStrainTensor();
////  if (getID() == 1055)
////     printScreenMatrix3D(eps, "InfinitesimalStrainTensor");

  double epsKK = eps(0, 0) + eps(1, 1) + eps(2, 2);


  Matrix3D sigma(0.0);
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      if (i == j)
      {
        sigma(i, j) = firstLamme*epsKK + 2*secondLamme*eps(i, j);
      }
      else
      {
        sigma(i, j) = 2*secondLamme*eps(i, j);
      }
    }
  }

  setCauchyStress(sigma);

////  if (getID() == 1055)
////     printScreenMatrix3D(getCauchyStress(), "CauchyStress");


}



void
StateMaterialPoint::hypoelasticRotatedCauchyStressTensor(const double& firstLamme, const double& secondLamme)  // Green-Naghdi
{
  updateUnrotatedCauchyStress(firstLamme, secondLamme);
  rotatedCauchyStress();

}


void
StateMaterialPoint::firstPiolaKirchhoffStress(const Matrix3D& CauchyStressTensor)
{
  double J = getDeformationGradient().Determinant();
  Matrix3D Ftranspose = getDeformationGradient().Transpose();
//  if (getID() == 700)
//  {
//    std::cout << "Jacobian= " << J << std::endl;
//    printScreenMatrix3D(getDeformationGradient(), "deformationGradient");
//    printScreenMatrix3D(Ftranspose, "deformationGradientTranspose");
//  }
//  std::cout << Ftranspose(0,0) << " " << Ftranspose(0,1) << " " << Ftranspose(0,2) << std::endl;
//  std::cout << Ftranspose(1,0) << " " << Ftranspose(1,1) << " " << Ftranspose(1,2) << std::endl;
//  std::cout << Ftranspose(2,0) << " " << Ftranspose(2,1) << " " << Ftranspose(2,2) << std::endl;
//  Matrix3D FtransposeInverse = Ftranspose.Inverse();
  if (Ftranspose.Determinant() == 0)
  {
    std::cout << "Hazard 0: Node " << getID() << " Deformation gradient has a zero determinant." << std::endl;
  }
  Matrix3D FtransposeInverse = invMatrix3D(Ftranspose);
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(FtransposeInverse, "deformationGradientTransposeInverse");
//    printScreenMatrix3D(FtransposeInverse*Ftranspose, "identity");
//  }
  setFirstPiolaKirchhoffStress((CauchyStressTensor*J)*(FtransposeInverse));


}  


// END OF (3) *************************************************************************************************************************************************************



// ************************************************************************************************************************************************************************
// (4) ************************************************************************ CalculateInternalForce ********************************************************************



Vector3D
StateMaterialPoint::divergenceMatrix3D(tensor f)  // State Based GDO divergance
{
  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  int i = 0;
  double influence = 0.0;
  Matrix3_1 eigenIntegrand = Matrix3_1::Zero();
  Vector3D integrand(0.0, 0.0, 0.0);
  Vector3D sai(0.0, 0.0, 0.0);
  Vector3D zero(0.0, 0.0, 0.0);
  Matrix3D matrix(0.0);
  Matrix3D fSecond(0.0);

  Matrix3D invShapeTensor(0.0);
  Matrix9  inverseShapeTensor = Matrix9::Zero();
  Matrix9_1 polynomialEigen = Matrix9_1::Zero();
  Matrix9_1 twinPolynomialEigen = Matrix9_1::Zero();

  std::array <double, 9>  polynomial{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  
  Vector3D bondIntegrand(0.0, 0.0,0.0);

  if (getOrderFlag() == 1)
  {
    shapeTensor();
    invShapeTensor = invMatrix3D(getShapeTensor());  

  }

  if (getOrderFlag() == 2)
  {
    shapeTensorPolynomialOrder2();
    if (getShapeTensorPolynomialOrder2().determinant() != 0)
    {
      inverseShapeTensor = getShapeTensorPolynomialOrder2().inverse();
    }
    else
    {  
      inverseShapeTensor = pseudoInverseMoorePenrose9by9(getShapeTensorPolynomialOrder2());
    }
  }

  
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    if (bond->getBrokenBond() == false)
    {
      influence = bond->getInfluenceFunction();
      sai = bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
      StateMaterialPointP horizonPoint = bond->getSecondPoint();
      StateMaterialPoint point = *horizonPoint;
      fSecond = (point.*f)();
      Matrix3D fSecondMinusfFirst = fSecond - (this->*f)();

      if (getShapeTensorPolynomialOrder2().determinant() == 0)
      {
//        std::cout << "Hazard 1: Node " << getID() << " Shape Tensor has a zero determinant." << std::endl;
      }

      if (getOrderFlag() == 1)
      {       
        matrix = fSecondMinusfFirst*invShapeTensor;
        integrand = integrand + (matrix*sai)*(influence*bond->getVolume());
      }

      else if (getOrderFlag() == 2)
      {
        bond->fillPolynomialBasisOrder2();
        polynomial = bond->getPolynomialBasisOrder2();

        for (int j = 0; j < 9; j++)
        {
          polynomialEigen(j) = polynomial[j];

        }
  

        Matrix3_9 d(3, 9);

        d << 1, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 1, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 1, 0, 0, 0, 0, 0, 0;

        eigenIntegrand = (matrix3dToEigenMatrix(fSecondMinusfFirst)*d*inverseShapeTensor*polynomialEigen)*influence;

        for (int j = 0; j < 3; j++)
        {
          bondIntegrand[j] = eigenIntegrand(j);

        }

        
        integrand = integrand + bondIntegrand*bond->getVolume();
        
      } // end of else if (getOrderFlag() == 2)
    
    } //end of if (bond->getBrokenBond() == false)


/*    else
    {
      bond->setForceVectorState(zero);
      bond->setBrokenBond(true);
      bond->setInfluenceFunction(0.0);
      StateBondP twin_bond = twinBond(bond);
      twin_bond->setBrokenBond(true);
      twin_bond->setForceVectorState(zero);
      twin_bond->setInfluenceFunction(0.0);

    }*/

  } // end of for (i = 0; i < size; i++)

  return (integrand);
}



Vector3D
StateMaterialPoint::nonLocalDivergenceMatrix3D(tensor f) // state based governing equation
{
  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  int i = 0;
  double influence = 0.0;
  Matrix3_1 eigenIntegrand = Matrix3_1::Zero();
  Vector3D integrand(0.0, 0.0, 0.0);
  Vector3D sai(0.0, 0.0, 0.0);
  Vector3D zero(0.0, 0.0, 0.0);
  Matrix3  eigenF = Matrix3::Zero();
  Matrix3D fSecond(0.0);
  Matrix3 eigenFSecond= Matrix3::Zero();

  Vector3D sum(0.0, 0.0, 0.0);


  Matrix3D invShapeTensor(0.0);
  Matrix3D twinInvShapeTensor(0.0);
  Matrix9  inverseShapeTensor = Matrix9::Zero();
  Matrix9  twinInverseShapeTensor = Matrix9::Zero();
  Matrix9_1 polynomialEigen = Matrix9_1::Zero();
  Matrix9_1 twinPolynomialEigen = Matrix9_1::Zero();

  std::array <double, 9>  polynomial{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  std::array <double, 9>  twinPolynomial{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  
  Vector3D bondIntegrand(0.0, 0.0,0.0);

  eigenF = matrix3dToEigenMatrix((this->*f)());

  if (getOrderFlag() == 1)
  {
    shapeTensor();
    invShapeTensor = invMatrix3D(getShapeTensor());  

  }

  if (getOrderFlag() == 2)
  {
    shapeTensorPolynomialOrder2();
    inverseShapeTensor = getShapeTensorPolynomialOrder2().inverse();  

  }

   Matrix3_9 d(3, 9);

   d << 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0;

  
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    if (bond->getBrokenBond() == false)
    {
      influence = bond->getInfluenceFunction();
      sai = bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
      StateMaterialPointP horizonPoint = bond->getSecondPoint();
      StateMaterialPoint point = *horizonPoint;
      fSecond = (point.*f)();
      eigenFSecond = matrix3dToEigenMatrix(fSecond);
//      Matrix3D fSecondMinusfFirst = fSecond - (this->*f)();


      if (getOrderFlag() == 1)
      {
        horizonPoint->shapeTensor();
        twinInvShapeTensor = invMatrix3D(horizonPoint->getShapeTensor());  

      }

      if (getOrderFlag() == 2)
      {
        horizonPoint->shapeTensorPolynomialOrder2();
        twinInverseShapeTensor = horizonPoint->getShapeTensorPolynomialOrder2().inverse();
      }

//      if (getShapeTensorPolynomialOrder2().determinant() == 0)
//      {
//        std::cout << "Hazard 1: Node " << getID() << " Shape Tensor has a zero determinant." << std::endl;
//      }


      if (getOrderFlag() == 1)
      {
        
        integrand = ((this->*f)()*invShapeTensor*bond->getSai())*influence;
        horizonPoint->shapeTensor();        
        integrand = integrand - (fSecond*twinInvShapeTensor*twinBond(bond)->getSai())*twinBond(bond)->getInfluenceFunction();
        sum = sum + integrand*bond->getVolume();
      }

      if (getOrderFlag() == 2)
      {
        bond->fillPolynomialBasisOrder2();
        polynomial = bond->getPolynomialBasisOrder2();

        for (int j = 0; j < 9; j++)
        {
          polynomialEigen(j) = polynomial[j];

        }

        eigenIntegrand = (eigenF*d*inverseShapeTensor*polynomialEigen)*influence;

        StateBondP twin_bond = twinBond(bond);
        twin_bond->fillPolynomialBasisOrder2();
        twinPolynomial = twin_bond->getPolynomialBasisOrder2();

        for (int j = 0; j < 9; j++)
        {
          twinPolynomialEigen(j) = twinPolynomial[j];

        }


        eigenIntegrand = eigenIntegrand - (eigenFSecond*d*twinInverseShapeTensor*twinPolynomialEigen)*twin_bond->getInfluenceFunction();

        for (int j = 0; j < 3; j++)
        {
          integrand[j] = eigenIntegrand(j);

        }

        sum = sum + integrand*bond->getVolume();      
      }

    }

/*    else
    {
      bond->setForceVectorState(zero);
      bond->setBrokenBond(true);
      bond->setInfluenceFunction(0.0);
      StateBondP twin_bond = twinBond(bond);
      twin_bond->setBrokenBond(true);
      twin_bond->setForceVectorState(zero);
      twin_bond->setInfluenceFunction(0.0);

    }*/
                     
 }

 return sum;  

}



void
StateMaterialPoint::internalForce()
{
////  if (getID() == 1055)
////    printScreenMatrix3D(getFirstPiolaKirchhoffStress(), "FirstPiloa");
  tensor firstPiolaKirchhoffStress = &StateMaterialPoint::getFirstPiolaKirchhoffStress;
  Vector3D intForce(0.0, 0.0, 0.0);

//  if ((getShapeTensorPolynomialOrder2().determinant() != 0.0) || (getShapeTensor().Determinant() != 0.0))
//  {  
    intForce = divergenceMatrix3D(firstPiolaKirchhoffStress);
//  }
//  else
//  {
//    setBondsInfluenceFunctionsZero();
//  }
    
  setInternalForce(intForce);

////  if (getID() == 1055)
////     std::cout << "InternalForce= " << getInternalForce() << std::endl;

}



void
StateMaterialPoint::internalForceUsingStateBasedPeridynamics()
{
  tensor firstPiolaKirchhoffStress = &StateMaterialPoint::getFirstPiolaKirchhoffStress;
  Vector3D intForce(0.0, 0.0, 0.0);

//  if ((getShapeTensorPolynomialOrder2().determinant() != 0.0) || (getShapeTensor().Determinant() != 0.0))
//  {  
    intForce = nonLocalDivergenceMatrix3D(firstPiolaKirchhoffStress);
//  }
//  else
//  {
//    setBondsInfluenceFunctionsZero();
//  }

  setInternalForce(intForce);

}



void
StateMaterialPoint::setBondsInfluenceFunctionsZero()
{
  double zero = 0.0;
  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  for (int i = 0; i < size; i++)
  {
    StateBondP cur_bond = bonds[i];
    cur_bond->setInfluenceFunction(zero);
  }



}



// END OF (4) **********************************************************************************************************************************************************************



// *********************************************************************************************************************************************************************************
// (5) ******************************************************** ObtainNextCurrentPositionSolvingConservationOfMomentum *************************************************************  



void
StateMaterialPoint::currentPositionVelocityAcceleration(const double& density, const double& delT)
{
  Vector3D zero(0.0, 0.0, 0.0);
  Vector3D intForce = getInternalForce();
/*  if (getID() == 1055)
  {
    std::cout << "InternalForce= (" << intForce.x() << ", " << intForce.y() << ", " << intForce.z() << ")" << std::endl;
    std::cout << "Acceleration= " << getAccele() << std::endl;
    std::cout << "Old Velocity= " << getVelOld() << std::endl;
    std::cout << "New Velocity= " << getVelNew() << std::endl;
    std::cout << "Displacement= " << getDisp() << std::endl << std::endl;
    std::cout << "Old Position= " << getPosOld() << std::endl;
    std::cout << "New Position= " << getPosNew() << std::endl;


  } */


  Vector3D acceleration = (intForce + getBodyForce())/density; // Conservation Of Momentum
  if (getVelocityBoundaryFlag() == true)
  {
    setAccele(zero);
  }
  else
  {
    setAccele(acceleration);
  }
  setVelNew(getVelNew() + (getAccele())*delT);
  setDisp(getDisp() + (getVelNew())*delT);
  setPosNew(getPosOld() + getDisp());

}


// END OF (5) **********************************************************************************************************************************************************************


// *********************************************************************************************************************************************************************************
// (Gathering the functions) *******************************************************************************************************************************************************



void
StateMaterialPoint::cauchyStress(const double& firstLamme, const double& secondLamme)
{
  deformationGradientUsingPolynomials();
  infinitesimalStrainTensor();
  isotropicLinearElasticityCauchyStressTensor(firstLamme, secondLamme);
  Matrix3D stress = getCauchyStress();
  firstPiolaKirchhoffStress(stress);

}


void
StateMaterialPoint::updatePosition(const double& density, const double& delT)
{
  internalForce();
  currentPositionVelocityAcceleration(density, delT);

}


void
StateMaterialPoint::updatePositionUsingStateBasedPeridynamics(const double& density, const double& delT)
{
  internalForceUsingStateBasedPeridynamics();
  currentPositionVelocityAcceleration(density, delT);

}


// **********************************************************************************************************************************************************************************


// ********************************************************* Implicit Time Integration ********************************************


void
StateMaterialPoint::updateBroydenMatrix(const double& alpha)
{
  setBroydenMatrix(getBroydenMatrix()*alpha);

}


Vector3D
StateMaterialPoint::internalForceLinearElasticImplicit(const double& firstLamme,
                                                       const double& secondLamme,
                                                       const int& governingType)
{
  int id = 1;
  Matrix3D deformationGradient(0.0);
  Matrix3D ident(0.0);
  ident.Identity();
  field dispNew = &StateMaterialPoint::returnDisplacementNew;
  deformationGradient = derivativeOperatorApproximationUsingSecondOrderPolynomialBasis(dispNew) + ident;  
  setDeformationGradient(deformationGradient);
//  if (getID() == id)
//  {
//    std::cout << "DeformationGradient= " << std::endl;
//    getDeformationGradient().printMatrix();
//  }
  infinitesimalStrainTensor();

  isotropicLinearElasticityCauchyStressTensor(firstLamme, secondLamme);

//  if (getID() == id)
//  {
//    std::cout << std::endl << "CauchyStress= " << std::endl;
//    getCauchyStress().printMatrix();
//  }


  Matrix3D stress = getCauchyStress();
  firstPiolaKirchhoffStress(stress);
  if (governingType == 0)
  {
    internalForce();  // Divergance Method
  }
  else
  {
    internalForceUsingStateBasedPeridynamics();
  }
  return (getInternalForce());
  

}



void
StateMaterialPoint::updateAccelerationDisplacementImplicitly(const double& firstLamme,
                                                             const double& secondLamme,
                                                             const double& alpha,
                                                             const double& beta,
                                                             const double& density,
                                                             const double& delT,
                                                             const double& tolerance,
                                                             const int& governingType)
{

  int id = 1050; 
  int max_iteration = 15000;
  int count = 0;

//  bool condition = false;
  bool condition1 = false;

  Matrix3D identity;
  identity.Identity();
  setBroydenMatrix(identity*alpha);

//  setBroydenMatrix(getBroydenMatrix()*alpha);

//  setDispNew(getDisp() - getBroydenMatrix()*getResidualVector());
//  double invDensity = 1/density;
  setDisp(getDispOld());
  do
  {
//    if (getID() == id)
//    {
//      std::cout << "id= " << getID() << " count= " << count << std::endl << std::endl;
//      std::cout << "Residual Vector before= " << getResidualVector() << std::endl << std::endl;
//    }

//  if (getID() == id)
//  {
//    std::cout << /*"id= " << getID() <<*/ "Broyden's Matrix before= " << std::endl;
//    getBroydenMatrix().printMatrix();
//    std::cout << std::endl;

//  }

    setDispNew(getDisp() - getBroydenMatrix()*getResidualVector());
//    if (getID() == id)
//    {
//      std::cout << "DispNew= " << getDispNew() << std::endl << std::endl;
//    }
    Vector3D intForce = internalForceLinearElasticImplicit(firstLamme, secondLamme, governingType);
//    if (getID() == id)
//    {
//      std::cout << "intForce= " << intForce << std::endl;
//      std::cout << "bodyForce= " << getBodyForce() << std::endl;
//      std::cout << "density= " << density << std::endl;
//    }
         
    setAcceleNew((intForce + getBodyForce())/density);    
    Vector3D eta = getDispNew() - getDispOld();
    Vector3D currentEta = getDispNew() - getDisp();
//    if (getID() == id)
//    {
//      std::cout << "dispNew - disp= " << currentEta << std::endl << std::endl;
//    } 
    Vector3D residual = eta - getVelOld()*delT - getAccele()*(delT*delT*(1-2*beta)*0.5) - getAcceleNew()*(delT*delT*beta); 
//    Vector3D residual = eta - getVelOld()*delT;
//    residual = residual - getAccele()*(delT*delT*(1-2*beta)*0.5);
//    residual = residual - getAcceleNew()*(delT*delT*beta);

    Vector3D du = getBroydenMatrix()*(residual - getResidualVector());
//    if (getID() == id)
//    {
//      std::cout << "F_new - F_current= " << residual - getResidualVector() << std::endl;
//      std::cout << "Acelec= " << getAccele() << std::endl;
//      std::cout << "AcelecNew= " << getAcceleNew() << std::endl;      
//      std::cout << "du= " << du << std::endl << std::endl;
//    }
    double dotDuAndCurrentEta = currentEta.dot(du);
//    if (getID() == id)
//    {    
//      std::cout << "denominator= " << dotDuAndCurrentEta << std::endl << std::endl;
//    }
    if (dotDuAndCurrentEta != 0)
    {
      Matrix3D tensorCurrentEtaMinusDuAndEta(Matrix3D(currentEta - du, currentEta));  
//      if (getID() == id)
//      {
//        std::cout << "tensor product of currentEtaMinusDu and du= " << std::endl;
//        tensorCurrentEtaMinusDuAndEta.printMatrix();
//        std::cout << std::endl << "Numerator= " << std::endl;
//        (tensorCurrentEtaMinusDuAndEta*getBroydenMatrix()).printMatrix();
//      }
      setBroydenMatrix(getBroydenMatrix() + (tensorCurrentEtaMinusDuAndEta*getBroydenMatrix())/dotDuAndCurrentEta);
    }
    setResidualVector(residual);

//    if (getID() == id)
//    {
//      std::cout << "id= " << getID() << " count= " << count << std::endl << std::endl;
//      std::cout << "Residual Vector after= " << getResidualVector() << std::endl << std::endl;
//    }


//  if (getID() == id)
//  {
//    std::cout << /*"id= " << getID() <<*/ "Broyden's Matrix after= " << std::endl;
//    getBroydenMatrix().printMatrix();
//    std::cout << std::endl;

//  }

//    condition = (getResidualVector().length() > tolerance);
    condition1 = ((getDispNew() - getDisp()).length() > tolerance);
    condition1 = condition1 && (count < max_iteration);

//    if (getID() == id)
//    {
//      std::cout << std::endl << "residualVector= " << getResidualVector() << std::endl;

//      std::cout << "Tolerance= " << tolerance << std::endl;
//    }
    count++;
//    std::cout << std::endl << "residualVector= " << getResidualVector() << std::endl;
//    std::cout << std::endl << "residualVectorLength= " << getResidualVector().length() << std::endl;
    setDisp(getDispNew());

  } while (condition1);

  setDispNew(getDisp());

  if ((getShapeTensorPolynomialOrder2().determinant() != 0.0) || (getShapeTensor().Determinant() != 0.0))
  {

    setBondsInfluenceFunctionsZero();

  }


}










void
StateMaterialPoint::updateAccelerationDisplacementImplicitly(const double& firstLamme,
                                                             const double& secondLamme,
                                                             const double& alpha,
                                                             const double& beta,
                                                             const double& density,
                                                             const double& delT,
                                                             const double& tolerance,
                                                             const int& governingType,
                                                             const int& num_iter,
                                                             const std::string& currentFolder)
{

//  int id = 1050; 
  int id = 3909;
  int id2 = 4184;
 
  int max_iteration = 15000;
  int count = 0;

//  bool condition = false;
  bool condition1 = false;


  Matrix3D identity;
  identity.Identity();
  setBroydenMatrix(identity*alpha);


//  bool num_condition = (num_iter > 10 && num_iter < 45);
//  bool num_condition = (num_iter > 10 && num_iter < 64);
  bool num_condition = (num_iter < 490);


//  setBroydenMatrix(getBroydenMatrix()*alpha);

//  setDispNew(getDisp() - getBroydenMatrix()*getResidualVector());
//  double invDensity = 1/density;

//  std::string fileName = "details_" + std::to_string(num_iter) + ".com";
  std::string fileName = "details_" + std::to_string(getID()) + "_" + std::to_string(num_iter) + ".com";


  int ret;


//  bool writing_condition = false;

  bool writing_condition = ((getID() == id) && num_condition);
  writing_condition = (writing_condition) || (getID() == id2);

  if (writing_condition)
  {
//    std::string fileName = "details_" + std::to_string(num_iter) + ".com";
    std::ofstream myfile(fileName);
    myfile << "num_iter= " << num_iter << std::endl;  
    myfile << "Node id= " << getID() << std::endl;
    myfile << "Referance Position= " << getPosOld() << std::endl;
    myfile << "Velocity Old[" << num_iter << "]= " << getVelOld() << std::endl;
    myfile << "Accele[" << num_iter << "]= " << getAccele() << std::endl;
    myfile << "DispOld[" << num_iter << "]= " << getDispOld() << std::endl;
    myfile << "NewPosition[" << num_iter << "]= " << getPosNew() << std::endl;
    myfile.close();
    ret = chdir(currentFolder.c_str());
    
  }


  setDisp(getDispOld());
  do
  {
    if (writing_condition)
    {
//      std::cout << "id= " << getID() << " count= " << count << std::endl << std::endl;
//      std::cout << "Residual Vector before= " << getResidualVector() << std::endl << std::endl;
      std::ofstream myfile;
      myfile.open(fileName, std::ofstream::app);
      myfile << std::endl << "id= " << getID() << " count= " << count << std::endl;
      myfile << "Residual Vector: " << std::endl << "F[" << count << "]= " << getResidualVector() << std::endl;
      myfile << "Broyden's Matrix: " << std::endl << "B[" << count << "]= " << std::endl << getBroydenMatrix();
      myfile << "Disp[" << count << "]= " << getDisp() << std::endl;
//      myfile << std::endl;
      myfile.close();
      ret = chdir(currentFolder.c_str());
    }

//    if (getID() == id)
//    {
//    std::cout << /*"id= " << getID() <<*/ "Broyden's Matrix before= " << std::endl;
//    getBroydenMatrix().printMatrix();
//    std::cout << std::endl;

//    }

    setDispNew(getDisp() - getBroydenMatrix()*getResidualVector());


//    if (getID() == id)
//    {
//      std::cout << "DispNew= " << getDispNew() << std::endl << std::endl;
//    }

    if (writing_condition)
    {
      std::ofstream myfile;
      myfile.open(fileName, std::ofstream::app);
      myfile << "DispNew[" << count + 1 << "]= " << getDispNew() << std::endl;
      myfile.close();
      ret = chdir(currentFolder.c_str());
    }

    Vector3D intForce = internalForceLinearElasticImplicit(firstLamme, secondLamme, governingType);

//    if (getID() == id)
//    {
//      std::cout << "intForce= " << intForce << std::endl;
//      std::cout << "bodyForce= " << getBodyForce() << std::endl;
//      std::cout << "density= " << density << std::endl;
//    }

    if (writing_condition)
    {
      std::ofstream myfile;
      myfile.open(fileName, std::ofstream::app);
      myfile << "density= " << density << std::endl;
      myfile << "intForce[" << count + 1 << "]= " << intForce << std::endl;
      myfile << "bodyForce[" << num_iter << "]= " << getBodyForce() << std::endl;
      myfile.close();
      ret = chdir(currentFolder.c_str());

    }
         
    setAcceleNew((intForce + getBodyForce())/density);    
    Vector3D eta = getDispNew() - getDispOld();
    Vector3D currentEta = getDispNew() - getDisp();

//    if (getID() == id)
//    {
//      std::cout << "dispNew - disp= " << currentEta << std::endl << std::endl;
//    } 

    Vector3D residual = eta - getVelOld()*delT - getAccele()*(delT*delT*(1-2*beta)*0.5) - getAcceleNew()*(delT*delT*beta); 
    Vector3D du = getBroydenMatrix()*(residual - getResidualVector());

//    if (getID() == id)
//    {
//      std::cout << "F_new - F_current= " << residual - getResidualVector() << std::endl;
//      std::cout << "Acelec= " << getAccele() << std::endl;
//      std::cout << "AcelecNew= " << getAcceleNew() << std::endl;      
//      std::cout << "du= " << du << std::endl << std::endl;
//    }

    if (writing_condition)
    {
      std::ofstream myfile;
      myfile.open(fileName, std::ofstream::app);
      myfile << "Acelec[" << count << "]= " << getAccele() << std::endl;
      myfile << "AcelecNew[" << count + 1 << "]= " << getAcceleNew() << std::endl;
      myfile << "eta = DispNew[" << count + 1 << "]- DispOld[" << num_iter << "]= " << getDispNew() - getDispOld() << std::endl;     
      myfile << "VelOld*delT= " << getVelOld()*delT << ", Accele*(delT*delT*(1-2*beta)*0.5)= " << getAccele()*(delT*delT*(1-2*beta)*0.5)
             << ", AcceleNew*(delT*delT*beta)= " << getAcceleNew()*(delT*delT*beta) << std::endl;
      myfile << "Current Residual Vector= " << "F[" << count + 1 << "]= " << residual << std::endl;
      myfile << "F[" << count + 1 << "]- F[" << count << "]= " << residual - getResidualVector() << std::endl;
      myfile << "B[" << count << "]*{F[" << count + 1 << "]- F[" << count << "]}= du= " << getBroydenMatrix()*(residual - getResidualVector()) << std::endl;
      myfile << "currentEta = DispNew[" << count + 1 << "]- Disp[" << count << "]= " << getDispNew() - getDisp() << std::endl; 

      myfile.close();
      ret = chdir(currentFolder.c_str());

     }



    double dotDuAndCurrentEta = currentEta.dot(du);

//    if (getID() == id)
//    {    
//      std::cout << "denominator= " << dotDuAndCurrentEta << std::endl << std::endl;
//    }

    if (dotDuAndCurrentEta != 0)
    {
      Matrix3D tensorCurrentEtaMinusDuAndEta(Matrix3D(currentEta - du, currentEta));
  
//      if (getID() == id)
//      {
//        std::cout << "tensor product of currentEtaMinusDu and du= " << std::endl;
//        tensorCurrentEtaMinusDuAndEta.printMatrix();
//        std::cout << std::endl << "Numerator= " << std::endl;
//        (tensorCurrentEtaMinusDuAndEta*getBroydenMatrix()).printMatrix();
//      }

      setBroydenMatrix(getBroydenMatrix() + (tensorCurrentEtaMinusDuAndEta*getBroydenMatrix())/dotDuAndCurrentEta);

      if (writing_condition)
      {
        std::ofstream myfile;
        myfile.open(fileName, std::ofstream::app);
        myfile << "{DispNew[" << count + 1 << "]- Disp[" << count << "] dotProduct B[" 
               << count << "]*{F[" << count + 1 << "]- F[" << count << "]}}= dotDuAndCurrentEta= " << dotDuAndCurrentEta << std::endl;
        myfile << "{DispNew[" << count + 1 << "]- Disp[" << count << "] tensorProduct B[" 
               << count << "]*{F[" << count + 1 << "]- F[" << count << "]}}= tensorCurrentEtaMinusDuAndEta= " << std::endl << tensorCurrentEtaMinusDuAndEta;
        myfile << "tensorCurrentEtaMinusDuAndEta*getBroydenMatrix())/dotDuAndCurrentEta= " 
               << std::endl << (tensorCurrentEtaMinusDuAndEta*getBroydenMatrix())/dotDuAndCurrentEta;
        myfile << "Current Broyden's Matrix= " << std::endl << "B[" << count + 1 << "]= " << std::endl << getBroydenMatrix();         
        myfile.close();
        ret = chdir(currentFolder.c_str());

      }


    }
    setResidualVector(residual);

//    if (getID() == id)
//    {
//      std::cout << "id= " << getID() << " count= " << count << std::endl << std::endl;
//      std::cout << "Residual Vector after= " << getResidualVector() << std::endl << std::endl;
//    }


//  if (getID() == id)
//  {
//    std::cout << /*"id= " << getID() <<*/ "Broyden's Matrix after= " << std::endl;
//    getBroydenMatrix().printMatrix();
//    std::cout << std::endl;

//  }

//    condition = (getResidualVector().length() > tolerance);
    condition1 = ((getDispNew() - getDisp()).length() > tolerance);
    condition1 = condition1 && (count < max_iteration);

//    if (getID() == id)
//    {
//      std::cout << std::endl << "residualVector= " << getResidualVector() << std::endl;

//      std::cout << "Tolerance= " << tolerance << std::endl;
//    }
    count++;
//    std::cout << std::endl << "residualVector= " << getResidualVector() << std::endl;
//    std::cout << std::endl << "residualVectorLength= " << getResidualVector().length() << std::endl;
    setDisp(getDispNew());

  } while (condition1);

  setDispNew(getDisp());

  if (writing_condition)
  {
    std::ofstream myfile;
    myfile.open(fileName, std::ofstream::app);

    myfile << std::endl << "FinalDispNew[" << num_iter + 1 << "]= " << getDispNew() << std::endl;
    myfile.close();
    ret = chdir(currentFolder.c_str());

  }

  if ((getShapeTensorPolynomialOrder2().determinant() == 0.0) /*|| (getShapeTensor().Determinant() == 0.0)*/)
  {

    setBondsInfluenceFunctionsZero();

  }

//  myfile.close();
//  ret = chdir(currentFolder.c_str());
 

}




void
StateMaterialPoint::getInformation(Matrix3D& BroydenCurrent,
                                   Matrix3D& BroydenNew,
                                   Vector3D& residualVectorCurrent,
                                   Vector3D& residualVectorNew)
{
  










}

// ********************************************************************************************************************************





StateBondP
StateMaterialPoint::twinBond(const StateBondP& bond)
{
  StateMaterialPointP second_point = bond->getSecondPoint();
  StateMaterialPointP center_point = bond->getCenterPoint();
  StateBondPArray bondsArray = second_point->getStateBonds();    
  int boundSize = bondsArray.size();
  int j = 0;
  for (j = 0; j < boundSize; j++)
  { 
    StateBondP twinBond = bondsArray[j];
    if (twinBond->getSecondPoint() == center_point)
    {
      return twinBond;
    }

  }    

}






Matrix3D
StateMaterialPoint::materialTimeDerivative()
{
  Matrix3D timeDerivative(0.0);
  Matrix3D tensor(0.0);
  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  int i = 0;  
//  Vector3D vecVelBond(0.0);
//  Vector3D vecBond(0.0);
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    if (bond->getBrokenBond() == false)
    {
//      vecBond = bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      vecVelBond = bond->getSecondPoint()->getVelNew() - bond->getCenterPoint()->getVelNew();
//      Matrix3D tensor(vecVelBond, vecBond);
      tensor = bond->velocityVecTensorSai();
      timeDerivative = timeDerivative + tensor*(bond->getVolume()*bond->getInfluenceFunction());
//      if ((bond->getCenterPoint()->getID() == 700)) // && (bond->getSecondPoint()->getID() == 1626))
//      {
//        std::cout << bond->getSecondPoint()->getID() << std::endl;
//        std::cout << bond->getSai() << std::endl;
//        std::cout << bond->getVolume() << std::endl;
////        std::cout << vecVelBond << std::endl;
//        std::cout << bond->getVelocityVector() << std::endl;
//        printScreenMatrix3D(tensor, "velocityTensorSai");
//        printScreenMatrix3D(timeDerivative, "timeDerivative");
//      } 
    }


  }

//  timeDerivative = timeDerivative*getShapeTensor().Inverse();
  if (getShapeTensor().Determinant() == 0)
  {
    std::cout << "Hazard 2: Node " << getID() << " Shape Tensor has a zero determinant." << std::endl;
  }
  timeDerivative = timeDerivative*invMatrix3D(getShapeTensor());
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(timeDerivative, "finalTimeDerivative");
//  }

  return timeDerivative;
  
}


Matrix3D
StateMaterialPoint::materialTimeDerivativeUsingPolynomials()
{

  Matrix3D velocity(0.0);

  field currentVelocity = &StateMaterialPoint::returnCurrentVelocity;

  if (getOrderFlag() == 1)
  {

    velocity = derivativeOperatorApproximationUsingFirstOrderPolynomialBasis(currentVelocity);

  }
  
  if (getOrderFlag() == 2)
  {

    velocity = derivativeOperatorApproximationUsingSecondOrderPolynomialBasis(currentVelocity);

  }

  return velocity;

}




void
StateMaterialPoint::spatialVelocityGradient()
{

//  setVelocityGradient(materialTimeDerivative()*getDeformationGradient().Inverse());
//  setVelocityGradient(materialTimeDerivative()*invMatrix3D(getDeformationGradient()));
  setVelocityGradient(materialTimeDerivativeUsingPolynomials()*invMatrix3D(getDeformationGradient()));

}


void
StateMaterialPoint::spatialVelocityGradient(const int& numIter,
                                            const Matrix3D& initialVelocityGradient)
{

//  setVelocityGradient(materialTimeDerivative()*getDeformationGradient().Inverse());
//  int num_initial = 10;
//  int num_initial = 21;
  int num_initial = 3;
//  if (numIter == 0)
//  if (numIter < 10)
  if (numIter < num_initial)
  {
//    setVelocityGradient(initialVelocityGradient*ramp(numIter, num_initial));
    setVelocityGradient(initialVelocityGradient*triangle(numIter, num_initial));    
  }
  else
  {
    if (getDeformationGradient().Determinant() == 0)
    {
      std::cout << "Hazard 3: Node " << getID() << " Deformation gradient has a zero determinant." << std::endl;
    }
//    setVelocityGradient(materialTimeDerivative()*invMatrix3D(getDeformationGradient()));
    setVelocityGradient(materialTimeDerivativeUsingPolynomials()*invMatrix3D(getDeformationGradient()));
  }
//  if (getID() == 4224)
//  {
//    printScreenMatrix3D(getVelocityGradient());
//    std::cout << "numIter= " << numIter << " num_initial= " << num_initial << std::endl;
//    std::cout << "ramp= " << ramp(numIter, num_initial) << std::endl;
//  }

//  if (getID() == 4224)
//  {
//    printScreenMatrix3D(getVelocityGradient());
//    std::cout << "numIter= " << numIter << " num_initial= " << num_initial << std::endl;
//    std::cout << "triangle= " << triangle(numIter, num_initial) << std::endl;
//  }
}



void
StateMaterialPoint::deformationRateTensor()
{
  setDeformationRateTensor((getVelocityGradient()+getVelocityGradient().Transpose())*0.5);

}


void
StateMaterialPoint::spinTensor()
{
  setSpinTensor((getVelocityGradient()-getVelocityGradient().Transpose())*0.5);

}


int
StateMaterialPoint::permutationTensor(const int& i, const int& j, const int& k)
{
  int p = 0;
  bool right = (((i == 0) && (j == 1) && (k == 2)) ||
                ((i == 1) && (j == 2) && (k == 0)) ||
                ((i == 2) && (j == 0) && (k == 1)));

  bool left = (((i == 0) && (j == 2) && (k == 1)) ||
               ((i == 2) && (j == 1) && (k == 0)) ||
               ((i == 1) && (j == 0) && (k == 2)));

  if (right)
  {
    p = 1;
  }

  if (left)
  {
    p = -1;
  }
  
  return p;
}


Vector3D
StateMaterialPoint::dualVector(Matrix3D skew)  // 2w(i) = -e(i,j,k)W(j,k), where W is an antisymmetric tensor rank 2
{
  Vector3D vec(0.0);
  int i = 0;
  int j = 0;
  int k = 0;
  double sum = 0.0;
  for (i = 0; i < 3; i++)
  {
    sum = 0.0;
    for (j = 0; j < 3; j++)
    {
       for (k = 0; k < 3; k++)
       {
         sum = sum + (-1)*permutationTensor(i, j, k)*skew(j, k);  
       }
    }
//    vec[i] = sum;
    vec[i] = sum*(0.5);
  }

  return vec;

}


Matrix3D
StateMaterialPoint::dualSkewMatrix(Vector3D vector)  // W(i,j) = e(i,k,j)w(k)
{
  Matrix3D matrix(0.0);
  int i = 0;
  int j = 0;
  int k = 0;
  double sum = 0.0;  
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      sum = 0.0;
      for (k = 0; k < 3; k++)
      {
        sum = sum + permutationTensor(i, k, j)*vector[k];  
      }
      matrix(i, j) = sum;
    }

  }

  return matrix;

}



MatrixXd
StateMaterialPoint::matrix3dToEigenMatrix(const Matrix3D& matrix)
{
  int i = 0;
  int j = 0;
  MatrixXd m(3,3);
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      m(i,j) = matrix(i,j);
    }
  }
  return m;
}
  

Matrix3D
StateMaterialPoint::eigenMatrixToMatrix3D(const MatrixXd& matrix)
{
  int i = 0;
  int j = 0;
  Matrix3D m;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      m(i,j) = matrix(i,j);
    }
  }
  return m;
}



void
StateMaterialPoint::polarDecomposition(const Matrix3D& deformation,
                                             Matrix3D& rotation,
                                             Matrix3D& right,
                                             Matrix3D& left)
{
  MatrixXd def(3,3);
  def = matrix3dToEigenMatrix(deformation);
  MatrixXd rightCauchyGreenDefomation = def.transpose()*def;
  SelfAdjointEigenSolver<MatrixXd> es(rightCauchyGreenDefomation);  
  MatrixXd U(3,3);
  U = es.operatorSqrt();
  MatrixXd R(3,3);
  R = def*U.inverse();
  MatrixXd V(3,3);    
  V = R*U*R.transpose();
  rotation = eigenMatrixToMatrix3D(R);
  right = eigenMatrixToMatrix3D(U);
  left = eigenMatrixToMatrix3D(V);

}


void
StateMaterialPoint::polarDecompositionOfDeformationGradient()
{
  Matrix3D deformation = getDeformationGradient();  
//  Matrix3D rotation = getRotationTensor();
  Matrix3D rotation(0.0);
//  Matrix3D right = getRightStretch();
  Matrix3D right(0.0);
//  Matrix3D left = getLeftStretch();
  Matrix3D left(0.0);
  polarDecomposition(deformation, rotation, right, left);
/*  if (getID() == 700)
  {
    printScreenMatrix3D(getShapeTensor(), "shapeTensor");
    printScreenMatrix3D(invMatrix3D(getShapeTensor()), "shapeTensorInverse");
    printScreenMatrix3D(getShapeTensor()*invMatrix3D(getShapeTensor()), "identity");
    printScreenMatrix3D(invMatrix3D(getShapeTensor())*getShapeTensor(), "identity");
    printScreenMatrix3D(deformation, "Deformation");
    printScreenMatrix3D(rotation, "RotationTensor");
    printScreenMatrix3D(rotation*rotation.Transpose(), "identity");
    printScreenMatrix3D(rotation.Transpose()*rotation, "identity");
    printScreenMatrix3D(right, "Right Stretch Tensor"); 
    printScreenMatrix3D(left, "Leftt Stretch Tensor");
  }  */ 
  setRotationTensor(rotation);
  setRightStretch(right);
  setLeftStretch(left);

}


void
StateMaterialPoint::angularVelocity()
{
  Vector3D w = dualVector(getSpinTensor());
  Matrix3D V = getLeftStretch();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getSpinTensor(), "spinTensor");
//    std::cout << "dual vector of the spin tensor= " << w << std::endl << std::endl;
//    printScreenMatrix3D(getLeftStretch(), "leftStretch V"); 
//  }       
  double traceV = V(0,0)+V(1,1)+V(2,2);
  Matrix3D coeffI(traceV, 0.0, 0.0,
                  0.0, traceV, 0.0,
                  0.0, 0.0, traceV);
//  Matrix3D matrix = (coeffI - V).Inverse();
  Matrix3D matrix = invMatrix3D((coeffI - V));
  Vector3D z(0.0, 0.0, 0.0);
  int i = 0;
  int j = 0;
  int k = 0;
  int m = 0;
  double sum = 0.0;
  for (i = 0; i < 3; i++)
  { 
    sum = 0.0;
    for (j = 0; j < 3; j++)
    {
      for (k = 0; k < 3; k++)
      {
        for (m = 0; m < 3; m++)
        {
          sum = sum + permutationTensor(i,k,j)*getDeformationRateTensor()(j,m)*V(m,k);
        }
      }
    }
    z[i] = sum;
  }

  Vector3D omega = w + matrix*z;
  setAngularVelocity(omega);

//  if (getID() == 700)
//  {
//    std::cout << "dual vector of the omega= " << getAngularVelocity() << std::endl << std::endl;
//  }         

}



void
StateMaterialPoint::rateRigidBodyRotation()
{
  angularVelocity();
  setRateRigidBodyRotation(dualSkewMatrix(getAngularVelocity()));

//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getRateRigidBodyRotation(), "omega"); 
//  }   
  
}


void
StateMaterialPoint::updateRotationTensorMoreEfficient(const double& delT)
{
  Matrix3D I;
  I.Identity();
  rateRigidBodyRotation(); 
  Matrix3D omega = getRateRigidBodyRotation();
  omega = omega*(0.5*delT);
//  Matrix3D Q = (I - omega).Inverse()*(I + omega);
  Matrix3D Q = invMatrix3D((I - omega))*(I + omega);
  setRotationTensor(Q*getRotationTensor());

}


void
StateMaterialPoint::updateRotationTensorMoreStable(const double& delT)
{
  Matrix3D I;
  I.Identity();
  rateRigidBodyRotation();   
  Matrix3D omega = getRateRigidBodyRotation();
  Vector3D vectorOmega = getAngularVelocity();
//  double squareLengthOmega = vectorOmega.lengthSq();
  double lengthOmega = vectorOmega.length();
  double coeff1 = delT;
  double coeff2 = 0.5*delT*delT;
//  if (!(squareLengthOmega == 0))
  if (!(lengthOmega == 0))
  {
//    coeff1 = sin(delT*squareLengthOmega)/(squareLengthOmega);
//    coeff2 = (1 - cos(delT*squareLengthOmega))/(squareLengthOmega*squareLengthOmega);
    coeff1 = sin(delT*lengthOmega)/(lengthOmega);
    coeff2 = (1 - cos(delT*lengthOmega))/(lengthOmega*lengthOmega);

  }
  Matrix3D Q = I + omega*coeff1 - (omega*omega)*coeff2;
  setRotationTensor(Q*getRotationTensor());

}


void
StateMaterialPoint::leftStretchRateTensor()
{
  Matrix3D Vdot = getVelocityGradient()*getLeftStretch() -
                  getLeftStretch()*getRateRigidBodyRotation();

  setLeftStretchRate(Vdot);

}


void
StateMaterialPoint::updateLeftStretch(const double& delT)
{
  setLeftStretch(getLeftStretch() + getLeftStretchRate()*delT);

}


void
StateMaterialPoint::unrotatedDeformationRate()
{
  Matrix3D matrix = getDeformationRateTensor()*getRotationTensor();
//  setUnrotatedDeformationGradient(getRotationTensor().Transpose()*
//                                  getDeformationRateTensor()*getRotationTensor());
  setUnrotatedDeformationGradient((getRotationTensor().Transpose())*matrix);


}


void
StateMaterialPoint::strainIncrement(const double& delT)
{
  setTotalElasticStrianIncrement(getUnrotatedDeformationGradient()*delT);
  Matrix3D I;
  I.Identity();
  double trace = getUnrotatedDeformationGradient().Trace();    
  I = I*trace;
  I = I*(1.0/3.0);
  Matrix3D Deviatoric = getUnrotatedDeformationGradient() - I;
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getUnrotatedDeformationGradient(), "UnrotatedDeformationGradient");
//    printScreenMatrix3D(I);    
//    printScreenMatrix3D(Deviatoric, "Deviatoric");
//  }   
  setDeviatoricStrianIncrement(Deviatoric*delT);
//  setDeviatoricStrianIncrement(getTotalElasticStrianIncrement() - 
//                               I*((1/3)*getTotalElasticStrianIncrement().Trace()));

}



void
StateMaterialPoint::updateUnrotatedCauchyStress(const double& firstLamme,
                                                const double& secondLamme)
{
  Matrix3D I;
  I.Identity();  
  setUnrotatedCauchyStress(getUnrotatedCauchyStress() +
                           I*(firstLamme*getTotalElasticStrianIncrement().Trace()) +
//                           getDeviatoricStrianIncrement()*(2*secondLamme));
                           getTotalElasticStrianIncrement()*(2*secondLamme));
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(I*(firstLamme*getTotalElasticStrianIncrement().Trace()));
//    printScreenMatrix3D(getTotalElasticStrianIncrement(), "TotalElasticStrianIncrement");
//    printScreenMatrix3D(getDeviatoricStrianIncrement(), "DeviatoricStrianIncrement");
//    printScreenMatrix3D(getDeviatoricStrianIncrement()*(2*secondLamme));    
//  }    
  

}


void
StateMaterialPoint::rotatedCauchyStress()
{
  Matrix3D matrix = getUnrotatedCauchyStress()*(getRotationTensor().Transpose());
//  setRotatedCauchyStress(getRotationTensor()*getUnrotatedCauchyStress()*
//                         getRotationTensor().Transpose());
  setRotatedCauchyStress(getRotationTensor()*matrix);


}


void
StateMaterialPoint::firstPiolaKirchhoffStress()
{
  double J = getDeformationGradient().Determinant();
  Matrix3D Ftranspose = getDeformationGradient().Transpose();
//  if (getID() == 700)
//  {
//    std::cout << "Jacobian= " << J << std::endl;
//    printScreenMatrix3D(getDeformationGradient(), "deformationGradient");
//    printScreenMatrix3D(Ftranspose, "deformationGradientTranspose");
//  }
//  std::cout << Ftranspose(0,0) << " " << Ftranspose(0,1) << " " << Ftranspose(0,2) << std::endl;
//  std::cout << Ftranspose(1,0) << " " << Ftranspose(1,1) << " " << Ftranspose(1,2) << std::endl;
//  std::cout << Ftranspose(2,0) << " " << Ftranspose(2,1) << " " << Ftranspose(2,2) << std::endl;
//  Matrix3D FtransposeInverse = Ftranspose.Inverse();
  if (Ftranspose.Determinant() == 0)
  {
    std::cout << "Hazard 4: Node " << getID() << " Deformation gradient has a zero determinant." << std::endl;
  }
  Matrix3D FtransposeInverse = invMatrix3D(Ftranspose);
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(FtransposeInverse, "deformationGradientTransposeInverse");
//    printScreenMatrix3D(FtransposeInverse*Ftranspose, "identity");
//  }
  setFirstPiolaKirchhoffStress((getRotatedCauchyStress()*J)*(FtransposeInverse));

}


/*void
StateMaterialPoint::forceVectorStateOfBonds()
{
  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  int i = 0;
  double influence = 0.0;
  Vector3D sai(0.0, 0.0, 0.0);
  Vector3D zero(0.0, 0.0, 0.0);
  Matrix3D matrix(0.0);  
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    if (bond->getBrokenBond() == false)
    {
      influence = bond->getInfluenceFunction();
      sai = bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      matrix = getFirstPiolaKirchhoffStress()*getShapeTensor().Inverse();
      if (getShapeTensor().Determinant() == 0)
      {
        std::cout << "Hazard: Node " << getID() << " Shape Tensor has a zero determinant." << std::endl;
      }

      matrix = getFirstPiolaKirchhoffStress()*invMatrix3D(getShapeTensor());  
      bond->setForceVectorState((matrix*sai)*influence);
    }
    else
    {
      bond->setForceVectorState(zero);
      bond->setBrokenBond(true);
      bond->setInfluenceFunction(0.0);
    }

//    if (getID() == 700)
//    {
//      std::cout << bond->getCenterPoint()->getID() << ", "
//                << bond->getSecondPoint()->getID() << std::endl;
//      std::cout << "sai= " << bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      std::cout << " BondForceState= " << bond->getForceVectorState() << std::endl;    
//    }

//    bool condition1 = ((bond->getCenterPoint()->getID() == 700) &&      
//                      (bond->getSecondPoint()->getID() == 613));

//    bool condition2 = ((bond->getCenterPoint()->getID() == 613) &&      
//                      (bond->getSecondPoint()->getID() == 700));
//    if (condition1)
//    {
//      printScreenMatrix3D(getShapeTensor(), "ShapeTensor");
//      printScreenMatrix3D(invMatrix3D(getShapeTensor()), "ShapeTensorInverse");
//      printScreenMatrix3D(getFirstPiolaKirchhoffStress()*invMatrix3D(getShapeTensor()), "PKinverse");
//      std::cout << bond->getCenterPoint()->getID() << ", "
//                << bond->getSecondPoint()->getID() << std::endl;
//      std::cout << "sai= " << bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      std::cout << " BondForceState= " << bond->getForceVectorState() << std::endl;    
//    }
    
//    if (condition2)
//    {
//      printScreenMatrix3D(getShapeTensor(), "ShapeTensor");
//      printScreenMatrix3D(invMatrix3D(getShapeTensor()), "ShapeTensorInverse");
//      printScreenMatrix3D(getFirstPiolaKirchhoffStress()*invMatrix3D(getShapeTensor()), "PKinverse");
//      std::cout << bond->getCenterPoint()->getID() << ", "
//                << bond->getSecondPoint()->getID() << std::endl;
//      std::cout << "sai= " << bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      std::cout << " BondForceState= " << bond->getForceVectorState() << std::endl;    
//    }

  }


}*/


void
StateMaterialPoint::forceVectorStateOfBonds()
{
  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  int i = 0;
  double influence = 0.0;
  Vector3D sai(0.0, 0.0, 0.0);
  Vector3D zero(0.0, 0.0, 0.0);
  Matrix3D matrix(0.0);  
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    if (bond->getBrokenBond() == false)
    {
      influence = bond->getInfluenceFunction();
      sai = bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      matrix = getFirstPiolaKirchhoffStress()*getShapeTensor().Inverse();
//      if (getShapeTensor().Determinant() == 0)
//      {
//        std::cout << "Hazard 5: Node " << getID() << " Shape Tensor has a zero determinant." << std::endl;
//      }

      if (getOrderFlag() == 1)
      {
        matrix = getFirstPiolaKirchhoffStress()*invMatrix3D(getShapeTensor());  
        bond->setForceVectorState((matrix*sai)*influence);

////        Matrix3D piolaSecondMinusPiolaFirst = bond->getSecondPoint()->getFirstPiolaKirchhoffStress()-getFirstPiolaKirchhoffStress();
////        matrix = piolaSecondMinusPiolaFirst*invMatrix3D(getShapeTensor());
////        bond->setFirstPiolaDivergenceIntegrand((matrix*sai)*influence);
      }
      else if (getOrderFlag() == 2)
      {
        Matrix9_3 d(9,3);
        d << 1, 0, 0,
             0, 1, 0,
             0, 0, 1,
             0, 0, 0,
             0, 0, 0,
             0, 0, 0,
             0, 0, 0,
             0, 0, 0,
             0, 0, 0;

         Matrix3 piola = matrix3dToEigenMatrix(getFirstPiolaKirchhoffStress());

////         Matrix3 piolaSecondMinusPiolaFirst = matrix3dToEigenMatrix(bond->getSecondPoint()->getFirstPiolaKirchhoffStress()-getFirstPiolaKirchhoffStress());  // P(X+delX) - P(X)
         
         Matrix3_9 W = piola*d.transpose()*getShapeTensorPolynomialOrder2().inverse();

////         Matrix3_9 WPrime = piolaSecondMinusPiolaFirst*d.transpose()*getShapeTensorPolynomialOrder2().inverse(); // [P(X+delX) - P(X)] * d^T * K_inverse

         std::array <double, 9>  polynomial = bond->getPolynomialBasisOrder2();

         Matrix9_1 polynomialEigen = Matrix9_1::Zero();

         for (int j = 0; j < 9; j++)
         {
           polynomialEigen(j) = polynomial[j];

         }
         Matrix3_1 Teigen = W*polynomialEigen;

////         Matrix3_1 TeigenPrime = WPrime*polynomialEigen;  // [P(X+delX) - P(X)] * d^T * K_inverse * PolyBase < (X+delX) - X >

         Vector3D T(0.0);

         for (int k = 0; k < 3; k++)
         {
           T[k] = Teigen(k);
         }
         T = T*influence;


////         Vector3D TPrime(0.0);

////         for (int k = 0; k < 3; k++)
////         {
////           TPrime[k] = TeigenPrime(k);
////         }
////         TPrime = TPrime*influence; // influence * [P(X+delX) - P(X)] * d^T * K_inverse * PolyBase < (X+delX) - X > 

         bond->setForceVectorState(T);

////         bond->setFirstPiolaDivergenceIntegrand(TPrime);
                 
      }
    }
    else
    {
      bond->setForceVectorState(zero);
      bond->setBrokenBond(true);
      bond->setInfluenceFunction(0.0);
//      StateBondP twin_bond = twinBond(bond);
//      twin_bond->setBrokenBond(true);
//      twin_bond->setForceVectorState(zero);
//      twin_bond->setInfluenceFunction(0.0);

    }

//    if (getID() == 700)
//    {
//      std::cout << bond->getCenterPoint()->getID() << ", "
//                << bond->getSecondPoint()->getID() << std::endl;
//      std::cout << "sai= " << bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      std::cout << " BondForceState= " << bond->getForceVectorState() << std::endl;    
//    }

//    bool condition1 = ((bond->getCenterPoint()->getID() == 700) &&      
//                      (bond->getSecondPoint()->getID() == 613));

//    bool condition2 = ((bond->getCenterPoint()->getID() == 613) &&      
//                      (bond->getSecondPoint()->getID() == 700));
//    if (condition1)
//    {
//      printScreenMatrix3D(getShapeTensor(), "ShapeTensor");
//      printScreenMatrix3D(invMatrix3D(getShapeTensor()), "ShapeTensorInverse");
//      printScreenMatrix3D(getFirstPiolaKirchhoffStress()*invMatrix3D(getShapeTensor()), "PKinverse");
//      std::cout << bond->getCenterPoint()->getID() << ", "
//                << bond->getSecondPoint()->getID() << std::endl;
//      std::cout << "sai= " << bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      std::cout << " BondForceState= " << bond->getForceVectorState() << std::endl;    
//    }
    
//    if (condition2)
//    {
//      printScreenMatrix3D(getShapeTensor(), "ShapeTensor");
//      printScreenMatrix3D(invMatrix3D(getShapeTensor()), "ShapeTensorInverse");
//      printScreenMatrix3D(getFirstPiolaKirchhoffStress()*invMatrix3D(getShapeTensor()), "PKinverse");
//      std::cout << bond->getCenterPoint()->getID() << ", "
//                << bond->getSecondPoint()->getID() << std::endl;
//      std::cout << "sai= " << bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      std::cout << " BondForceState= " << bond->getForceVectorState() << std::endl;    
//    }

  }


}



void
StateMaterialPoint::forceVectorStateOfBonds(const double& constant)
{
//  std::cout << "ID= " << getID() << "constant= " << constant << std::endl;
  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  int i = 0;
  double influence = 0.0;
  Vector3D sai(0.0, 0.0, 0.0);
  Vector3D zero(0.0, 0.0, 0.0);
  Matrix3D matrix(0.0);

  Vector3D eta(0.0, 0.0, 0.0);

  double volume = 0.0; 
  
  Vector3D sum(0.0, 0.0, 0.0);
  Vector3D averageDisplacementState(0.0, 0.0, 0.0);

  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    if (bond->getBrokenBond() == false)
    {
      influence = bond->getInfluenceFunction();
//      std::cout << "(" << bond->getCenterPoint()->getID() << "," << bond->getSecondPoint()->getID() << ")"; //<< " influence= " << influence << std::endl;
      sai = bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      std::cout << " sai= " << sai << std::endl;
      bond->calculateEta();
      eta = bond->getEta();
//      if (eta != zero)
//      {
//        std::cout << "(" << bond->getCenterPoint()->getID() << "," << bond->getSecondPoint()->getID() << ")";
//        std::cout << " eta= " << eta << std::endl;
//      }
      volume = bond->getVolume();
      sum = sum + eta*(influence*volume);
    }
  }

  averageDisplacementState = sum*constant;

//  std::cout << "ID= " << getID() << " averageDisplacementState= " << averageDisplacementState << std::endl;

 
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    if (bond->getBrokenBond() == false)
    {
      influence = bond->getInfluenceFunction();
      sai = bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
      eta = bond->getEta();
//      matrix = getFirstPiolaKirchhoffStress()*getShapeTensor().Inverse();
      if (getShapeTensor().Determinant() == 0)
      {
        std::cout << "Hazard 6: Node " << getID() << " Shape Tensor has a zero determinant." << std::endl;
      }

      matrix = getFirstPiolaKirchhoffStress()*invMatrix3D(getShapeTensor());  
      bond->setForceVectorState((matrix*sai)*influence);
      bond->setForceVectorState(bond->getForceVectorState() + averageDisplacementState);
    }
    else
    {
      bond->setForceVectorState(zero);
      bond->setBrokenBond(true);
      bond->setInfluenceFunction(0.0);
//      StateBondP twin_bond = twinBond(bond);
//      twin_bond->setBrokenBond(true);
//      twin_bond->setForceVectorState(zero);
//      twin_bond->setInfluenceFunction(0.0);

    }

//    if (getID() == 700)
//    {
//      std::cout << bond->getCenterPoint()->getID() << ", "
//                << bond->getSecondPoint()->getID() << std::endl;
//      std::cout << "sai= " << bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      std::cout << " BondForceState= " << bond->getForceVectorState() << std::endl;    
//    }

//    bool condition1 = ((bond->getCenterPoint()->getID() == 700) &&      
//                      (bond->getSecondPoint()->getID() == 613));

//    bool condition2 = ((bond->getCenterPoint()->getID() == 613) &&      
//                      (bond->getSecondPoint()->getID() == 700));
//    if (condition1)
//    {
//      printScreenMatrix3D(getShapeTensor(), "ShapeTensor");
//      printScreenMatrix3D(invMatrix3D(getShapeTensor()), "ShapeTensorInverse");
//      printScreenMatrix3D(getFirstPiolaKirchhoffStress()*invMatrix3D(getShapeTensor()), "PKinverse");
//      std::cout << bond->getCenterPoint()->getID() << ", "
//                << bond->getSecondPoint()->getID() << std::endl;
//      std::cout << "sai= " << bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      std::cout << " BondForceState= " << bond->getForceVectorState() << std::endl;    
//    }
    
//    if (condition2)
//    {
//      printScreenMatrix3D(getShapeTensor(), "ShapeTensor");
//      printScreenMatrix3D(invMatrix3D(getShapeTensor()), "ShapeTensorInverse");
//      printScreenMatrix3D(getFirstPiolaKirchhoffStress()*invMatrix3D(getShapeTensor()), "PKinverse");
//      std::cout << bond->getCenterPoint()->getID() << ", "
//                << bond->getSecondPoint()->getID() << std::endl;
//      std::cout << "sai= " << bond->getSecondPoint()->getPosOld() - bond->getCenterPoint()->getPosOld();
//      std::cout << " BondForceState= " << bond->getForceVectorState() << std::endl;    
//    }

  }


}




void
StateMaterialPoint::calculateForceVectorStateOfBonds(const double& delT,
                                                     const double& firstLamme,
                                                     const double& secondLamme,
                                                     const double& yield)
{

//  if (getID() == 700)
//  {
//    std::cout << "delT= " << delT << " firstLamme= " << firstLamme 
//                          << std::endl << "secondLamme= " << secondLamme
//                          << " criticalEquivalentStrian= " << yield << std::endl; 
//  } 

////  shapeTensor();

//  if (getID() == 700)
//  {
//  printScreenMatrix3D(getShapeTensor(), "shapeTensor");
//  printScreenMatrix3D(getShapeTensor().Inverse(), "shapeTensorInverse");
//  printScreenMatrix3D(getShapeTensor()*getShapeTensor().Inverse(), "Check Identity");
//  printScreenMatrix3D(getShapeTensor().Inverse()*getShapeTensor(), "Check Identity");

//  printScreenMatrix3D(invMatrix3D(getShapeTensor()), "shapeTensorInverse");
//  printScreenMatrix3D(getShapeTensor()*invMatrix3D(getShapeTensor()), "Check Identity");
//  printScreenMatrix3D(invMatrix3D(getShapeTensor())*getShapeTensor(), "Check Identity");
//  }

////  deformationGradient();
//  if (getID() == 700)
//    printScreenMatrix3D(getDeformationGradient(), "deformationGradient");


////  LagrangianStrainTensor();
//  if (getID() == 700)
//    printScreenMatrix3D(getLagrangianStrainTensor(), "LagrangianStrainTensor");
//  std::cout << "After calculating LagrangianStrainTensor of point number " << getID() << std::endl;

 
////  checkBrokenBonds(yield);  // Damage modeling, compare equivalent strain with critical equivalent strain
//  std::cout << "After checking broken bonds of point number " << getID() << std::endl;

 
////  calculateStateLocalDamage();  // find the local damage
//  std::cout << "After calculating StateLocalDamage of point number " << getID() << std::endl;
 

  spatialVelocityGradient();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getDeformationGradient(), "DeformationGradient");
//    printScreenMatrix3D(invMatrix3D(getDeformationGradient()), "inverseOfDeformationGradient");
//    printScreenMatrix3D(invMatrix3D(getDeformationGradient())*getDeformationGradient(), "Identity");
//    printScreenMatrix3D(materialTimeDerivative(), "materialTimeDerivativeOfDeformationGradientTensor");
//    printScreenMatrix3D(getVelocityGradient(), "spatialVelocityGradient");
//    printScreenMatrix3D(materialTimeDerivative()*invMatrix3D(getDeformationGradient()), "spatialVelocityGradient");
//  }
//  std::cout << "After calculating spatialVelocityGradient of point number " << getID() << std::endl;


  deformationRateTensor();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getVelocityGradient(), "spatialVelocityGradient");
//    printScreenMatrix3D(getVelocityGradient().Transpose(), "spatialVelocityGradientTranspose");
//    printScreenMatrix3D(getDeformationRateTensor(), "deformationRateTensor");

//  }
//  std::cout << "After calculating deformationRateTensor of point number " << getID() << std::endl;
 

  spinTensor();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getSpinTensor(), "spinTensor");
//    printScreenMatrix3D(getSpinTensor() + getDeformationRateTensor());
//  }
//  std::cout << "After calculating spinTensor of point number " << getID() << std::endl; 
////  polarDecompositionOfDeformationGradient();


  updateRotationTensorMoreStable(delT);


//  updateRotationTensorMoreEfficient(delT);
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getRateRigidBodyRotation(), "omega");
//    printScreenMatrix3D(getRotationTensor(), "rotationTensor");
//    printScreenMatrix3D(getRotationTensor()*getRotationTensor().Transpose(), "checkOrthogonality"); 
//  }   
//  std::cout << "After updating RotationTensorMoreEfficient of point number " << getID() << std::endl;

 
  leftStretchRateTensor();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getLeftStretchRate(), "LeftStretchRate");
//  }     
//  std::cout << "After calculating leftStretchRateTensor of point number " << getID() << std::endl;

 
  updateLeftStretch(delT);
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getLeftStretch(), "LeftStretch");
//  }     
//  std::cout << "After calculating updateLeftStretch of point number " << getID() << std::endl;


  unrotatedDeformationRate();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getUnrotatedDeformationGradient(), "UnrotatedDeformationGradient");
//    printScreenMatrix3D(getDeformationRateTensor(), "DeformationRateTensor");
//    printScreenMatrix3D(getRotationTensor(), "RotationTensor");
//  }     
//  std::cout << "After calculating unrotatedDeformationRate of point number " << getID() << std::endl;


  strainIncrement(delT);
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getTotalElasticStrianIncrement(), "TotalElasticStrianIncrement");
//    printScreenMatrix3D(getDeviatoricStrianIncrement(), "DeviatoricStrianIncrement");
//  } 
//  std::cout << "After calculating strainIncrement of point number " << getID() << std::endl;


  updateUnrotatedCauchyStress(firstLamme, secondLamme);
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getUnrotatedCauchyStress(), "UnrotatedCauchyStress");
//  }
//  std::cout << "After calculating updateUnrotatedCauchyStress of point number " << getID() << std::endl;


  rotatedCauchyStress();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getRotatedCauchyStress(), "RotatedCauchyStress");
//  }
//  std::cout << "After calculating rotatedCauchyStress of point number " << getID() << std::endl;



  firstPiolaKirchhoffStress();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getFirstPiolaKirchhoffStress(), "FirstPiolaKirchhoffStress");
//  }
//  std::cout << "After calculating firstPiolaKirchhoffStress of point number " << getID() << std::endl;


  forceVectorStateOfBonds();
//  std::cout << "After calculating forceVectorStateOfBonds of point number " << getID() << std::endl;

}



void
StateMaterialPoint::calculateForceVectorStateOfBonds(const double& delT,
                                                     const double& firstLamme,
                                                     const double& secondLamme,
                                                     const double& yield,
                                                     const int& numIter, 
                                                     const Matrix3D& initialVelocityGradient)
{

//  if (getID() == 700)
//  {
//    std::cout << "delT= " << delT << " firstLamme= " << firstLamme 
//                          << std::endl << "secondLamme= " << secondLamme
//                          << " criticalEquivalentStrian= " << yield << std::endl; 
//  } 

////  shapeTensor();

// if (getID() == 15604)
// {
//   printScreenMatrix3D(getShapeTensor(), "shapeTensor");
//   printScreenMatrix3D(getShapeTensor().Inverse(), "shapeTensorInverse");
//   printScreenMatrix3D(getShapeTensor()*getShapeTensor().Inverse(), "Check Identity");
//   printScreenMatrix3D(getShapeTensor().Inverse()*getShapeTensor(), "Check Identity");

//   printScreenMatrix3D(invMatrix3D(getShapeTensor()), "shapeTensorInverse");
//   printScreenMatrix3D(getShapeTensor()*invMatrix3D(getShapeTensor()), "Check Identity");
//   printScreenMatrix3D(invMatrix3D(getShapeTensor())*getShapeTensor(), "Check Identity");
// }

////  deformationGradient();
//  if (getID() == 15604)
//  {
//    printScreenMatrix3D(getDeformationGradient(), "deformationGradient");
//  }

////  LagrangianStrainTensor();
//  if (getID() == 700)
//    printScreenMatrix3D(getLagrangianStrainTensor(), "LagrangianStrainTensor");
//  std::cout << "After calculating LagrangianStrainTensor of point number " << getID() << std::endl;

 
////  checkBrokenBonds(yield);  // Damage modeling, compare equivalent strain with critical equivalent strain
//  std::cout << "After checking broken bonds of point number " << getID() << std::endl;

 
////  calculateStateLocalDamatwinBondge();  // find the local damage
//  std::cout << "After calculating StateLocalDamage of point number " << getID() << std::endl;
 

  spatialVelocityGradient(numIter, initialVelocityGradient);  // L = (Fdot)(Finverce)


//  if (getID() == 15604)
//  {
//    printScreenMatrix3D(getDeformationGradient(), "DeformationGradient");
//    printScreenMatrix3D(invMatrix3D(getDeformationGradient()), "inverseOfDeformationGradient");
//    printScreenMatrix3D(invMatrix3D(getDeformationGradient())*getDeformationGradient(), "Identity");
//    printScreenMatrix3D(materialTimeDerivative(), "materialTimeDerivativeOfDeformationGradientTensor");
//    printScreenMatrix3D(getVelocityGradient(), "spatialVelocityGradient");
//    printScreenMatrix3D(materialTimeDerivative()*invMatrix3D(getDeformationGradient()), "spatialVelocityGradient");
//  }
//  std::cout << "After calculating spatialVelocityGradient of point number " << getID() << std::endl;


  deformationRateTensor();  //  D = sym(L) = (L + Ltranspose)/2


//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getVelocityGradient(), "spatialVelocityGradient");
//    printScreenMatrix3D(getVelocityGradient().Transpose(), "spatialVelocityGradientTranspose");
//    printScreenMatrix3D(getDeformationRateTensor(), "deformationRateTensor");

//  }
//  std::cout << "After calculating deformationRateTensor of point number " << getID() << std::endl;
 

  spinTensor();   //  W = skew(L) = (L - Ltranspose)/2          "Antisymmetric tensor"


//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getSpinTensor(), "spinTensor");
//    printScreenMatrix3D(getSpinTensor() + getDeformationRateTensor());
//  }
//  std::cout << "After calculating spinTensor of point number " << getID() << std::endl; 
////  polarDecompositionOfDeformationGradient();


  updateRotationTensorMoreStable(delT);  // Rubinstein and Atluri (1983)


//  updateRotationTensorMoreEfficient(delT);  // Hughes and Winget (1980) 
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getRateRigidBodyRotation(), "omega");
//    printScreenMatrix3D(getRotationTensor(), "rotationTensor");
//    printScreenMatrix3D(getRotationTensor()*getRotationTensor().Transpose(), "checkOrthogonality"); 
//  }   
//  std::cout << "After updating RotationTensorMoreEfficient of point number " << getID() << std::endl;

 
  leftStretchRateTensor();    //   Vdot = (L)(V) - (V)(omega)


//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getLeftStretchRate(), "LeftStretchRate");
//  }     
//  std::cout << "After calculating leftStretchRateTensor of point number " << getID() << std::endl;

 
  updateLeftStretch(delT);    //  V(t+delT) = V(t) + delT*Vdot


//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getLeftStretch(), "LeftStretch");
//  }     
//  std::cout << "After calculating updateLeftStretch of point number " << getID() << std::endl;


  unrotatedDeformationRate();   //  d = (Rtranspose)(D)(R)


//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getUnrotatedDeformationGradient(), "UnrotatedDeformationGradient");
//    printScreenMatrix3D(getDeformationRateTensor(), "DeformationRateTensor");
//    printScreenMatrix3D(getRotationTensor(), "RotationTensor");
//  }     
//  std::cout << "After calculating unrotatedDeformationRate of point number " << getID() << std::endl;


  strainIncrement(delT);   //  delE = d(delT), delEdev = delE - (1/3)(trace(delE))(I)


//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getTotalElasticStrianIncrement(), "TotalElasticStrianIncrement");
//    printScreenMatrix3D(getDeviatoricStrianIncrement(), "DeviatoricStrianIncrement");
//  } 
//  std::cout << "After calculating strainIncrement of point number " << getID() << std::endl;


  updateUnrotatedCauchyStress(firstLamme, secondLamme);  // taw(t+delT) = taw(t) + (firstLamme)(trace(delE))I + 2(secondLamme)(delEdev)


//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getUnrotatedCauchyStress(), "UnrotatedCauchyStress");
//  }
//  std::cout << "After calculating updateUnrotatedCauchyStress of point number " << getID() << std::endl;


  rotatedCauchyStress();  //  sigma = (R)(taw)(Rtranspose)


//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getRotatedCauchyStress(), "RotatedCauchyStress");
//  }
//  std::cout << "After calculating rotatedCauchyStress of point number " << getID() << std::endl;



  firstPiolaKirchhoffStress();  // P = (J)(sigma)(FminusTtranspose)


//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getFirstPiolaKirchhoffStress(), "FirstPiolaKirchhoffStress");
//  }
//  std::cout << "After calculating firstPiolaKirchhoffStress of point number " << getID() << std::endl;


  forceVectorStateOfBonds();  //  T<sai> = (influenceFunction<sai>)(P)(Kinverse)(sai)


//  std::cout << "After calculating forceVectorStateOfBonds of point number " << getID() << std::endl;

}



void
StateMaterialPoint::calculateForceVectorStateOfBonds(const double& delT,
                                                     const double& firstLamme,
                                                     const double& secondLamme,
                                                     const double& yield,
                                                     const int& numIter, 
                                                     const Matrix3D& initialVelocityGradient,
                                                     const double& constant)
{

//  if (getID() == 700)
//  {
//    std::cout << "delT= " << delT << " firstLamme= " << firstLamme 
//                          << std::endl << "secondLamme= " << secondLamme
//                          << " criticalEquivalentStrian= " << yield << std::endl; 
//  } 

////  shapeTensor();

// if (getID() == 15604)
// {
//   printScreenMatrix3D(getShapeTensor(), "shapeTensor");
//   printScreenMatrix3D(getShapeTensor().Inverse(), "shapeTensorInverse");
//   printScreenMatrix3D(getShapeTensor()*getShapeTensor().Inverse(), "Check Identity");
//   printScreenMatrix3D(getShapeTensor().Inverse()*getShapeTensor(), "Check Identity");

//   printScreenMatrix3D(invMatrix3D(getShapeTensor()), "shapeTensorInverse");
//   printScreenMatrix3D(getShapeTensor()*invMatrix3D(getShapeTensor()), "Check Identity");
//   printScreenMatrix3D(invMatrix3D(getShapeTensor())*getShapeTensor(), "Check Identity");
// }

////  deformationGradient();
//  if (getID() == 15604)
//  {
//    printScreenMatrix3D(getDeformationGradient(), "deformationGradient");
//  }

////  LagrangianStrainTensor();
//  if (getID() == 700)
//    printScreenMatrix3D(getLagrangianStrainTensor(), "LagrangianStrainTensor");
//  std::cout << "After calculating LagrangianStrainTensor of point number " << getID() << std::endl;

 
////  checkBrokenBonds(yield);  // Damage modeling, compare equivalent strain with critical equivalent strain
//  std::cout << "After checking broken bonds of point number " << getID() << std::endl;

 
////  calculateStateLocalDamage();  // find the local damage
//  std::cout << "After calculating StateLocalDamage of point number " << getID() << std::endl;
 

  spatialVelocityGradient(numIter, initialVelocityGradient);
//  if (getID() == 15604)
//  {
//    printScreenMatrix3D(getDeformationGradient(), "DeformationGradient");
//    printScreenMatrix3D(invMatrix3D(getDeformationGradient()), "inverseOfDeformationGradient");
//    printScreenMatrix3D(invMatrix3D(getDeformationGradient())*getDeformationGradient(), "Identity");
//    printScreenMatrix3D(materialTimeDerivative(), "materialTimeDerivativeOfDeformationGradientTensor");
//    printScreenMatrix3D(getVelocityGradient(), "spatialVelocityGradient");
//    printScreenMatrix3D(materialTimeDerivative()*invMatrix3D(getDeformationGradient()), "spatialVelocityGradient");
//  }
//  std::cout << "After calculating spatialVelocityGradient of point number " << getID() << std::endl;


  deformationRateTensor();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getVelocityGradient(), "spatialVelocityGradient");
//    printScreenMatrix3D(getVelocityGradient().Transpose(), "spatialVelocityGradientTranspose");
//    printScreenMatrix3D(getDeformationRateTensor(), "deformationRateTensor");

//  }
//  std::cout << "After calculating deformationRateTensor of point number " << getID() << std::endl;
 

  spinTensor();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getSpinTensor(), "spinTensor");
//    printScreenMatrix3D(getSpinTensor() + getDeformationRateTensor());
//  }
//  std::cout << "After calculating spinTensor of point number " << getID() << std::endl; 
////  polarDecompositionOfDeformationGradient();


  updateRotationTensorMoreStable(delT);


//  updateRotationTensorMoreEfficient(delT);
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getRateRigidBodyRotation(), "omega");
//    printScreenMatrix3D(getRotationTensor(), "rotationTensor");
//    printScreenMatrix3D(getRotationTensor()*getRotationTensor().Transpose(), "checkOrthogonality"); 
//  }   
//  std::cout << "After updating RotationTensorMoreEfficient of point number " << getID() << std::endl;

 
  leftStretchRateTensor();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getLeftStretchRate(), "LeftStretchRate");
//  }     
//  std::cout << "After calculating leftStretchRateTensor of point number " << getID() << std::endl;

 
  updateLeftStretch(delT);
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getLeftStretch(), "LeftStretch");
//  }     
//  std::cout << "After calculating updateLeftStretch of point number " << getID() << std::endl;


  unrotatedDeformationRate();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getUnrotatedDeformationGradient(), "UnrotatedDeformationGradient");
//    printScreenMatrix3D(getDeformationRateTensor(), "DeformationRateTensor");
//    printScreenMatrix3D(getRotationTensor(), "RotationTensor");
//  }     
//  std::cout << "After calculating unrotatedDeformationRate of point number " << getID() << std::endl;


  strainIncrement(delT);
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getTotalElasticStrianIncrement(), "TotalElasticStrianIncrement");
//    printScreenMatrix3D(getDeviatoricStrianIncrement(), "DeviatoricStrianIncrement");
//  } 
//  std::cout << "After calculating strainIncrement of point number " << getID() << std::endl;


  updateUnrotatedCauchyStress(firstLamme, secondLamme);
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getUnrotatedCauchyStress(), "twinBondUnrotatedCauchyStress");
//  }
//  std::cout << "After calculating updateUnrotatedCauchyStress of point number " << getID() << std::endl;


  rotatedCauchyStress();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getRotatedCauchyStress(), "RotatedCauchyStress");
//  }
//  std::cout << "After calculating rotatedCauchyStress of point number " << getID() << std::endl;



  firstPiolaKirchhoffStress();
//  if (getID() == 700)
//  {
//    printScreenMatrix3D(getFirstPiolaKirchhoffStress(), "FirstPiolaKirchhoffStress");
//  }
//  std::cout << "After calculating firstPiolaKirchhoffStress of point number " << getID() << std::endl;


  forceVectorStateOfBonds(constant);
//  std::cout << "After calculating forceVectorStateOfBonds of point number " << getID() << std::endl;

}




void
StateMaterialPoint::calculateForceVectorStateOfBonds(const double& delT,
                                                     const double& firstLamme,
                                                     const double& secondLamme,
                                                     const double& yield,
                                                     const int& numIter, 
                                                     const Matrix3D& initialVelocityGradient,
                                                     const int& orderFlag)
{
  setOrderFlag(orderFlag);
//  calculateForceVectorStateOfBonds(delT, firstLamme, secondLamme, yield, numIter, initialVelocityGradient, orderFlag);
  deformationGradientUsingPolynomials();
  spatialVelocityGradient(numIter, initialVelocityGradient);
  deformationRateTensor();
  spinTensor();
  updateRotationTensorMoreStable(delT);
  leftStretchRateTensor();
  updateLeftStretch(delT); 
  unrotatedDeformationRate();
  strainIncrement(delT);
  updateUnrotatedCauchyStress(firstLamme, secondLamme);
  rotatedCauchyStress();
  firstPiolaKirchhoffStress();
  forceVectorStateOfBonds();
//  forceVectorStateOfBonds(constant);  

}



void
StateMaterialPoint::calculateFirstPiolaKirchhoffStress(const double& delT,
                                                       const double& firstLamme,
                                                       const double& secondLamme,
                                                       const double& yield,
                                                       const int& numIter, 
                                                       const Matrix3D& initialVelocityGradient,
                                                       const int& orderFlag)
{
  setOrderFlag(orderFlag);
//  calculateForceVectorStateOfBonds(delT, firstLamme, secondLamme, yield, numIter, initialVelocityGradient, orderFlag);
  deformationGradientUsingPolynomials();
  spatialVelocityGradient(numIter, initialVelocityGradient);
  deformationRateTensor();
  spinTensor();
  updateRotationTensorMoreStable(delT);
  leftStretchRateTensor();
  updateLeftStretch(delT); 
  unrotatedDeformationRate();
  strainIncrement(delT);
  updateUnrotatedCauchyStress(firstLamme, secondLamme);
  rotatedCauchyStress();
  firstPiolaKirchhoffStress();
//  forceVectorStateOfBonds();
//  forceVectorStateOfBonds(constant);  

}





void
StateMaterialPoint::setForceBondsZero()
{
  StateBondPArray bonds = getStateBonds();
  int size = bonds.size();
  int i = 0;
//  double influence = 0.0;
//  Vector3D sai(0.0, 0.0, 0.0);
  Vector3D zero(0.0, 0.0, 0.0);
//  Matrix3D matrix(0.0);  
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bonds[i];
    bond->setForceVectorState(zero);
    bond->setBrokenBond(true);
  }

}


/*void
StateMaterialPoint::LagrangianStrainTensor()
{
//  shapeTensor();
//  deformationGradient();

  deformationGradientUsingPolynomials();

  Matrix3D F = getDeformationGradient();
  Matrix3D I;
  I.Identity(); 
  setLagrangianStrainTensor((F.Transpose()*F - I)*(0.5));
    
}*/



void 
StateMaterialPoint::calculateStateLocalDamage()
{
  StateBondPArray bondsArray = getStateBonds();
  int size = bondsArray.size();
  double numerator = 0.0;
  double denominator = 0.0;
  int i = 0;
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bondsArray[i];
    StateMaterialPointP secondPoint = bond->getSecondPoint();
////    double volume = secondPoint->getInitialVolume();

    double volume = bond->getVolume();

    if (!bond->getBrokenBond())
    {
      numerator += volume;
    }
    denominator += volume;
    
  }
  double damage = numerator/denominator;
  setLocalDamage(damage);
//  if (getID() == 700)
//    std::cout << "Damage of point number " << getID() << " is: " << getLocalDamage() << std::endl;
  
}


void
StateMaterialPoint::checkBrokenBonds(const double& yield)
{
  Vector3D zero(0.0, 0.0, 0.0);
  StateBondPArray bondsArray = getStateBonds();
  int size = bondsArray.size();
//  double numerator = 0.0;
//  double denominator = 0.0;
  int i = 0;
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bondsArray[i];
    bond->equivalentStrain();

//    if (getID() == 700)
//    {
//      std::cout << "Equivalent strain of bond " << bond->getCenterPoint()->getID() 
//                << ", " << bond->getSecondPoint()->getID() << " is: "
//                << bond->getEquivalentStrain() << std::endl;
//      std::cout << "Critical equivalent strain is: " << yield << std::endl;  
//    }

    if (bond->getEquivalentStrain() > yield)
    {
      bond->setBrokenBond(true);
      bond->setForceVectorState(zero);
      bond->setInfluenceFunction(0.0);
//      StateBondP twin_bond = twinBond(bond);
//      twin_bond->setBrokenBond(true);
//      twin_bond->setForceVectorState(zero);
//      twin_bond->setInfluenceFunction(0.0);

    }
  }

}



void
StateMaterialPoint::checkBrokenBonds(const double& yield, bool& flag)
{

  Vector3D zero(0.0, 0.0, 0.0);
  StateBondPArray bondsArray = getStateBonds();
  int size = bondsArray.size();
//  double numerator = 0.0;
//  double denominator = 0.0;
  int i = 0;
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bondsArray[i];
    if (bond->getBrokenBond() == false)
    {
      bond->equivalentStrain();

  //    if (getID() == 700)
  //    {
  //      std::cout << "Equivalent strain of bond " << bond->getCenterPoint()->getID() 
  //                << ", " << bond->getSecondPoint()->getID() << " is: "
  //                << bond->getEquivalentStrain() << std::endl;
  //      std::cout << "Critical equivalent strain is: " << yield << std::endl;  
  //    }

      if (bond->getEquivalentStrain() > yield)
      {
        bond->setBrokenBond(true);
        bond->setForceVectorState(zero);
        bond->setInfluenceFunction(0.0);
        StateBondP twin_bond = twinBond(bond);
//        twin_bond->setBrokenBond(true);
//        twin_bond->setForceVectorState(zero);
//        twin_bond->setInfluenceFunction(0.0);


        flag = true;
      } 
    }
//    else
//    {
//      flag = false;
//    }

  }

}





/*void
StateMaterialPoint::checkBrokenBonds(const double& yield, bool& flag)
{
  Vector3D zero(0.0, 0.0, 0.0);
  StateBondPArray bondsArray = getStateBonds();
  int size = bondsArray.size();
//  double numerator = 0.0;
//  double denominator = 0.0;
  int i = 0;
  for (i = 0; i < size; i++)
  {
    StateBondP bond = bondsArray[i];
    if (bond->getBrokenBond() == false)
    {
      bond->equivalentStrain();

  //    if (getID() == 700)
  //    {
  //      std::cout << "Equivalent strain of bond " << bond->getCenterPoint()->getID() 
  //                << ", " << bond->getSecondPoint()->getID() << " is: "
  //                << bond->getEquivalentStrain() << std::endl;
  //      std::cout << "Critical equivalent strain is: " << yield << std::endl;  
  //    }

      if (bond->getEquivalentStrain() > yield)
      {
        bond->setBrokenBond(true);
        bond->setForceVectorState(zero);
        bond->setInfluenceFunction(0.0);
        StateBondP twin_bond = twinBond(bond);
        twin_bond->setBrokenBond(true);
        twin_bond->setForceVectorState(zero);
        twin_bond->setInfluenceFunction(0.0);


        flag = true;
      } 
    }
//    else
//    {
//      flag = false;
//    }

  }

}*/







void
StateMaterialPoint::printScreenMatrix3D(const Matrix3D& matrix)
{
  std::cout << getID() << std::endl;
  std::cout << matrix(0,0) << " " << matrix(0,1) << " " << matrix(0,2) << std::endl;
  std::cout << matrix(1,0) << " " << matrix(1,1) << " " << matrix(1,2) << std::endl;
  std::cout << matrix(2,0) << " " << matrix(2,1) << " " << matrix(2,2) << std::endl;

}


void
StateMaterialPoint::printScreenMatrix3D(const Matrix3D& matrix, const std::string& name)
{
  std::cout << "The " << name << " of point number " << getID() << ":" << std::endl;
  std::cout << matrix(0,0) << " " << matrix(0,1) << " " << matrix(0,2) << std::endl;
  std::cout << matrix(1,0) << " " << matrix(1,1) << " " << matrix(1,2) << std::endl;
  std::cout << matrix(2,0) << " " << matrix(2,1) << " " << matrix(2,2) << std::endl << std::endl;

}



Matrix3D
StateMaterialPoint::invMatrix3D(const Matrix3D& matrix)
{
  MatrixXd mat(3,3);
  mat = matrix3dToEigenMatrix(matrix);
  MatrixXd inv = mat.inverse();
  return eigenMatrixToMatrix3D(inv); 

}




double 
StateMaterialPoint::ramp(const int& num_iter, const int& num_initial)
{
  
  if (num_initial == 1)
  {
    return 1.0;
  }
  else
  {
    double denominator = (double) (num_initial - 1);    
    double coeff =  1/denominator;
//    if (getID() == 4224)
//    {
//      std::cout << "coeff= " << coeff << std::endl;
//    } 
    return coeff*num_iter;
  }


}




double
StateMaterialPoint::triangle(const int& num_iter, const int& num_initial)
{
  double result = 0.0;
  if (num_initial == 1)
  {
    result = 1.0;
  }
  else
  {
    double denominator = (double) (num_initial - 1);    
    double coeff =  (2.0)/denominator;
//    if (getID() == 4224)
//    {
//      std::cout << "coeff= " << coeff << std::endl;
//    } 
    double half = denominator/(2.0);
    double num = (double) (num_iter);
    if ((0.0 <= num) && (num <= half))
    {
      result = coeff*num;
    }
    else if ((half <= num) && (num <= denominator))
    {
      result = 2 - coeff*num;
    }
    else
    {
      result = 0.0;
    }
  }

  return result;
  
}








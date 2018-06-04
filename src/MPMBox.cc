#include <Box.h>
#include <MPMBox.h>
#include <NodeP.h>
#include <Node.h>
#include <SymmetricMaterialTrilinearElementSP.h>
#include <SymmetricMaterialTrilinearElement.h>

#include <Tensor.h>

#include<GeometryMath/Vector3D.h>

#include <MaterialPointP.h>  // for MPM
#include <MaterialPoint.h>  // for MPM
#include <math.h>


using namespace FiniteElement;

  MPMBox::MPMBox()
   : d_box(new Box()), d_basis_function(LinearHexahedral)
  {
    d_box_flag = true;
  }


  MPMBox::MPMBox(BoxSP box)
   : d_box(box), d_basis_function(LinearHexahedral)
  {
    d_box_flag = true;
  }

  MPMBox::MPMBox(BoxSP box, basisFunction basis)
   : d_box(box), d_basis_function(basis)
  {
    d_box_flag = true;
  }

  MPMBox::MPMBox(ComplicatedGeometrySP geometry, basisFunction basis)
   : d_geometry(geometry), d_basis_function(basis)
  {
    d_box_flag = false;
      
  }



  BoxSP
  MPMBox::Boxx()
  {
    if (d_box_flag == true)
    {
      return d_box;
    }
    else
    {
      return d_geometry->getBox();
    }

  }  


  std::vector<MaterialPointP>
  MPMBox::getPoints()
  {
    if (d_box_flag == true)
    {
      return d_box->getMaterialPoints();
    }
    else
    {
      return d_geometry->getMaterialPointsArray();
    }

  }  


  double
  MPMBox::interpolate(NodeP node, MaterialPointP material)
  {
    double result = 0.0;
    switch (d_basis_function)
    {
      case LinearHexahedral:
        result = linearHexahedral(node, material);
      break;
      case PiecewiseLinear:
        result = tentFunction(node, material);
      break;
      case GIMPfunction:
        result = tentContiguousParticlesGIMPFunction(node, material);
//        result = tentContiguousParticlesGIMPFunction(node, material, 0.05);
      break;
      case BSpline:
        result = bSplineFunction(node, material);
      break;
      default:
           std::cout << "Four basis functions: LinearHexahedral, PiecewiseLinear, GIMPfunction, BSpline" << std::endl;
      break;
    }
    

    return result;

  }


  Vector3D
  MPMBox::gradientInterpolate(NodeP node, MaterialPointP material)
  {
    Vector3D result(0.0);
    switch (d_basis_function)
    {
      case LinearHexahedral:
        result = gradientLinearHexahedral(node, material);
      break;
      case PiecewiseLinear:
        result = grtadientTentFunction(node, material);
      break;
      case GIMPfunction:
        result = gradientTentContiguousParticlesGIMPFunction(node, material);
//        result = gradientTentContiguousParticlesGIMPFunction(node, material, 0.5);
      break;
      case BSpline:
        result = gradientBSplineFunction(node, material);
      break;
      default:
           std::cout << "Four basis functions: LinearHexahedral, PiecewiseLinear, GIMPfunction, BSpline" << std::endl;
      break;
    }
    

    return result;

  }


  bool
  MPMBox::isPointInElement(MaterialPointP material, SymmetricMaterialTrilinearElementSP element)
  {
    TrilinearVolumeElementSP triElem = element->triElement();
    bool x = ((triElem->getNodeP1()->x() <= material->getPosOld().x()) && 
              (material->getPosOld().x() <= triElem->getNodeP8()->x()));
    bool y = ((triElem->getNodeP1()->y() <= material->getPosOld().y()) && 
              (material->getPosOld().y() <= triElem->getNodeP8()->y()));
    bool z = ((triElem->getNodeP1()->z() <= material->getPosOld().z()) && 
              (material->getPosOld().z() <= triElem->getNodeP8()->z()));


    return ((x && y) && z);
  }


  
  double
  MPMBox::linearHexahedral(NodeP node, MaterialPointP material)
  { 
    double result = 0.0;
  
//    std::cout << "Number of iteration: " << num;
//    std::cout << "   node in interpolation: " << node->getID();

/****    std::vector<MaterialPointP> pointsArray = d_box->getMaterialPoints();
    int pointSize = pointsArray.size();

//  std::cout << pointsArray[i]->getPosOld().x() << ",  "
//            << pointsArray[i]->getPosOld().y() << ",  "
//            << pointsArray[i]->getPosOld().z() << std::endl; 

//  std::cout << pointsArray[i]->getMass() << std::endl; 

    for (int i = 0; i < pointSize; i++)
    {
      MaterialPointP cur_point = pointsArray[i];
      std::cout << "Point: " << cur_point->getID() << std::endl;
      std::cout << cur_point->getPosOld().x() << ",  "
                << cur_point->getPosOld().y() << ",  "
                << cur_point->getPosOld().z() << std::endl; 
    }
//  std::cout << pointsArray[i]->getMass() << std::endl; ****/
    

    bool condition = false;
//    std::vector<SymmetricMaterialTrilinearElementSP> elementsArray = d_box->getElements();
    std::vector<SymmetricMaterialTrilinearElementSP> elementsArray = Boxx()->getElements();
    int elementSize = elementsArray.size();
    for (int jj = 0; jj < elementSize; jj++)
    { 
      SymmetricMaterialTrilinearElementSP cur_element = elementsArray[jj];
//      if ((isPointInElement(material, cur_element)) && (d_box->isNodeInElement(node, cur_element)))
      if ((isPointInElement(material, cur_element)) && (Boxx()->isNodeInElement(node, cur_element)))
      {
//        std::cout << "Point number " << material->getID() << " is in element number " << cur_element->id()
//                  << "." << std::endl;
        TrilinearVolumeElementSP triElem = cur_element->triElement();          
        double xi1 = triElem->symXtoXi1(material->getPosOld().x());
        double xi2 = triElem->symYtoXi2(material->getPosOld().y());
        double xi3 = triElem->symZtoXi3(material->getPosOld().z());

//        std::cout << "Point number " << material->getID() << " (" << xi1 << ", " << xi2 << ", " << xi3 << ")" 
//                  << "------> (" << material->getPosOld().x() << ", " <<
//                                    material->getPosOld().y() << ", " <<
//                                    material->getPosOld().z() << ")" <<std::endl;

        std::array <NodeP, 8> elemNodes = triElem->getArray();
 
        for (int kk = 1; kk < 9; kk++)
        {

//          std::cout << "num= " << num << " kk= " << kk << std::endl; 

          NodeP cur_node = elemNodes[kk-1];
//          std::cout << "The ID of the element node: " << cur_node->getID() << std::endl;          
//          if (node == cur_node)
          if (node->getID() == cur_node->getID())
//          if (node->operator==(cur_node))
          {
            result = triElem->symTestFunction(kk, xi1, xi2, xi3);
//            return result;
//            return triElem->symTestFunction(kk, xi1, xi2, xi3);
            condition = true; 
            break;
          }       
        }
      }
//      std::cout << "  " << condition << std::endl;
      if (condition == true)  break;
    }
    if (condition == false) 
    { 
      result = 0.0;   
//      return 0.0;
    }
//    std::cout << "  result= " << result << std::endl;
   return result;
  }

  
  Vector3D
  MPMBox::gradientLinearHexahedral(NodeP node, MaterialPointP material)
  { 
    bool condition = false;
    Vector3D vec(0.0, 0.0, 0.0);  

//    std::cout << "   node: " << node->getID() << std::endl;
  
//    std::vector<SymmetricMaterialTrilinearElementSP> elementsArray = d_box->getElements();
    std::vector<SymmetricMaterialTrilinearElementSP> elementsArray = Boxx()->getElements();
    int elementSize = elementsArray.size();
    for (int jj = 0; jj < elementSize; jj++)
    { 
      SymmetricMaterialTrilinearElementSP cur_element = elementsArray[jj];
//      if ((isPointInElement(material, cur_element)) && (d_box->isNodeInElement(node, cur_element)))
      if ((isPointInElement(material, cur_element)) && (Boxx()->isNodeInElement(node, cur_element)))
      {
//        std::cout << "Point number " << material->getID() << " is in element number " << cur_element->id()
//                  << "." << std::endl;
        TrilinearVolumeElementSP triElem = cur_element->triElement();          
        double xi1 = triElem->symXtoXi1(material->getPosOld().x());
        double xi2 = triElem->symYtoXi2(material->getPosOld().y());
        double xi3 = triElem->symZtoXi3(material->getPosOld().z());

//        std::cout << "Point number " << material->getID() << " (" << xi1 << ", " << xi2 << ", " << xi3 << ")" 
//                  << "------> (" << material->getPosOld().x() << ", " <<
//                                    material->getPosOld().y() << ", " <<
//                                    material->getPosOld().z() << ")" <<std::endl;



        std::array <NodeP, 8> elemNodes = triElem->getArray();
        for (int kk = 1; kk < 9; kk++)
        {
          NodeP cur_node = elemNodes[kk-1];
//          std::cout << "The ID of the element node: " << cur_node->getID(); // << std::endl;
          if (node == cur_node)
//          if (node->operator==(cur_node))
          {
//            std::cout << "  kk= " << kk; // << std::endl;

//            double xi1Node = triElem->symXtoXi1(cur_node->x());
//            double xi2Node = triElem->symYtoXi2(cur_node->y());
//            double xi3Node = triElem->symZtoXi3(cur_node->z());
//            std::cout << "  The ID of the element node: " << cur_node->getID();
//            std::cout << "  (" << xi1Node << ", " << xi2Node << ", " << xi3Node << ")" << std::endl;


            double xx = triElem->symDerivativeTestFunctionToPosition(kk, 1, xi1, xi2, xi3);
            double yy = triElem->symDerivativeTestFunctionToPosition(kk, 2, xi1, xi2, xi3);
            double zz = triElem->symDerivativeTestFunctionToPosition(kk, 3, xi1, xi2, xi3);
            vec.x(xx); vec.y(yy); vec.z(zz);
//            return vec;
            condition = true;
            break;
          }       

        }
      }
      if (condition == true) break;
    }
    if (condition == false)
    {
      Vector3D zero(0.0, 0.0, 0.0);
      vec = zero;
    }

//    std::cout << "  vec= " << vec << std::endl;


    return vec;


  }
  

  double
  MPMBox::tentFunction(NodeP node, MaterialPointP material)
  {
    double interpulatedValue = 0.0;
//    double xLength = (d_box->getMaxPoint().x() - d_box->getMinPoint().x())/(d_box->getNumElementsX());
//    double yLength = (d_box->getMaxPoint().y() - d_box->getMinPoint().y())/(d_box->getNumElementsY());
//    double zLength = (d_box->getMaxPoint().z() - d_box->getMinPoint().z())/(d_box->getNumElementsZ());

    double xLength = (Boxx()->getMaxPoint().x() - Boxx()->getMinPoint().x())/(Boxx()->getNumElementsX());
    double yLength = (Boxx()->getMaxPoint().y() - Boxx()->getMinPoint().y())/(Boxx()->getNumElementsY());
    double zLength = (Boxx()->getMaxPoint().z() - Boxx()->getMinPoint().z())/(Boxx()->getNumElementsZ());

    double xDiff = material->getPosOld().x() - node->x();
    double yDiff = material->getPosOld().y() - node->y();
    double zDiff = material->getPosOld().z() - node->z();
    interpulatedValue = tent(xDiff, xLength)*tent(yDiff, yLength)*tent(zDiff, zLength);
    return interpulatedValue;

  }


  double
  MPMBox::tent(double diff, double length)
  {
    double result = 0.0;
    if (diff <= (-1)*length)
    {
      result = 0.0;
    }
    else if (((-1)*length < diff) && (diff <= 0.0 )) 
    {
      result = 1 + diff/length;
    }
    else if ((0.0 < diff) && (diff <= length))
    {
      result = 1 - diff/length;
    }
    else
    {
      result = 0.0;
    }
   
    return result;      
  }


  Vector3D
  MPMBox::grtadientTentFunction(NodeP node, MaterialPointP material)
  {
    Vector3D vec;
//    double interpulatedValue = 0.0;
//    double xLength = (d_box->getMaxPoint().x() - d_box->getMinPoint().x())/(d_box->getNumElementsX());
//    double yLength = (d_box->getMaxPoint().y() - d_box->getMinPoint().y())/(d_box->getNumElementsY());
//    double zLength = (d_box->getMaxPoint().z() - d_box->getMinPoint().z())/(d_box->getNumElementsZ());

    double xLength = (Boxx()->getMaxPoint().x() - Boxx()->getMinPoint().x())/(Boxx()->getNumElementsX());
    double yLength = (Boxx()->getMaxPoint().y() - Boxx()->getMinPoint().y())/(Boxx()->getNumElementsY());
    double zLength = (Boxx()->getMaxPoint().z() - Boxx()->getMinPoint().z())/(Boxx()->getNumElementsZ());

    double xDiff = material->getPosOld().x() - node->x();
    double yDiff = material->getPosOld().y() - node->y();
    double zDiff = material->getPosOld().z() - node->z();
    vec.x(gradientTent(xDiff, xLength)*tent(yDiff, yLength)*tent(zDiff, zLength));
    vec.y(gradientTent(yDiff, yLength)*tent(xDiff, xLength)*tent(zDiff, zLength));
    vec.z(gradientTent(zDiff, zLength)*tent(xDiff, xLength)*tent(yDiff, yLength));
    return vec;

  }


  double
  MPMBox::gradientTent(double diff, double length)
  {
    double result = 0.0;
    if (diff < (-1)*length)
    {
      result = 0.0;
    }
    else if (((-1)*length < diff) && (diff < 0.0 )) 
    {
      result = 1/length;
    }
    else if ((0.0 < diff) && (diff < length))
    {
      result = (-1)/length;
    }
    else
    {
      result = 0.0;
    }
   
    return result;     
  }



  double
  MPMBox::tentContiguousParticlesGIMPFunction(NodeP node, MaterialPointP material)
//  MPMBox::tentContiguousParticlesGIMPFunction(NodeP node, MaterialPointP material, double particleHalfSize)
  {
    double oneThird = 1.0/3.0;
    double interpulatedValue = 0.0;
//    double lp = particleHalfSize;
//    double lp = 0.5*pow(material->getIVolume(), 1/3);
//    if (node->getID() == 19)
//    {
//      std::cout << "PointID= " << material->getID() << " Volume= " << material->getInitialVolume() << " Volume^(1/3)= " << std::pow(material->getInitialVolume(), oneThird) << std::endl;
//    }     
    double lp = 0.5*std::pow(material->getInitialVolume(), oneThird);
//    double xLength = (d_box->getMaxPoint().x() - d_box->getMinPoint().x())/(d_box->getNumElementsX());
//    double yLength = (d_box->getMaxPoint().y() - d_box->getMinPoint().y())/(d_box->getNumElementsY());
//    double zLength = (d_box->getMaxPoint().z() - d_box->getMinPoint().z())/(d_box->getNumElementsZ());

    double xLength = (Boxx()->getMaxPoint().x() - Boxx()->getMinPoint().x())/(Boxx()->getNumElementsX());
    double yLength = (Boxx()->getMaxPoint().y() - Boxx()->getMinPoint().y())/(Boxx()->getNumElementsY());
    double zLength = (Boxx()->getMaxPoint().z() - Boxx()->getMinPoint().z())/(Boxx()->getNumElementsZ());

    double xDiff = material->getPosOld().x() - node->x();
    double yDiff = material->getPosOld().y() - node->y();
    double zDiff = material->getPosOld().z() - node->z();
    interpulatedValue = GIMP(xDiff, xLength, lp)*GIMP(yDiff, yLength, lp)*GIMP(zDiff, zLength, lp);
//    if (node->getID() == 19)
//    {
//      std::cout << "PointID= " << material->getID() << " lp= " << lp << " xLength= " << xLength << " xDiff= " << xDiff << std::endl;
//    }
    return interpulatedValue;

  }


  double
  MPMBox::GIMP(double diff, double length, double lp)
  {
    double result = 0.0;
    if (diff <= (-1)*length-lp)
    {
      result = 0.0;
    }
    else if (((-1)*length-lp < diff) && (diff <= (-1)*length+lp)) 
    {
      result = (1/(4*length*lp))*(length+lp+diff)*(length+lp+diff);
    }
    else if (((-1)*length+lp < diff) && (diff <= (-1)*lp))
    {
      result = 1 + diff/length;
    }
    else if (((-1)*lp < diff) && (diff <= lp))
    {
      result = 1 - (1/(2*length*lp))*(diff*diff+lp*lp);
    }
    else if ((lp < diff) && (diff <= length - lp))
    {
      result = 1 - diff/length;

    }
    else if ((length - lp < diff) && (diff <= length + lp))
    {
      result = (1/(4*length*lp))*(length+lp-diff)*(length+lp-diff); 
    }
    else
    {
      result = 0.0;
    }
    
    return result;
    
  }


  Vector3D
  MPMBox::gradientTentContiguousParticlesGIMPFunction(NodeP node, MaterialPointP material)
//  MPMBox::gradientTentContiguousParticlesGIMPFunction(NodeP node, MaterialPointP material, double particleHalfSize)
  {
    Vector3D vec(0.0);
//    double lp = particleHalfSize;
    double lp = 0.5*pow(material->getVolume(), 1/3);
//    double lp = 0.5*pow(material->getInitialVolume(), 1/3);
//    double xLength = (d_box->getMaxPoint().x() - d_box->getMinPoint().x())/(d_box->getNumElementsX());
//    double yLength = (d_box->getMaxPoint().y() - d_box->getMinPoint().y())/(d_box->getNumElementsY());
//    double zLength = (d_box->getMaxPoint().z() - d_box->getMinPoint().z())/(d_box->getNumElementsZ());

    double xLength = (Boxx()->getMaxPoint().x() - Boxx()->getMinPoint().x())/(Boxx()->getNumElementsX());
    double yLength = (Boxx()->getMaxPoint().y() - Boxx()->getMinPoint().y())/(Boxx()->getNumElementsY());
    double zLength = (Boxx()->getMaxPoint().z() - Boxx()->getMinPoint().z())/(Boxx()->getNumElementsZ());

    double xDiff = material->getPosOld().x() - node->x();
    double yDiff = material->getPosOld().y() - node->y();
    double zDiff = material->getPosOld().z() - node->z();
    vec.x(gradientGIMP(xDiff, xLength, lp)*GIMP(yDiff, yLength, lp)*GIMP(zDiff, zLength, lp));
    vec.y(gradientGIMP(yDiff, yLength, lp)*GIMP(xDiff, xLength, lp)*GIMP(zDiff, zLength, lp));
    vec.z(gradientGIMP(zDiff, zLength, lp)*GIMP(xDiff, xLength, lp)*GIMP(yDiff, yLength, lp));
    return vec;

  }


  double
  MPMBox::gradientGIMP(double diff, double length, double lp)
  {
    double result = 0.0;
    if (diff <= (-1)*length-lp)
    {
      result = 0.0;
    }
    else if (((-1)*length-lp < diff) && (diff <= (-1)*length+lp)) 
    {
      result = (1/(2*length*lp))*(length+lp+diff);
    }
    else if (((-1)*length+lp < diff) && (diff <= (-1)*lp))
    {
      result = 1/length;
    }
    else if (((-1)*lp < diff) && (diff <= lp))
    {
      result =  - (1/(length*lp))*diff;
    }
    else if ((lp < diff) && (diff <= length - lp))
    {
      result = -1/length;

    }
    else if ((length - lp < diff) && (diff <= length + lp))
    {
      result = (-1)*(1/(2*length*lp))*(length+lp-diff); 
    }
    else
    {
      result = 0.0;
    }
 
    return result;   
  }


  double
  MPMBox::bSplineFunction(NodeP node, MaterialPointP material)
  {
    double interpulatedValue = 0.0;
//    double xLength = (d_box->getMaxPoint().x() - d_box->getMinPoint().x())/(d_box->getNumElementsX());
//    double yLength = (d_box->getMaxPoint().y() - d_box->getMinPoint().y())/(d_box->getNumElementsY());
//    double zLength = (d_box->getMaxPoint().z() - d_box->getMinPoint().z())/(d_box->getNumElementsZ());

    double xLength = (Boxx()->getMaxPoint().x() - Boxx()->getMinPoint().x())/(Boxx()->getNumElementsX());
    double yLength = (Boxx()->getMaxPoint().y() - Boxx()->getMinPoint().y())/(Boxx()->getNumElementsY());
    double zLength = (Boxx()->getMaxPoint().z() - Boxx()->getMinPoint().z())/(Boxx()->getNumElementsZ());

    double xDiff = material->getPosOld().x() - node->x();
    double yDiff = material->getPosOld().y() - node->y();
    double zDiff = material->getPosOld().z() - node->z();
    interpulatedValue = bSpline(xDiff, xLength)*bSpline(yDiff, yLength)*bSpline(zDiff, zLength);
    return interpulatedValue;

  }


  double
  MPMBox::bSpline(double diff, double length)
  {
    double result = 0.0;
    if ((diff >= (-1.5)*length) && (diff <= (-0.5)*length))
    {
      result = (1/(2*length*length))*(diff*diff) + (3/(2*length))*diff + 9/8;
    }
    else if ((diff >= (-0.5)*length) && (diff <= 0.5*length)) 
    {
      result = 0.75 - ((1/(length*length))*(diff*diff));
    }
    else if ((diff >= 0.5*length) && (diff <= 1.5*length))
    {
      result = (1/(2*length*length))*(diff*diff) - (3/(2*length))*diff + 9/8;
    }
    else
    {
      result = 0.0;
    }
    
    return result;
    
  }


  Vector3D
  MPMBox::gradientBSplineFunction(NodeP node, MaterialPointP material)
  {
    Vector3D vec(0.0);
//    double xLength = (d_box->getMaxPoint().x() - d_box->getMinPoint().x())/(d_box->getNumElementsX());
//    double yLength = (d_box->getMaxPoint().y() - d_box->getMinPoint().y())/(d_box->getNumElementsY());
//    double zLength = (d_box->getMaxPoint().z() - d_box->getMinPoint().z())/(d_box->getNumElementsZ());

    double xLength = (Boxx()->getMaxPoint().x() - Boxx()->getMinPoint().x())/(Boxx()->getNumElementsX());
    double yLength = (Boxx()->getMaxPoint().y() - Boxx()->getMinPoint().y())/(Boxx()->getNumElementsY());
    double zLength = (Boxx()->getMaxPoint().z() - Boxx()->getMinPoint().z())/(Boxx()->getNumElementsZ());

    double xDiff = material->getPosOld().x() - node->x();
    double yDiff = material->getPosOld().y() - node->y();
    double zDiff = material->getPosOld().z() - node->z();
    vec.x(gradientBSpline(xDiff, xLength)*bSpline(yDiff, yLength)*bSpline(zDiff, zLength));
    vec.y(gradientBSpline(yDiff, yLength)*bSpline(xDiff, xLength)*bSpline(zDiff, zLength));
    vec.z(gradientBSpline(zDiff, zLength)*bSpline(xDiff, xLength)*bSpline(yDiff, yLength));
    return vec;

  }


  double
  MPMBox::gradientBSpline(double diff, double length)
  {
    double result = 0.0;
    if ((diff >= (-1.5)*length) && (diff <= (-0.5)*length))
    {
      result = (1/(length*length))*diff + 3/(2*length);
    }
    else if ((diff >= (-0.5)*length) && (diff <= 0.5*length)) 
    {
      result = - (2/(length*length))*diff;
    }
    else if ((diff >= 0.5*length) && (diff <= 1.5*length))
    {
      result = (1/(length*length))*diff - 3/(2*length);
    }
    else
    {
      result = 0.0;
    }
    
    return result;
  }

  void
  MPMBox::updatePointsVolume()
  {
    double determinent = 1.0;
//    std::vector<MaterialPointP> pointsArray = d_box->getMaterialPoints();
    std::vector<MaterialPointP> pointsArray = getPoints();
    int pointSize = pointsArray.size();
//    int elementSize = elementsArray.size();
//    Vector3D force(0.0, 0.0, 0.0);
    for (int ii = 0; ii < pointSize; ii++)
    {
      MaterialPointP cur_point = pointsArray[ii];
      determinent = cur_point->getDeformationGradientOld().Determinant();
      cur_point->setVolume(determinent*cur_point->getInitialVolume());
      determinent = 1.0;
//      if (cur_point->getID() == 50)
//        std::cout << "Volume= " << cur_point->getVolume() << std::endl;      
    }
    

  }

  Vector3D
  MPMBox::internalForce(NodeP node)
  {
    //std::vector<NodeP> nodesArray = d_box->getNodes();
//    std::vector<SymmetricMaterialTrilinearElementSP> elementsArray = d_box->getElements();

    //std::vector<NodeP> nodesArray = Boxx()->getNodes();
//    std::vector<SymmetricMaterialTrilinearElementSP> elementsArray = Boxx()->getElements();
    double density = 0.0;
    if (d_box_flag == true)
    {
      density = d_box->getDensity();
    }
    else
    { 
      density = d_geometry->getDensity();
    }
//    std::cout << density << std::endl;
    double determinent = 1.0;
//    std::vector<MaterialPointP> pointsArray = d_box->getMaterialPoints();
    std::vector<MaterialPointP> pointsArray = getPoints();
    int pointSize = pointsArray.size();
//    int elementSize = elementsArray.size();
    Vector3D force(0.0, 0.0, 0.0);
    for (int ii = 0; ii < pointSize; ii++)
    {
      MaterialPointP cur_point = pointsArray[ii];
//      determinent = cur_point->getDeformationGradientOld().Determinant();
//      cur_point->setVolume(determinent*cur_point->getVolume());
      force = force + (cur_point->getStressOld()*gradientInterpolate(node, cur_point))*
                      //(cur_point->getVolume());
                      (cur_point->getMass()/density);
//      cur_point->getStressOld().printMatrix(); 
   
    }
    return force*(-1);
  }



  double
  MPMBox::consistentMassComponent(NodeP node1, NodeP node2)
  {
//    std::vector<MaterialPointP> pointsArray = d_box->getMaterialPoints();
    std::vector<MaterialPointP> pointsArray = getPoints();
    int pointSize = pointsArray.size();
//    std::cout << std::endl << pointSize;
    double mass = 0.0;
    for (int ii = 0; ii < pointSize; ii++)
    {
      MaterialPointP cur_point = pointsArray[ii];
      mass += cur_point->getMass()*interpolate(node1, cur_point)*interpolate(node2, cur_point);
    }
   return mass;
  }



  double
  MPMBox::lumpedMassComponent(NodeP node)  // sigmaOverPoints(massOfPoint*interpolateOfPoint(node))
  { 
    
//    std::cout << "ID of input node of lumpedMassComponent function: " << node->getID() << std::endl;
//    std::vector<MaterialPointP> pointsArray = d_box->getMaterialPoints();
    std::vector<MaterialPointP> pointsArray = getPoints();
    int pointSize = pointsArray.size();
//    std::cout << "The number of points: " << pointSize << std::endl;
    double mass = 0.0;
    for (int ii = 0; ii < pointSize; ii++)
    {
      MaterialPointP cur_point = pointsArray[ii];

//      std::cout << "Mass= " << cur_point->getMass() << std::endl;
//      std::cout << "InitialVolume= " << cur_point->getInitialVolume()
//                << " Volume= " << cur_point->getVolume() << std::endl;
      mass += cur_point->getMass()*interpolate(node, cur_point);
//      if (node->getID() == 19)
//      {
//        std::cout << "PointID= " << cur_point->getID() << " Mass= " << cur_point->getMass() << 
//                     " Volume= " << cur_point->getInitialVolume() << "  interpolate= " << interpolate(node, cur_point) << std::endl;
//      }
    }

 //   std::cout << "Node's id: " << node->getID() << " Mass of the node: " << mass << std::endl;
    return mass;
  }



  Vector3D
  MPMBox::momentumComponent(NodeP node)
  {
//    std::vector<MaterialPointP> pointsArray = d_box->getMaterialPoints();
    std::vector<MaterialPointP> pointsArray = getPoints();
    int pointSize = pointsArray.size();
    Vector3D momentum(0.0);
    for (int ii = 0; ii < pointSize; ii++)
    {
      MaterialPointP cur_point = pointsArray[ii];
//      std::cout << "Point " << cur_point->getID() << " Velocity= " << cur_point->getVelOld() << std::endl;
      momentum = momentum + cur_point->getVelOld()*(interpolate(node, cur_point)*cur_point->getMass());
    }
   return momentum;
  }


  Vector3D
  MPMBox::newMomentumComponent(NodeP node)
  {
//    std::vector<MaterialPointP> pointsArray = d_box->getMaterialPoints();
    std::vector<MaterialPointP> pointsArray = getPoints();
    int pointSize = pointsArray.size();
    Vector3D momentum(0.0);
    for (int ii = 0; ii < pointSize; ii++)
    {
      MaterialPointP cur_point = pointsArray[ii];
      momentum = momentum + cur_point->getVelNew()*(interpolate(node, cur_point)*cur_point->getMass());
    }
   return momentum;
  }


  Vector3D
  MPMBox::bodyForceComponent(NodeP node)
  {
//    std::vector<MaterialPointP> pointsArray = d_box->getMaterialPoints();
    std::vector<MaterialPointP> pointsArray = getPoints();
    int pointSize = pointsArray.size();
    Vector3D force(0.0);
    Vector3D zero(0.0);
    for (int ii = 0; ii < pointSize; ii++)
    {
      MaterialPointP cur_point = pointsArray[ii]; 
      if (!(cur_point->getBodyForce() == zero))
      {
        std::cout << "ID= " << cur_point->getID() << " BodyForce= " << cur_point->getBodyForce() << std::endl;
      }
      force += cur_point->getBodyForce()*(interpolate(node, cur_point)*cur_point->getMass());
    }
//    std::cout << "NodeID= " << node->getID() << " Force= " << force << std::endl;
   return force;
  }


  Vector3D
  MPMBox::pointForceComponent(NodeP node)
  {
//    std::vector<MaterialPointP> pointsArray = d_box->getMaterialPoints();
    std::vector<MaterialPointP> pointsArray = getPoints();
    int pointSize = pointsArray.size();
    Vector3D force(0.0);
    Vector3D zero(0.0);
    for (int ii = 0; ii < pointSize; ii++)
    {
      MaterialPointP cur_point = pointsArray[ii];
//      if (!(cur_point->getPointForce() == zero))
//      {
//        std::cout << "ID= " << cur_point->getID() << " PointForce= " << cur_point->getPointForce() << std::endl;
//      }
//      std::cout << cur_point->getPointForce() << std::endl;
      force += cur_point->getPointForce()*interpolate(node, cur_point);
    }
   return force;
  }


  double
  MPMBox::matrixElement(NodeP node, MaterialPointP material, int t, int s)
  {
    double sum = 0.0;
//    Tensor tensor = d_box->getElements()[0]->tensor();
    Tensor tensor = Boxx()->getElements()[0]->tensor();
    Vector3D grad = gradientInterpolate(node, material);
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        sum = sum + tensor.get(t, i, j, s)*grad[i]*grad[j];
        sum = sum + tensor.get(t, i, s, j)*grad[i]*grad[j];
      } 
    }

//   std::cout << "(" << t << " " << s <<")  " << sum << std::endl;

   return sum;
  }


  Matrix3D
  MPMBox::forceIncrementComponent(NodeP node)
  {
    Matrix3D mat(0.0);
//    double density = d_box->getDensity();
//    std::vector<MaterialPointP> pointsArray = d_box->getMaterialPoints();

    double density = Boxx()->getDensity();
    std::vector<MaterialPointP> pointsArray = getPoints();

    int pointSize = pointsArray.size();
    for (int ii = 0; ii < pointSize; ii++)
    {
      MaterialPointP cur_point = pointsArray[ii];
      Matrix3D point_mat(0.0);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          point_mat.set(i, j, matrixElement(node, cur_point, i, j)); 
      mat += point_mat*(cur_point->getMass()/density);
                     // *(cur_point->getVolume());
 
    }
//   mat.printMatrix();
//   std::cout << std::endl;
   return mat*(-1); 
      

  }




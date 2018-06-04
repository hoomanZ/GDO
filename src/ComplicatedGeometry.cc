#include <ComplicatedGeometry.h>
#include <Node.h>
//#include <NodeP.h>
#include <SymmetricMaterialTrilinearElementSP.h>
#include <SymmetricMaterialTrilinearElement.h>
//#include <MaterialPointP.h>
#include <MaterialPoint.h>

#include <GeometryMath/Vector3D.h>
#include <GeometryMath/Matrix3D.h>
#include <GeometryMath/Point3D.h>

#include <array>
#include <vector>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>


using namespace FiniteElement;


ComplicatedGeometry::ComplicatedGeometry()
  :d_nodes_file_name("BoneDataNodes.txt"), d_elements_file_name("BoneDataElements.txt"),
   d_num_points_sai1(2), d_num_points_sai2(2), d_num_points_sai3(2),
   d_traction_force(0.0), d_boundary_displacement(0.0) 
{
  d_points_array.reserve(10000);
  createNodesArray();
  createElementsArray();
}


ComplicatedGeometry::ComplicatedGeometry(double density, double young, double poission,
                                         Vector3D tractionForce, Vector3D boundaryDisplacement)
  :d_nodes_file_name("BoneDataNodes.txt"), d_elements_file_name("BoneDataElements.txt"),
   d_num_points_sai1(2), d_num_points_sai2(2), d_num_points_sai3(2)
{
  Body::setProperties(density, young, poission);
  setTractionForce(tractionForce);
  setBoundaryDisplacement(boundaryDisplacement);
  d_points_array.reserve(10000);
  createNodesArray();
  createElementsArray();
}



ComplicatedGeometry::ComplicatedGeometry(std::string nodeFileName, std::string elementFileName)
  :d_nodes_file_name(nodeFileName), d_elements_file_name(elementFileName),
   d_num_points_sai1(2), d_num_points_sai2(2), d_num_points_sai3(2),
   d_traction_force(0.0), d_boundary_displacement(0.0)  
{
//  Body::setProperties(density, young, poission);
  d_points_array.reserve(10000);
  createNodesArray();
  createElementsArray();
}


ComplicatedGeometry::ComplicatedGeometry(std::string nodeFileName, std::string elementFileName,
                                         double density, double young, double poission,
                                         Vector3D tractionForce, Vector3D boundaryDisplacement)
  :d_nodes_file_name(nodeFileName), d_elements_file_name(elementFileName),
   d_num_points_sai1(2), d_num_points_sai2(2), d_num_points_sai3(2)  
{
  Body::setProperties(density, young, poission);
  setTractionForce(tractionForce);
  setBoundaryDisplacement(boundaryDisplacement);
  d_points_array.reserve(10000);
  createNodesArray();
  createElementsArray();
}


ComplicatedGeometry::ComplicatedGeometry(std::string nodeFileName, std::string elementFileName,
                                         std::string tractionForceNodesFile, std::string dirichletBoundaryNodesFile,
                                         double density, double young, double poission,
                                         Vector3D tractionForce, Vector3D boundaryDisplacement)
  :d_nodes_file_name(nodeFileName), d_elements_file_name(elementFileName),
   d_num_points_sai1(2), d_num_points_sai2(2), d_num_points_sai3(2)
{
  Body::setProperties(density, young, poission);
//  std::cout << "setProperties" << std::endl;
  setTractionForce(tractionForce);
  setBoundaryDisplacement(boundaryDisplacement);
  d_points_array.reserve(10000);
  createNodesArray();
  createElementsArray();
//  std::cout << "Nodes and Elements are created" << std::endl;
  boundaryConditionNodesID(tractionForceNodesFile, d_traction_force_nodes); 
  boundaryConditionNodesID(dirichletBoundaryNodesFile, d_dirichlet_boundary_nodes);
  createTractionForce();
//  std::cout << "traction forces" << std::endl;
//  Body::makeGlobalMatrices();

}


ComplicatedGeometry::ComplicatedGeometry(std::string nodeFileName, std::string elementFileName,
                                         double density, double young, double poission,
                                         Vector3D tractionForce, Vector3D boundaryDisplacement,
                                         int numSai1, int numSai2, int numSai3)
  :d_nodes_file_name(nodeFileName), d_elements_file_name(elementFileName),
   d_num_points_sai1(numSai1), d_num_points_sai2(numSai2), d_num_points_sai3(numSai3)  
{
  Body::setProperties(density, young, poission);
  setTractionForce(tractionForce);
  setBoundaryDisplacement(boundaryDisplacement);
  d_points_array.reserve(10000);
  createNodesArray();
  createElementsArray();
  createPointsArray();
  createRefinedElements();
  calculateNodesVolume();
  pointsFile();
}



ComplicatedGeometry::ComplicatedGeometry(std::string nodeFileName, std::string elementFileName,
                                         std::string tractionForceNodesFile, std::string dirichletBoundaryNodesFile,
                                         double density, double young, double poission,
                                         Vector3D tractionForce, Vector3D boundaryDisplacement,
                                         int numSai1, int numSai2, int numSai3)
  :d_nodes_file_name(nodeFileName), d_elements_file_name(elementFileName),
   d_num_points_sai1(numSai1), d_num_points_sai2(numSai2), d_num_points_sai3(numSai3)  
{
  Body::setProperties(density, young, poission);
  setTractionForce(tractionForce);
  setBoundaryDisplacement(boundaryDisplacement);
  d_points_array.reserve(10000);
  createNodesArray();
  createElementsArray();
  createPointsArray();
  createRefinedElements();
  calculateNodesVolume();
//  std::string modifiedTractionNodesFile = modifyNodesID(Body::getNodes(), tractionForceNodesFile);
  
  boundaryConditionNodesID(tractionForceNodesFile, d_point_force_points);
  modifyPointsID(d_point_force_points);
  boundaryConditionNodesID(tractionForceNodesFile, d_traction_force_nodes);
  modifyNodesID();//tractionForceNodesFile); 
  boundaryConditionNodesID(dirichletBoundaryNodesFile, d_dirichlet_boundary_nodes);
  modifyPointsID(d_dirichlet_boundary_nodes);  
//  modifyNodesID();
//  createTractionForce();
  createRefinedTractionForce();
  calculateNodesArea();
//  std::cout << "Successfull constructor" << std::endl;  
//  pointsFile();
}


ComplicatedGeometry::ComplicatedGeometry(std::string nodeFileName, std::string elementFileName,
                                         std::string tractionForceNodesFile, std::string dirichletBoundaryNodesFile,
                                         double density, double young, double poission,
                                         Vector3D tractionForce, Vector3D boundaryDisplacement,
                                         int numSai1, int numSai2, int numSai3,
                                         double xNumElem, double yNumElem, double zNumElem)
  :d_nodes_file_name(nodeFileName), d_elements_file_name(elementFileName),
   d_num_points_sai1(numSai1), d_num_points_sai2(numSai2), d_num_points_sai3(numSai3), d_box(new Box())  
{
  d_box->setNumElementsX(xNumElem);
  d_box->setNumElementsY(yNumElem);
  d_box->setNumElementsZ(zNumElem);

  Body::setProperties(density, young, poission);
//  setTractionForce(tractionForce);
  setPointForce(tractionForce);
  setBoundaryDisplacement(boundaryDisplacement);
  d_points_array.reserve(10000);
  createNodesArray();
  createElementsArray();
  createPointsArray();
  createRefinedElements();
  calculateNodesVolume();
//  std::string modifiedTractionNodesFile = modifyNodesID(Body::getNodes(), tractionForceNodesFile);
  
  boundaryConditionNodesID(tractionForceNodesFile, d_point_force_points);
  modifyPointsID(d_point_force_points);
  boundaryConditionNodesID(tractionForceNodesFile, d_traction_force_nodes);
  modifyNodesID();//tractionForceNodesFile); 
  boundaryConditionNodesID(dirichletBoundaryNodesFile, d_dirichlet_boundary_nodes);
  modifyPointsID(d_dirichlet_boundary_nodes);  
//  modifyNodesID();
//  createTractionForce();
////  createRefinedTractionForce();
  createPointForce(); 
  calculateNodesArea();
//  std::cout << "Successfull constructor" << std::endl;  
//  pointsFile();
}




ComplicatedGeometry::~ComplicatedGeometry()
{
} 


void
ComplicatedGeometry::createNodesArray()
{
  int nodeID = 1;
  std::vector <NodeP> nodesArray;
  std::ifstream infile(d_nodes_file_name);
//  std::ifstream infile("BoneData.txt");

 
  if (infile.is_open())
  {
    for (std::string line; std::getline(infile, line); )
    {  
      double x, y, z, s;
      if (!(std::istringstream(line) >> x >> y >> z >> s)) continue; 
      else 
      {
//        std::istringstream(line) >> x >> y >> z >> s;
        NodeP node(new Node(nodeID, x, y, z, false));
        node->densityNode(Body::getDensity());
//        std::cout << "NodeID= " << node->getID() << " Velocity= " << node->velocity() << std::endl;
        nodesArray.emplace_back(node);    
        nodeID += 1;
      }

    }
  }
  else
  {
    std::cout << "Unable to open the file" <<std::endl;
  }

  infile.close();
  Body::setNodes(nodesArray);

//  std::cout << d_nodes_file_name << std::endl;
//  std::cout << Body::getNodes().size() << std::endl; 

//  std::vector <NodeP> nodes = Body::getNodes();
//  int size = nodes.size();
//  for (int i = 0; i < size; i++)
//  {
//    NodeP cur_node = nodes[i];
//    std::cout << " Node ID= " << cur_node->getID() << " (" << cur_node->x() << ", " << cur_node->y() << ", " << cur_node->z() << ")" << std::endl;
//  }



////////////////////////////////////// for element #45
//  Vector3D v5(-52.7934, -23.82357, -61.67243);
//  Vector3D v6(-53.74731, -30.12011, -62.91138);
//  Vector3D v7(-50.0757, -31.48918, -64.13622);
//  Vector3D v8(-49.75115, -24.14222, -62.72066);

//  Vector3D f5(105.093, -109.347, 647.173);
//  Vector3D f6(98.8691, -102.871, 608.849);
//  Vector3D f7(106.795, -111.118, 657.656);
//  Vector3D f8(113.022, -117.597, 696.002);

//  double a = ((v5 - v6).cross(v7 - v6)).length();
//  double b = ((v5 - v8).cross(v7 - v8)).length();

//  std::cout << "S= " << 0.5*(a+b)*100 << std::endl;
//  std::cout << "s= " << f5.length() + f6.length() 
//          + f7.length() + f8.length() << std::endl;
///////////////////////////////////////////////////////

}


void
ComplicatedGeometry::createElementsArray()
{
  double density = Body::getDensity();
  double poission = Body::getPoission();
  double young = Body::getYoung();
  int elementID = 1;
  std::vector <SymmetricMaterialTrilinearElementSP>  elementsArray;   
  std::ifstream infile(d_elements_file_name);  
  for (std::string line; std::getline(infile, line); )
  {  
    int six, seven;
    int eight, five;
    int one, two;
    int three, four; 
    if (!(std::istringstream(line) >> six >> seven >> eight >> five >> one >> two >> three >> four)) continue; 
    SymmetricMaterialTrilinearElementSP element(new SymmetricMaterialTrilinearElement
                                                      (elementID, density, young, poission, Body::getNodes(), 
                                                      one, two, three, four, five, six, seven, eight));
    elementsArray.emplace_back(element);
    elementID++;

  }
  Body::setElements(elementsArray);

//  std::cout << d_elements_file_name << std::endl;
//  std::cout << Body::getElements().size() << std::endl; 

//  std::vector <SymmetricMaterialTrilinearElementSP> elements = Body::getElements();
//  int size = elements.size();
//  for (int i = 0; i < size; i++)
//  {
//    SymmetricMaterialTrilinearElementSP cur_element = elements[i];
//    std::cout << " Element ID= " << cur_element->id() << " (";
//    TrilinearVolumeElementSP element = cur_element->triElement();
//    std::array <NodeP, 8> nodes = element->getArray();
//    for (int j = 0; j < 8; j++)
//    {
//      NodeP node = nodes[j]; 
//      if (j < 7)
//      {   
//        std::cout << node->getID() << ", ";
//      }
//      else
//      {
//        std::cout << node->getID() << ")" << std::endl;
//      }
//    }
//    std::cout << ")" << std::endl;
//  }
}


void
ComplicatedGeometry::createPointsArray()
{
  std::vector <SymmetricMaterialTrilinearElementSP> elements = Body::getElements();
  int index = 1;
  int nodeIndex = 1;
  MaterialPointPArray pointArray;
  int size = elements.size();
  for (int ii = 0; ii < size; ii++)
  {
    SymmetricMaterialTrilinearElementSP cur_element = elements[ii];
    TrilinearVolumeElementSP element = cur_element->triElement();
    std::vector <std::vector <double> >  nodesInformation;
    nodesInformation.clear();
//    std::array <NodeP, 8> nodes = element->getArray();

    double nx = d_num_points_sai1 - 1; 
    double ny = d_num_points_sai2 - 1;
    double nz = d_num_points_sai3 - 1;
 
    double dsai1 = 2/nx; 
    double dsai2 = 2/ny;
    double dsai3 = 2/nz;
    for (int k = 0; k < nz+1; k++)
     {
       for (int j = 0; j < ny+1; j++)
       {
         for (int i = 0; i < nx+1; i++)
         {
           Vector3D p(-1 + i*dsai1, -1 + j*dsai2, -1 + k*dsai3);
           double x = 0;
           double y = 0;
           double z = 0;
           std::array <NodeP, 8> nodes = element->getArray();
           for (int s = 0; s < 8; s++)
           {
             x += nodes[s]->x()*element->symTestFunction(s+1, p[0], p[1], p[2]);
             y += nodes[s]->y()*element->symTestFunction(s+1, p[0], p[1], p[2]);
             z += nodes[s]->z()*element->symTestFunction(s+1, p[0], p[1], p[2]);
           }
           Vector3D position(x, y, z);
           std::vector <double>  nodesInf;
           nodesInf.clear();
           std::vector <NodeP>  refinedNodes = getRefinedNodesArray();
           int curIndex = 1;
           if (nodeRepeated(position, refinedNodes, curIndex))
           {
             nodesInf.emplace_back(curIndex);
             nodesInf.emplace_back(-1 + i*dsai1);
             nodesInf.emplace_back(-1 + j*dsai2);
             nodesInf.emplace_back(-1 + k*dsai3);
             nodesInformation.emplace_back(nodesInf);
             continue;              
           }
           else
           {
             NodeP node(new Node(nodeIndex, x, y, z, false));
             refinedNodes.emplace_back(node);
             setRefinedNodesArray(refinedNodes);

             nodesInf.emplace_back(nodeIndex);
             nodesInf.emplace_back(-1 + i*dsai1);
             nodesInf.emplace_back(-1 + j*dsai2);
             nodesInf.emplace_back(-1 + k*dsai3);
             nodesInformation.emplace_back(nodesInf);

             nodeIndex++;             
           }
           MaterialPointPArray pointArray = getMaterialPointsArray();
           if (pointRepeated(position, pointArray)) continue;
           else
           {
             MaterialPointP point(new MaterialPoint(index, x, y, z)); 
             index++;
             pointArray.emplace_back(point);
             setMaterialPointsArray(pointArray);
           }          
             
//           NodeP node(new Node(nodeID(i, j, k),
//                      d_min_point.x()+i*dx,
//                      d_min_point.y()+j*dy,
//                      d_min_point.z()+k*dz, false));
//           surfaceOnX = ((node->x() == d_min_point.x()) || (node->x() == d_max_point.x()));
//           surfaceOnY = ((node->y() == d_min_point.y()) || (node->y() == d_max_point.y()));
//           surfaceOnZ = ((node->z() == d_min_point.z()) || (node->z() == d_max_point.z()));
//           node->onSurface(surfaceOnX || surfaceOnY || surfaceOnZ);
//           node->densityNode(Body::getDensity());
//           nodesArray.emplace_back(node);
         }
       }
     }
     cur_element->setRefinedNodes(nodesInformation);


  }


//  std::vector <NodeP>  refinedNodes = getRefinedNodesArray();
//  int siz = refinedNodes.size();
//  for (int kk = 0; kk < siz; kk++)
//  {
//    NodeP cur_node = refinedNodes[kk];
//    std::cout << "Node ID " << cur_node->getID() << " (" << cur_node->x() << ", " << cur_node->y() 
//                                                   << ", " << cur_node->z() << ")" << std::endl;
//  }


//  MaterialPointPArray array = getMaterialPointsArray();
//  int siz = array.size();
//  for (int kk = 0; kk < siz; kk++)
//  {
//    MaterialPointP cur_point = array[kk];
//    std::cout << "Point ID " << cur_point->getID() << " (" << cur_point->getPosOld().x() 
//              << ", " << cur_point->getPosOld().y() 
//              << ", " << cur_point->getPosOld().z() << ")" << std::endl;
//  }

//  std::cout << getMaterialPointsArray().size() << std::endl;
//  std::cout << getRefinedNodesArray().size() << std::endl;
}



bool
ComplicatedGeometry::pointRepeated(const Vector3D& position, const MaterialPointPArray& pointArray)
{
  int size = pointArray.size();
  double Tol = 0.0000001;
  for (int k = 0; k < size; k++)
  {
    MaterialPointP cur_point = pointArray[k];
    Point3D p = cur_point->getPosOld();
    Vector3D pos(p[0], p[1], p[2]);
    double distance = (position.x() - pos.x())*(position.x() - pos.x())+
                      (position.y() - pos.y())*(position.y() - pos.y())+
                      (position.z() - pos.z())*(position.z() - pos.z());
    distance = std::sqrt(distance);
    if (distance < Tol) return true;

  }
  return false;

}


bool
ComplicatedGeometry::nodeRepeated(const Vector3D& position, const std::vector <NodeP>& refinedNodes, int& index)
{
  int size = refinedNodes.size();
  double Tol = 0.0000001;
  for (int k = 0; k < size; k++)
  {
    NodeP cur_node = refinedNodes[k];
    Point3D p = cur_node->position();
    Vector3D pos(p[0], p[1], p[2]);
    double distance = (position.x() - pos.x())*(position.x() - pos.x())+
                      (position.y() - pos.y())*(position.y() - pos.y())+
                      (position.z() - pos.z())*(position.z() - pos.z());
    distance = std::sqrt(distance);
    if (distance < Tol)
    {
      index = cur_node->getID();
      return true;
    }

  }

  return false;

}


int
ComplicatedGeometry::findNodeIndex(const SymmetricMaterialTrilinearElementSP& element,
                                  const double& sai1, const double& sai2, const double& sai3)
{
  std::vector <std::vector <double> >  nodesInformation = element->getRefinedNodes();
  int size = nodesInformation.size();
  double Tol = 0.0000001;
  for (int i = 0; i < size; i++)
  { 
    std::vector <double> nodeInf = nodesInformation[i];
    double distance = (nodeInf[1] - sai1)*(nodeInf[1] - sai1)+
                      (nodeInf[2] - sai2)*(nodeInf[2] - sai2)+
                      (nodeInf[3] - sai3)*(nodeInf[3] - sai3);
    distance = std::sqrt(distance);
    if (distance < Tol)
    {
      return nodeInf[0];
    }
    
  }
  std::cout << "There is no node in the element with this (" << sai1 << ", " << sai2 << ", " << sai3 << ")" << std::endl; 
  return 1;
}


void
ComplicatedGeometry::createRefinedElements()
{
  std::vector <SymmetricMaterialTrilinearElementSP> elements = Body::getElements();
  int size = elements.size();

  double density = Body::getDensity();
  double poission = Body::getPoission();
  double young = Body::getYoung();
  int elementID = 1;

  std::vector <SymmetricMaterialTrilinearElementSP>  refinedElementsArray;
  refinedElementsArray.clear();

  for (int ii = 0; ii < size; ii++)
  {   
    SymmetricMaterialTrilinearElementSP cur_element = elements[ii];
//    TrilinearVolumeElementSP element = cur_element->triElement();
    
    double nx = d_num_points_sai1 - 1; 
    double ny = d_num_points_sai2 - 1;
    double nz = d_num_points_sai3 - 1;
 
    double dsai1 = 2/nx; 
    double dsai2 = 2/ny;
    double dsai3 = 2/nz;

    int n1, n2, n3, n4;
    int n5, n6, n7, n8;
    for (int k = 0; k < nz; k++)
    {
      for (int j = 0; j < ny; j++)
      {
        for (int i = 0; i < nx; i++)
        {
//          n1 = nodeID(i, j, k);
          n1 = findNodeIndex(cur_element, -1 + i*dsai1, -1 + j*dsai2, -1 + k*dsai3); 
//          n2 = nodeID(i+1, j, k);
          n2 = findNodeIndex(cur_element, -1 + (i+1)*dsai1, -1 + j*dsai2, -1 + k*dsai3);           
//           n3 = nodeID(i+1, j+1, k);//nodeID(i, j+1, k);
          n3 = findNodeIndex(cur_element, -1 + (i+1)*dsai1, -1 + (j+1)*dsai2, -1 + k*dsai3); 
//           n4 = nodeID(i, j+1, k);//nodeID(i+1, j+1, k);
          n4 = findNodeIndex(cur_element, -1 + i*dsai1, -1 + (j+1)*dsai2, -1 + k*dsai3); 
//           n5 = nodeID(i, j+1, k+1);//nodeID(i, j, k+1);
          n5 = findNodeIndex(cur_element, -1 + i*dsai1, -1 + (j+1)*dsai2, -1 + (k+1)*dsai3); 
//           n6 = nodeID(i, j, k+1);//nodeID(i+1, j, k+1);
          n6 = findNodeIndex(cur_element, -1 + i*dsai1, -1 + j*dsai2, -1 + (k+1)*dsai3); 
//           n7 = nodeID(i+1, j, k+1);//nodeID(i, j+1, k+1);
          n7 = findNodeIndex(cur_element, -1 + (i+1)*dsai1, -1 + j*dsai2, -1 + (k+1)*dsai3); 
//           n8 = nodeID(i+1, j+1, k+1);
          n8 = findNodeIndex(cur_element, -1 + (i+1)*dsai1, -1 + (j+1)*dsai2, -1 + (k+1)*dsai3);
 
//          if (cur_element->id() < 20)
//            std::cout << n1 << " " << n2 << " " << n3 << " " << n4 << " "
//                      << n5 << " " << n6 << " " << n7 << " " << n8 << std::endl;

          SymmetricMaterialTrilinearElementSP element(new SymmetricMaterialTrilinearElement
                                                          (elementID, density, young, poission, getRefinedNodesArray(),          
                                                           n1, n2, n3, n4, n5, n6, n7, n8));
          refinedElementsArray.emplace_back(element);
          elementID++;
    
        }
      }
    }  

  }

  setRefinedElementsArray(refinedElementsArray);

/*  std::vector <SymmetricMaterialTrilinearElementSP> elementsArray = getRefinedElementsArray();
  int sizeArray = elementsArray.size();
  for (int i = 0; i < sizeArray; i++)
  {
    SymmetricMaterialTrilinearElementSP cur_elementArray = elementsArray[i];
    std::cout << " Element ID= " << cur_elementArray->id() << " (";
    TrilinearVolumeElementSP elementArray = cur_elementArray->triElement();
    std::array <NodeP, 8> nodes = elementArray->getArray();
    for (int j = 0; j < 8; j++)
    {
      NodeP node = nodes[j]; 
      if (j < 7)
      {   
        std::cout << node->getID() << ", ";
      }
      else
      {
        std::cout << node->getID() << ")" << std::endl;
      }
    }
//    std::cout << ")" << std::endl;
  }*/

}




void 
ComplicatedGeometry::calculateNodesVolumeOfElement(SymmetricMaterialTrilinearElementSP& element)
{
  TrilinearVolumeElementSP elem = element->triElement();
  std::array <NodeP, 8> nodes = elem->getArray();
  NodeP node1 = nodes[0]; 
  std::vector <double>  node1VolumeArray = node1->getVolumeNodes(); 
  Point3D pos1 = node1->position();
  Vector3D node1Pos(pos1.x(), pos1.y(), pos1.z());
  double volume1 = nodeVolume(node1Pos, findGlobalPosition(element, 0, -1, -1),
                                        findGlobalPosition(element, 0, 0, -1),  
                                        findGlobalPosition(element, -1, 0, -1),
                                        findGlobalPosition(element, -1, 0, 0),
                                        findGlobalPosition(element, -1, -1, 0),
                                        findGlobalPosition(element, 0, -1, 0),
                                        findGlobalPosition(element, 0, 0, 0));
  node1VolumeArray.emplace_back(volume1);
  node1->setVolumeNodes(node1VolumeArray);
   
  NodeP node2 = nodes[1]; 
  std::vector <double>  node2VolumeArray = node2->getVolumeNodes(); 
  Point3D pos2 = node2->position();
  Vector3D node2Pos(pos2.x(), pos2.y(), pos2.z());
  double volume2 = nodeVolume(findGlobalPosition(element, 0, -1, -1),
                              node2Pos,
                              findGlobalPosition(element, 1, 0, -1),  
                              findGlobalPosition(element, 0, 0, -1),
                              findGlobalPosition(element, 0, 0, 0),
                              findGlobalPosition(element, 0, -1, 0),
                              findGlobalPosition(element, 1, -1, 0),
                              findGlobalPosition(element, 1, 0, 0));
  node2VolumeArray.emplace_back(volume2);
  node2->setVolumeNodes(node2VolumeArray);

  NodeP node3 = nodes[2]; 
  std::vector <double>  node3VolumeArray = node3->getVolumeNodes(); 
  Point3D pos3 = node3->position();
  Vector3D node3Pos(pos3.x(), pos3.y(), pos3.z());
  double volume3 = nodeVolume(findGlobalPosition(element, 0, 0, -1),                                     
                              findGlobalPosition(element, 1, 0, -1),
                              node3Pos,  
                              findGlobalPosition(element, 0, 1, -1),
                              findGlobalPosition(element, 0, 1, 0),
                              findGlobalPosition(element, 0, 0, 0),
                              findGlobalPosition(element, 1, 0, 0),
                              findGlobalPosition(element, 1, 1, 0));
  node3VolumeArray.emplace_back(volume3);
  node3->setVolumeNodes(node3VolumeArray);

  NodeP node4 = nodes[3]; 
  std::vector <double>  node4VolumeArray = node4->getVolumeNodes(); 
  Point3D pos4 = node4->position();
  Vector3D node4Pos(pos4.x(), pos4.y(), pos4.z());
  double volume4 = nodeVolume(findGlobalPosition(element, -1, 0, -1),                                     
                              findGlobalPosition(element, 0, 0, -1),                                         
                              findGlobalPosition(element, 0, 1, -1),
                              node4Pos,
                              findGlobalPosition(element, -1, 1, 0),
                              findGlobalPosition(element, -1, 0, 0),
                              findGlobalPosition(element, 0, 0, 0),
                              findGlobalPosition(element, 0, 1, 0));
  node4VolumeArray.emplace_back(volume4);
  node4->setVolumeNodes(node4VolumeArray);

  NodeP node5 = nodes[4]; 
  std::vector <double>  node5VolumeArray = node5->getVolumeNodes(); 
  Point3D pos5 = node5->position();
  Vector3D node5Pos(pos5.x(), pos5.y(), pos5.z());
  double volume5 = nodeVolume(findGlobalPosition(element, -1, 0, 0),                                     
                              findGlobalPosition(element, 0, 0, 0),                                         
                              findGlobalPosition(element, 0, 1, 0),                                       
                              findGlobalPosition(element, -1, 1, 0),
                              node5Pos,
                              findGlobalPosition(element, -1, 0, 1),
                              findGlobalPosition(element, 0, 0, 1),
                              findGlobalPosition(element, 0, 1, 1));
  node5VolumeArray.emplace_back(volume5);
  node5->setVolumeNodes(node5VolumeArray);

  NodeP node6 = nodes[5]; 
  std::vector <double>  node6VolumeArray = node6->getVolumeNodes(); 
  Point3D pos6 = node6->position();
  Vector3D node6Pos(pos6.x(), pos6.y(), pos6.z());
  double volume6 = nodeVolume(findGlobalPosition(element, -1, -1, 0),                                     
                              findGlobalPosition(element, 0, -1, 0),                                         
                              findGlobalPosition(element, 0, 0, 0),                                       
                              findGlobalPosition(element, -1, 0, 0),                                       
                              findGlobalPosition(element, -1, 0, 1),
                              node6Pos,
                              findGlobalPosition(element, 0, -1, 1),
                              findGlobalPosition(element, 0, 0, 1));
  node6VolumeArray.emplace_back(volume6);
  node6->setVolumeNodes(node6VolumeArray);

  NodeP node7 = nodes[6]; 
  std::vector <double>  node7VolumeArray = node7->getVolumeNodes(); 
  Point3D pos7 = node7->position();
  Vector3D node7Pos(pos7.x(), pos7.y(), pos7.z());
  double volume7 = nodeVolume(findGlobalPosition(element, 0, -1, 0),                                     
                              findGlobalPosition(element, 1, -1, 0),                                         
                              findGlobalPosition(element, 1, 0, 0),                                       
                              findGlobalPosition(element, 0, 0, 0),                                       
                              findGlobalPosition(element, 0, 0, 1),                                       
                              findGlobalPosition(element, 0, -1, 1),
                              node7Pos,
                              findGlobalPosition(element, 1, 0, 1));

//  std::cout << findGlobalPosition(element, 0, -1, 0).x() << "  " <<findGlobalPosition(element, 0, -1, 0).y() << "  " << findGlobalPosition(element, 0, -1, 0).z() << std::endl;
//  std::cout << node7Pos.x() << " " << node7Pos.y() << " " << node7Pos.z() << std::endl;

  node7VolumeArray.emplace_back(volume7);
  node7->setVolumeNodes(node7VolumeArray);

  NodeP node8 = nodes[7]; 
  std::vector <double>  node8VolumeArray = node8->getVolumeNodes(); 
  Point3D pos8 = node8->position();
  Vector3D node8Pos(pos8.x(), pos8.y(), pos8.z());
  double volume8 = nodeVolume(findGlobalPosition(element, 0, 0, 0),                                     
                              findGlobalPosition(element, 1, 0, 0),                                         
                              findGlobalPosition(element, 1, 1, 0),                                       
                              findGlobalPosition(element, 0, 1, 0),                                       
                              findGlobalPosition(element, 0, 1, 1),                                       
                              findGlobalPosition(element, 0, 0, 1),                                      
                              findGlobalPosition(element, 1, 0, 1),
                              node8Pos);
//  std::cout << volume1 << " " << volume2 << " " << volume3 << " " << volume4 << " "
//            << volume5 << " " << volume6 << " " << volume7 << " " << volume8 << std::endl;
  node8VolumeArray.emplace_back(volume8);
  node8->setVolumeNodes(node8VolumeArray);

  elem->setArray(node1, node2, node3, node4,
                 node5, node6, node7, node8);
  element->triElement(elem);

}





Vector3D
ComplicatedGeometry::findGlobalPosition(const SymmetricMaterialTrilinearElementSP& element,
                                        const double& sai1, const double& sai2, const double& sai3)
{
  TrilinearVolumeElementSP trilinearElement = element->triElement();
  std::array <NodeP, 8> nodes = trilinearElement->getArray();
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  for (int s = 0; s < 8; s++)
  {
    x += nodes[s]->x()*trilinearElement->symTestFunction(s+1, sai1, sai2, sai3);
    y += nodes[s]->y()*trilinearElement->symTestFunction(s+1, sai1, sai2, sai3);
    z += nodes[s]->z()*trilinearElement->symTestFunction(s+1, sai1, sai2, sai3);
  }
  Vector3D v(x, y, z);
  return v;

}


double
ComplicatedGeometry::nodeVolume(Vector3D v1, Vector3D v2, Vector3D v3, Vector3D v4,
                                Vector3D v5, Vector3D v6, Vector3D v7, Vector3D v8)
{
  double volume = 0.0;

  Vector3D A1 = v3 - v6;
  Vector3D B1 = v7 - v6;
  Vector3D C1 = v2 - v8;
  Matrix3D M1(A1.x(), B1.x(), C1.x(),
              A1.y(), B1.y(), C1.y(),
              A1.z(), B1.z(), C1.z());

  Vector3D A2 = v3 - v6;
  Vector3D B2 = v5 - v6;
  Vector3D C2 = v8 - v4;
  Matrix3D M2(A2.x(), B2.x(), C2.x(),
              A2.y(), B2.y(), C2.y(),
              A2.z(), B2.z(), C2.z());

  Vector3D A3 = v3 - v6;
  Vector3D B3 = v1 - v6;
  Vector3D C3 = v4 - v2;
  Matrix3D M3(A3.x(), B3.x(), C3.x(),
              A3.y(), B3.y(), C3.y(),
              A3.z(), B3.z(), C3.z());

//  std::cout << M1.Determinant() << " " << M2.Determinant() << " " << M3.Determinant() << std::endl;

  volume = M1.Determinant() + M2.Determinant() + M3.Determinant();
  volume /= 6;

//  std::cout << volume << std::endl;
 
  return volume;

}



/*void
ComplicatedGeometry::calculateNodesVolume()
{
  std::vector <SymmetricMaterialTrilinearElementSP> elementsArray = getRefinedElementsArray();
  int size = elementsArray.size();
  for (int i = 0; i < size; i++)
  {
    SymmetricMaterialTrilinearElementSP cur_element = elementsArray[i];
    calculateNodesVolumeOfElement(cur_element);
  }

  std::vector <NodeP> nodes;
  nodes.clear();
  std::vector <NodeP>  nodesArray = getRefinedNodesArray();
  int nodeSize = nodesArray.size();  
  for (int j = 0; j < nodeSize; j++)
  {
    NodeP cur_node = nodesArray[j];
    std::vector <double> volumesArray = cur_node->getVolumeNodes();
    double volume = 0.0;
    int volumeSize = volumesArray.size();
    for (int k = 0; k < volumeSize; k++)
    {
      volume = volume + volumesArray[k];
    }
      
   cur_node->volume(volume);
   nodes.emplace_back(cur_node); 
  }
  
  setRefinedNodesArray(nodes);

 
  double totalVolume = 0.0;
  std::vector <NodeP>  nodess = getRefinedNodesArray(); 
  int siz = nodess.size();
  for (int kk = 0; kk < siz; kk++)
  {
    std::cout << "Node ID= " << nodess[kk]->getID() << " volume= " << nodess[kk]->volume() << std::endl;
    totalVolume += nodess[kk]->volume();
  }

  std::cout << "Total Volume= " << totalVolume << std::endl;

}*/


void
ComplicatedGeometry::calculateNodesVolume()
{
  std::vector <SymmetricMaterialTrilinearElementSP> elementsArray = getRefinedElementsArray();
  int size = elementsArray.size();
  for (int i = 0; i < size; i++)
  {
    SymmetricMaterialTrilinearElementSP cur_element = elementsArray[i];
    calculateNodesVolumeOfElement(cur_element);
  }

  std::vector <NodeP> nodes;
  MaterialPointPArray points;  
  nodes.clear();
  std::vector <NodeP>  nodesArray = getRefinedNodesArray();
  MaterialPointPArray pointsArray = getMaterialPointsArray();
  int nodeSize = nodesArray.size();  
  for (int j = 0; j < nodeSize; j++)
  {
    NodeP cur_node = nodesArray[j];
    MaterialPointP cur_point = pointsArray[j];
    std::vector <double> volumesArray = cur_node->getVolumeNodes();
    double volume = 0.0;
    int volumeSize = volumesArray.size();
    for (int k = 0; k < volumeSize; k++)
    {
      volume = volume + volumesArray[k];
    }
      
   cur_node->volume(volume);
   cur_point->setInitialVolume(volume);
   cur_point->setMass(cur_point->getInitialVolume()*getDensity());
   nodes.emplace_back(cur_node);
   points.emplace_back(cur_point); 
  }
  
  setRefinedNodesArray(nodes);
  setMaterialPointsArray(points);

 
//  double totalVolume = 0.0;
//  std::vector <NodeP>  nodess = getRefinedNodesArray(); 
//  int siz = nodess.size();
//  for (int kk = 0; kk < siz; kk++)
//  {
//    std::cout << "Node ID= " << nodess[kk]->getID() << " volume= " << nodess[kk]->volume() << std::endl;
//    totalVolume += nodess[kk]->volume();
//  }

//  std::cout << "Total Volume= " << totalVolume << std::endl;




//  double totalVolume = 0.0;
//  MaterialPointPArray  pointss = getMaterialPointsArray(); 
//  int siz = pointss.size();
//  for (int kk = 0; kk < siz; kk++)
//  {
//    std::cout << "Point ID= " << pointss[kk]->getID() << " volume= " << pointss[kk]->getInitialVolume() << std::endl;
//    totalVolume += pointss[kk]->getInitialVolume();
//  }

//  std::cout << "Total Volume= " << totalVolume << std::endl;

}



  void
  ComplicatedGeometry::nodesFile(std::string fileName)
  {
//    std::string pointFolderName = "pointsFolder";
//    int status2;
////    int ret2;
//    if (iter_num == 0)
//    {

//      makeOutputFolder(status2, ret2, pointFolderName);
//    }
//    else
//    {
//      std::string direcLocation = currentFolder + "/" + pointFolderName;
//      ret2 = chdir(direcLocation.c_str());
//    }
//    std::string fileName = "point.exdata";
    std::ofstream myfile(fileName);
    myfile << " Group name: node" // << "_" << std::to_string(iter_num) 
                                  << std::endl;
    myfile << " #Fields=1" << std::endl;
    myfile << " 1) coordinates, coordinate, rectangular cartesian, #Components=3" << std::endl;
    myfile << "   x.  Value index= 1, #Derivatives= 0" << std::endl;
    myfile << "   y.  Value index= 2, #Derivatives= 0" << std::endl;
    myfile << "   z.  Value index= 3, #Derivatives= 0" << std::endl;
//    BoxSP box = d_mpm_box->getBox();
//    std::vector<MaterialPointP> pointsArray = getMaterialPointsArray();
    std::vector<NodeP> nodesArray = Body::getNodes();
    int nodeSize = nodesArray.size();
    for (int k = 0; k < nodeSize; k++)
    {
      NodeP cur_node = nodesArray[k];
      myfile << " Node:         " << cur_node->getID() << std::endl;
      myfile << "  " << nodesArray[k]->position().x() << "  "
             << nodesArray[k]->position().y() << "  "
             << nodesArray[k]->position().z()  
             << std::endl;
    }
    
    myfile.close();
//    ret2 = chdir(currentFolder.c_str());

  }


  void
  ComplicatedGeometry::pointsFile()
  {
//    std::string pointFolderName = "pointsFolder";
//    int status2;
////    int ret2;
//    if (iter_num == 0)
//    {

//      makeOutputFolder(status2, ret2, pointFolderName);
//    }
//    else
//    {
//      std::string direcLocation = currentFolder + "/" + pointFolderName;
//      ret2 = chdir(direcLocation.c_str());
//    }
    std::string fileName = "point.exdata";
    std::ofstream myfile(fileName);
    myfile << " Group name: point" // << "_" << std::to_string(iter_num) 
                                  << std::endl;
    myfile << " #Fields=1" << std::endl;
    myfile << " 1) coordinates, coordinate, rectangular cartesian, #Components=3" << std::endl;
    myfile << "   x.  Value index= 1, #Derivatives= 0" << std::endl;
    myfile << "   y.  Value index= 2, #Derivatives= 0" << std::endl;
    myfile << "   z.  Value index= 3, #Derivatives= 0" << std::endl;
//    BoxSP box = d_mpm_box->getBox();
    std::vector<MaterialPointP> pointsArray = getMaterialPointsArray();
    int pointSize = pointsArray.size();
    for (int k = 0; k < pointSize; k++)
    {
      MaterialPointP cur_point = pointsArray[k];
      myfile << " Node:         " << cur_point->getID() << std::endl;
      myfile << "  " << pointsArray[k]->getPosOld().x() << "  "
             << pointsArray[k]->getPosOld().y() << "  "
             << pointsArray[k]->getPosOld().z()  
             << std::endl;
    }

    
    
    myfile.close();
//    ret2 = chdir(currentFolder.c_str());

  }



void
ComplicatedGeometry::boundaryConditionNodesID(const std::string& nodesFileName,
                                              std::vector <int>& nodesIDArray)
{
  nodesIDArray.clear();
  std::ifstream infile(nodesFileName); 
  if (infile.is_open())
  {
    for (std::string line; std::getline(infile, line); )
    {  
      std::string name;
      int nodeID;
      
      if (!(std::istringstream(line) >> name >> nodeID)) continue; 
      else 
      {
        nodesIDArray.emplace_back(nodeID);    
      }

    }
  }
  else
  {
    std::cout << "Unable to open the file" <<std::endl;
  }


//  if (nodesIDArray == getTractionForceNodes())
//  if (nodesIDArray == getPointForcePoints())
//  {
//    int size = getPointForcePoints().size();
//    for (int i = 0; i < size; i++)
//    {
//      std::cout << getPointForcePoints()[i] << " ";
//      if (i == size - 1) std::cout << std::endl << std::endl;
//    }
//  }

//  if (nodesIDArray == getDirichletBoundaryNodes())
//  {
//    for (int j = 0; j < getDirichletBoundaryNodes().size(); j++)
//    {
//      std::cout << getDirichletBoundaryNodes()[j] << " ";
//      if (j == getDirichletBoundaryNodes().size() - 1) std::cout << std::endl << std::endl;
//    } 
//  }


}



bool
ComplicatedGeometry::isNodeNumberInArray(const NodeP& node, const std::vector <int>& nodesNumber)
{
  int size = nodesNumber.size();
  for (int i = 0; i < size; i++)
  {
    if (node->getID() == nodesNumber[i]) return true;

  }
  
  return false;
}


bool
ComplicatedGeometry::isPointNumberInArray(const MaterialPointP& point, const std::vector <int>& pointsNumber)
{
  int size = pointsNumber.size();
  for (int i = 0; i < size; i++)
  {
    if (point->getID() == pointsNumber[i]) return true;

  }
  
  return false;
}



void
ComplicatedGeometry::createTractionForce()
{
  NodePArray nodes;
  nodes.clear();
  NodePArray nodesArray = Body::getNodes();
  int size = nodesArray.size();
  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];
    if (isNodeNumberInArray(cur_node, getTractionForceNodes()))
    {
      cur_node->surfaceForce(getTractionForce());
      cur_node->onSurface(true);
//      std::cout << "Node" << cur_node->getID() << " traction surface= ("
//                << cur_node->surfaceForce().x() << ", " << cur_node->surfaceForce().y() << ", "
//                << cur_node->surfaceForce().z() << ")" << std::endl;
    }
               
    nodes.emplace_back(cur_node);

  }
  Body::setNodes(nodes);
  std::cout << Body::getNodes().size() << std::endl;

//      std::cout << "Node " << Body::getNodes()[64]->getID() << " traction surface= ("
//                << Body::getNodes()[64]->surfaceForce().x() << ", " 
//                << Body::getNodes()[64]->surfaceForce().y() << ", "
//                << Body::getNodes()[64]->surfaceForce().z() << ")" << std::endl;

//      std::cout << "Node " << Body::getNodes()[123]->getID() << " traction surface= ("
//                << Body::getNodes()[123]->surfaceForce().x() << ", " 
//                << Body::getNodes()[123]->surfaceForce().y() << ", "
//                << Body::getNodes()[123]->surfaceForce().z() << ")" << std::endl;

}


void
ComplicatedGeometry::createRefinedTractionForce()
{
  NodePArray nodes;
  nodes.clear();
  NodePArray nodesArray = getRefinedNodesArray();
  int size = nodesArray.size();
  for (int i = 0; i < size; i++)
  {
    NodeP cur_node = nodesArray[i];
    if (isNodeNumberInArray(cur_node, getTractionForceNodes()))
    {
      cur_node->surfaceForce(getTractionForce());
//      std::cout << "Node" << cur_node->getID() << " traction surface= ("
//                << cur_node->surfaceForce().x() << ", " << cur_node->surfaceForce().y() << ", "
//                << cur_node->surfaceForce().z() << ")" << std::endl;
    }
               
    nodes.emplace_back(cur_node);

  }
  setRefinedNodesArray(nodes);
//  std::cout << getRefinedNodesArray().size() << std::endl;

//      std::cout << "Node " << getRefinedNodesArray()[64]->getID() << " traction surface= ("
//                << getRefinedNodesArray()[64]->surfaceForce().x() << ", " 
//                << getRefinedNodesArray()[64]->surfaceForce().y() << ", "
//                << getRefinedNodesArray()[64]->surfaceForce().z() << ")" << std::endl;

//      std::cout << "Node " << getRefinedNodesArray()[123]->getID() << " traction surface= ("
//                << getRefinedNodesArray()[123]->surfaceForce().x() << ", " 
//                << getRefinedNodesArray()[123]->surfaceForce().y() << ", "
//                << getRefinedNodesArray()[123]->surfaceForce().z() << ")" << std::endl;

}


void
ComplicatedGeometry::createPointForce()
{
  Vector3D zero(0.0);
  MaterialPointPArray pointsArray = getMaterialPointsArray();
  int size = pointsArray.size();
  for (int i = 0; i < size; i++)
  {
    MaterialPointP cur_point = pointsArray[i];
    if (isPointNumberInArray(cur_point, getTractionForceNodes()))
    {
      cur_point->setPointForce(getPointForce());

//      if (!(cur_point->getPointForce() == zero))
//      {
//        std::cout << "ID= " << cur_point->getID() << " PointForce= " << cur_point->getPointForce() << std::endl;
//      }  
    
    }
  }

}


void 
ComplicatedGeometry::calculateNodesAreaOfElement(SymmetricMaterialTrilinearElementSP& element)
{
  Vector3D zero(0.0, 0.0, 0.0);
  TrilinearVolumeElementSP elem = element->triElement();
  std::array <NodeP, 8> nodes = elem->getArray();
  NodeP node1 = nodes[0];
//  if (!(node1->surfaceForce() == zero))
//     std::cout << element->id() << " node1 " << node1->getID() << std::endl; 
  std::vector <double>  node1AreaArray = node1->getAreaNodes(); 
  Point3D pos1 = node1->position();
  Vector3D node1Pos(pos1.x(), pos1.y(), pos1.z());
  NodeP node2 = nodes[1];
//  if (!(node2->surfaceForce() == zero))
//     std::cout << element->id() << " node2 " << node2->getID() << std::endl; 
  std::vector <double>  node2AreaArray = node2->getAreaNodes(); 
  Point3D pos2 = node2->position();
  Vector3D node2Pos(pos2.x(), pos2.y(), pos2.z());
  NodeP node3 = nodes[2];
//  if (!(node3->surfaceForce() == zero))
//     std::cout << element->id() << " node3 " << node3->getID() << std::endl; 
  std::vector <double>  node3AreaArray = node3->getAreaNodes(); 
  Point3D pos3 = node3->position();
  Vector3D node3Pos(pos3.x(), pos3.y(), pos3.z());
  NodeP node4 = nodes[3];
//  if (!(node4->surfaceForce() == zero))
//     std::cout << element->id() << " node4 " << node4->getID() << std::endl; 
  std::vector <double>  node4AreaArray = node4->getAreaNodes(); 
  Point3D pos4 = node4->position();
  Vector3D node4Pos(pos4.x(), pos4.y(), pos4.z());
  NodeP node5 = nodes[4];
//  if (!(node5->surfaceForce() == zero))
//     std::cout << element->id() << " node5 " << node5->getID() << std::endl; 
  std::vector <double>  node5AreaArray = node5->getAreaNodes(); 
  Point3D pos5 = node5->position();
  Vector3D node5Pos(pos5.x(), pos5.y(), pos5.z());
  NodeP node6 = nodes[5]; 
//  if (!(node6->surfaceForce() == zero))
//     std::cout << element->id() << " node6 " << node6->getID() << std::endl;
  std::vector <double>  node6AreaArray = node6->getAreaNodes(); 
  Point3D pos6 = node6->position();
  Vector3D node6Pos(pos6.x(), pos6.y(), pos6.z());
  NodeP node7 = nodes[6];
//  if (!(node7->surfaceForce() == zero))
//     std::cout << element->id() << " node7 " << node7->getID() << std::endl; 
  std::vector <double>  node7AreaArray = node7->getAreaNodes(); 
  Point3D pos7 = node7->position();
  Vector3D node7Pos(pos7.x(), pos7.y(), pos7.z());
  NodeP node8 = nodes[7];
//  if (!(node8->surfaceForce() == zero))
//     std::cout << element->id() << " node8 " << node8->getID() << std::endl << std::endl; 
  std::vector <double>  node8AreaArray = node8->getAreaNodes(); 
  Point3D pos8 = node8->position();
  Vector3D node8Pos(pos8.x(), pos8.y(), pos8.z());

  double area1 = 0.0; double area2 = 0.0;
  double area3 = 0.0; double area4 = 0.0;
  double area5 = 0.0; double area6 = 0.0;
  double area7 = 0.0; double area8 = 0.0;


  if (elem->face(5, 6, 7, 8)) // xi3 = 1
  {
    area5 = nodeArea(findGlobalPosition(element, -1, 0, 1),
                       findGlobalPosition(element, 0, 0, 1),  
                       findGlobalPosition(element, 0, 1, 1),
                       node5Pos);

  
    area6 = nodeArea(node6Pos,
                       findGlobalPosition(element, 0, -1, 1),  
                       findGlobalPosition(element, 0, 0, 1),
                       findGlobalPosition(element, -1, 0, 1));


    area7 = nodeArea(findGlobalPosition(element, 0, -1, 1),
                       node7Pos,  
                       findGlobalPosition(element, 1, 0, 1),
                       findGlobalPosition(element, 0, 0, 1));


    area8 = nodeArea(findGlobalPosition(element, 0, 0, 1),  
                       findGlobalPosition(element, 1, 0, 1),
                       node8Pos,
                       findGlobalPosition(element, 0, 1, 1));
  }


  if (elem->face(1, 2, 3, 4)) // xi3 = -1
  {
    area1 = nodeArea(node1Pos,
                       findGlobalPosition(element, 0, -1, -1),
                       findGlobalPosition(element, 0, 0, -1),  
                       findGlobalPosition(element, -1, 0, -1));

  
    area2 = nodeArea(findGlobalPosition(element, 0, -1, -1),
                       node2Pos,  
                       findGlobalPosition(element, 1, 0, -1),
                       findGlobalPosition(element, 0, 0, -1));


    area3 = nodeArea(findGlobalPosition(element, 0, 0, -1),  
                       findGlobalPosition(element, 1, 0, -1),
                       node3Pos,
                       findGlobalPosition(element, 0, 1, -1));


    area4 = nodeArea(findGlobalPosition(element, -1, 0, 1),  
                       findGlobalPosition(element, 0, 0, -1),
                       findGlobalPosition(element, 0, 1, -1),
                       node4Pos);
  }


  if (elem->face(3, 4, 5, 8)) // xi2 = 1
  {
    area3 += nodeArea(findGlobalPosition(element, 0, 1, 0),
                        findGlobalPosition(element, 1, 1, 0),
                        node3Pos,  
                        findGlobalPosition(element, 0, 1, -1));

  
    area4 += nodeArea(findGlobalPosition(element, -1, 1, 0),  
                        findGlobalPosition(element, 0, 1, 0),
                        findGlobalPosition(element, 0, 1, -1),
                        node4Pos);


    area5 += nodeArea(node5Pos,
                        findGlobalPosition(element, 0, 1, 1),  
                        findGlobalPosition(element, 0, 1, 0),
                        findGlobalPosition(element, -1, 1, 0));


    area8 += nodeArea(findGlobalPosition(element, 0, 1, 1),
                        node8Pos,  
                        findGlobalPosition(element, 1, 1, 0),
                        findGlobalPosition(element, 0, 1, 0));
  }


  if (elem->face(1, 2, 6, 7)) // xi2 = -1
  {
    area1 += nodeArea(findGlobalPosition(element, -1, -1, 0),
                        findGlobalPosition(element, 0, -1, 0), 
                        findGlobalPosition(element, 0, -1, -1),
                        node1Pos);

  
    area2 += nodeArea(findGlobalPosition(element, 0, -1, 0),  
                        findGlobalPosition(element, 1, -1, 0),
                        node2Pos,
                        findGlobalPosition(element, 0, -1, -1));


    area6 += nodeArea(node6Pos,
                        findGlobalPosition(element, 0, -1, 1),  
                        findGlobalPosition(element, 0, -1, 0),
                        findGlobalPosition(element, -1, -1, 0));


    area7 += nodeArea(findGlobalPosition(element, 0, -1, 1),
                        node7Pos,  
                        findGlobalPosition(element, 1, -1, 0),
                        findGlobalPosition(element, 0, -1, 0));
  }


  if (elem->face(2, 3, 7, 8)) // xi1 = 1
  {
    area2 += nodeArea(findGlobalPosition(element, 1, -1, 0),
                        node2Pos,
                        findGlobalPosition(element, 1, 0, -1), 
                        findGlobalPosition(element, 1, 0, 0));

  
    area3 += nodeArea(findGlobalPosition(element, 1, 0, 0),  
                        findGlobalPosition(element, 1, 0, -1),
                        node3Pos,
                        findGlobalPosition(element, 1, 1, 0));


    area7 += nodeArea(node7Pos,
                        findGlobalPosition(element, 1, -1, 0),  
                        findGlobalPosition(element, 1, 0, 0),
                        findGlobalPosition(element, 1, 0, 1));


    area8 += nodeArea(findGlobalPosition(element, 1, 0, 1),  
                        findGlobalPosition(element, 1, 0, 0),
                        findGlobalPosition(element, 1, 1, 0),
                        node8Pos);
  }


  if (elem->face(1, 4, 5, 6)) // xi1 = -1
  {
    area1 += nodeArea(findGlobalPosition(element, -1, -1, 0),
                        node1Pos,
                        findGlobalPosition(element, -1, 0, -1), 
                        findGlobalPosition(element, -1, 0, 0));

  
    area4 += nodeArea(findGlobalPosition(element, -1, 0, 0),  
                        findGlobalPosition(element, -1, 0, -1),
                        node4Pos,
                        findGlobalPosition(element, -1, 1, 0));


    area5 += nodeArea(findGlobalPosition(element, -1, 0, 1),  
                        findGlobalPosition(element, -1, 0, 0),
                        findGlobalPosition(element, -1, 1, 0),
                        node5Pos);


    area6 += nodeArea(node6Pos,
                        findGlobalPosition(element, -1, -1, 0),  
                        findGlobalPosition(element, -1, 0, 0),
                        findGlobalPosition(element, -1, 0, 1));
  }

  if (area1 != 0)
    node1AreaArray.emplace_back(area1);
  node1->setAreaNodes(node1AreaArray);

  if (area2 != 0)
    node2AreaArray.emplace_back(area2);
  node2->setAreaNodes(node2AreaArray);

  if (area3 != 0)
    node3AreaArray.emplace_back(area3);
  node3->setAreaNodes(node3AreaArray);

  if (area4 != 0)
    node4AreaArray.emplace_back(area4);
  node4->setAreaNodes(node4AreaArray);

  if (area5 != 0)
    node5AreaArray.emplace_back(area5);
  node5->setAreaNodes(node5AreaArray);

  if (area6 != 0)
     node6AreaArray.emplace_back(area6);
  node6->setAreaNodes(node6AreaArray);

  if (area7 != 0)
    node7AreaArray.emplace_back(area7);
  node7->setAreaNodes(node7AreaArray);

  if (area8 != 0)
    node8AreaArray.emplace_back(area8);
  node8->setAreaNodes(node8AreaArray);

  elem->setArray(node1, node2, node3, node4,
                 node5, node6, node7, node8);
  element->triElement(elem);

}


double
ComplicatedGeometry::nodeArea(Vector3D v1, Vector3D v2, Vector3D v3, Vector3D v4)
{
  double a = ((v4 - v1).cross(v2 - v1)).length();
  double b = ((v2 - v3).cross(v4 - v3)).length();

  return (0.5*(a + b));

}



void
ComplicatedGeometry::calculateNodesArea()
{
  std::vector <SymmetricMaterialTrilinearElementSP> elementsArray = getRefinedElementsArray();
  int size = elementsArray.size();
  for (int i = 0; i < size; i++)
  {
    SymmetricMaterialTrilinearElementSP cur_element = elementsArray[i];
    calculateNodesAreaOfElement(cur_element);
  }

  std::vector <NodeP> nodes;
  MaterialPointPArray points;  
  nodes.clear();
  points.clear();
  std::vector <NodeP>  nodesArray = getRefinedNodesArray();
  MaterialPointPArray pointsArray = getMaterialPointsArray();
  int nodeSize = nodesArray.size();  
  for (int j = 0; j < nodeSize; j++)
  {
    NodeP cur_node = nodesArray[j];
    MaterialPointP cur_point = pointsArray[j];
    std::vector <double> areasArray = cur_node->getAreaNodes();
//    if (areasArray.size() != 0)
//      std::cout << cur_node->getID() << " " << areasArray.size() << std::endl;
    double area = 0.0;
    int areaSize = areasArray.size();
    for (int k = 0; k < areaSize; k++)
    {
      area = area + areasArray[k];
    }

//    if (area != 0)
//      std::cout << cur_node->getID() << " " << area << std::endl;
      
   cur_node->area(area);
   cur_point->setInitialArea(area);
   nodes.emplace_back(cur_node);
   points.emplace_back(cur_point); 
  }
  
  setRefinedNodesArray(nodes);
  setMaterialPointsArray(points);

 
//  double totalVolume = 0.0;
//  std::vector <NodeP>  nodess = getRefinedNodesArray(); 
//  int siz = nodess.size();
//  for (int kk = 0; kk < siz; kk++)
//  {
//    std::cout << "Node ID= " << nodess[kk]->getID() << " volume= " << nodess[kk]->volume() << std::endl;
//    totalVolume += nodess[kk]->volume();
//  }

//  std::cout << "Total Volume= " << totalVolume << std::endl;




  double totalVolume = 0.0;
  MaterialPointPArray  pointss = getMaterialPointsArray(); 
  int siz = pointss.size();
  for (int kk = 0; kk < siz; kk++)
  {
//   if (pointss[kk]->getInitialArea() != 0.0)
//    {
//      std::cout << "Point ID= " << pointss[kk]->getID() << " area= " << pointss[kk]->getInitialArea()
//                                << "   volume= " << pointss[kk]->getInitialVolume() << std::endl;
//    }
    totalVolume += pointss[kk]->getInitialVolume();
  }

//  std::cout << "Total Volume= " << totalVolume << std::endl;

}



void 
ComplicatedGeometry::modifyNodesID()//const std::string& nodesFileName)
{
//  int size = Body::getNodes().size();
  std::vector <int> newNodesID;
  std::vector <int> tractionNodesID = getTractionForceNodes();
//  std::vector <int> tractionNodesID = getDirichletBoundaryNodes();
  std::vector <NodeP> nodes = Body::getNodes();
  std::vector <NodeP> refinedNodes = getRefinedNodesArray();
  for (int i = 0; i < getRefinedNodesArray().size(); i++)
  {
    NodeP refined_node = refinedNodes[i];
    for (int j = 0; j < getTractionForceNodes().size(); j++)
//    for (int j = 0; j < getDirichletBoundaryNodes().size(); j++)
    {
//      std::cout << nodes[tractionNodesID[j]-1]->getID() << " ";
      if (dist(nodes[tractionNodesID[j]-1], refined_node) < 0.000001)
      {
        newNodesID.emplace_back(refined_node->getID());
        continue;
      }     

    }
//    std::cout << std::endl;
  }
  
//  for (int j = 0; j < d_traction_force_nodes.size(); j++)
//  {
//    std::cout << d_traction_force_nodes[j] << "  ";
//  }

//  for (int j = 0; j < d_dirichlet_boundary_nodes.size(); j++)
//  {
//    std::cout << d_dirichlet_boundary_nodes[j] << "  ";
//  }


  d_traction_force_nodes.clear();
//  d_dirichlet_boundary_nodes.clear();


  setTractionForceNodes(newNodesID);
//  setDirichletBoundaryNodes(newNodesID);

//  std::cout << std::endl;

//  for (int k = 0; k < d_traction_force_nodes.size(); k++)
//  {
//    std::cout << d_traction_force_nodes[k] << "  ";
//  }
//  std::cout << std::endl;

//  for (int k = 0; k < d_dirichlet_boundary_nodes.size(); k++)
//  {
//    std::cout << d_dirichlet_boundary_nodes[k] << "  ";
//  }
//  std::cout << std::endl;

//  for (int l = 0; l < refinedNodes.size(); l++)
//  {
//    std::cout << refinedNodes[l]->getID() << " " << refinedNodes[l]->x() << " "
//                                          << refinedNodes[l]->y() << " "
//                                          << refinedNodes[l]->z() << std::endl;

//  }
//  std::cout << std::endl; 

}


void
ComplicatedGeometry::modifyPointsID(std::vector <int>& pointsID)
{
  std::vector <int> newPointsID;
  std::vector <MaterialPointP> refinedPoints = getMaterialPointsArray();
  std::vector <NodeP> nodes = Body::getNodes();
  for (int i = 0; i < getMaterialPointsArray().size(); i++)
  {
    MaterialPointP refined_point = refinedPoints[i];
    for (int j = 0; j < pointsID.size(); j++)
    {
//      std::cout << nodes[tractionNodesID[j]-1]->getID() << " ";
      if (dist(nodes[pointsID[j]-1], refined_point) < 0.000001)
      {
        newPointsID.emplace_back(refined_point->getID());
        continue;
      }     

    }

  }





//  for (int j = 0; j < d_dirichlet_boundary_nodes.size(); j++)
//  {
//    std::cout << d_dirichlet_boundary_nodes[j] << "  ";
//  }


  pointsID.clear();
  int size = newPointsID.size();
  for (int k =0; k < size; k++)
  {
    pointsID.emplace_back(newPointsID[k]);
  }


//  if (pointsID == getPointForcePoints())
//  {
//    for (int j = 0; j < getPointForcePoints().size(); j++)
//    {
//      std::cout << getPointForcePoints()[j] << " ";
//    }
//    std::cout << std::endl; 
//  }


//  for (int j = 0; j < d_dirichlet_boundary_nodes.size(); j++)
//  {
//    std::cout << d_dirichlet_boundary_nodes[j] << "  ";
//  }


}


  



double
ComplicatedGeometry::dist(NodeP node1, NodeP node2)
{
  double x = node1->x() - node2->x();
  double y = node1->y() - node2->y();
  double z = node1->z() - node2->z();
  double distance = std::sqrt(x*x + y*y + z*z);
  return distance;


}


double
ComplicatedGeometry::dist(NodeP node, MaterialPointP point)
{
  double x = node->x() - point->x();
  double y = node->y() - point->y();
  double z = node->z() - point->z();
  double distance = std::sqrt(x*x + y*y + z*z);
  return distance;


}

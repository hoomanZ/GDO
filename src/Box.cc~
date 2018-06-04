#include <Box.h>
#include <Node.h>
#include <SymmetricMaterialTrilinearElementSP.h>
#include <SymmetricMaterialTrilinearElement.h>


#include <MaterialPointP.h>  // for MPM
#include <MaterialPoint.h>   // for MPM


using namespace FiniteElement;


  Box::Box()
    : Body::Body(), d_num_elements_x(1),  
     d_num_elements_y(1), d_num_elements_z(1), 
     d_min_point(0.0), d_max_point(0.0),
     d_surface_force(0.0), d_body_force(0.0),
     d_point_force(0.0)

 //    d_num_points_x(0), d_num_points_y(0), d_num_points_z(0)  //  for implementing MPM
  {
//    d_surface = none;
//    d_point_surface = none;
//    createNodesArray();
//    createElementsArray();
//    Body::makeGlobalMatrices();
  }

 Box::Box(Vector3D minPoint, Vector3D maxPoint,
         int xNumElements, int yNumElements, int zNumElements)
    : Body::Body(), d_surface_force(0.0), d_body_force(0.0), d_point_force(0.0),
      d_num_points_x(0), d_num_points_y(0), d_num_points_z(0)  //  for implementing MPM
  {
    setTractionSurface(none);
    setPointSurface(none);
    setMinPoint(minPoint);
    setMaxPoint(maxPoint);
    setNumElementsX(xNumElements);
    setNumElementsY(yNumElements);
    setNumElementsZ(zNumElements);
    createNodesArray();
    createElementsArray();
    Body::makeGlobalMatrices();
  }

 Box::Box(Vector3D minPoint, Vector3D maxPoint,
         int xNumElements, int yNumElements, int zNumElements,
         double density, double poission, double young)
     : d_surface_force(0.0), d_body_force(0.0), d_point_force(0.0),
       d_num_points_x(0), d_num_points_y(0), d_num_points_z(0)  //  for implementing MPM
  {
    Body::setProperties(density, young, poission);
    setTractionSurface(none);
    setPointSurface(none);
    setMinPoint(minPoint);
    setMaxPoint(maxPoint);
    setNumElementsX(xNumElements);
    setNumElementsY(yNumElements);
    setNumElementsZ(zNumElements);
    createNodesArray();
    createElementsArray();
    Body::makeGlobalMatrices();
  }



 Box::Box(const Vector3D& minPoint, const Vector3D& maxPoint,
          const int& xNumElements, const int& yNumElements, const int& zNumElements,
          const double& density, const double& poission, const double& young,
          const double& bodyForce, const double& surForce, const ForceSurface& tracSurface,
          const FixedDirichletBoundary& boundary, const Vector3D& boundary_position)
   : d_num_points_x(0), d_num_points_y(0), d_num_points_z(0), d_point_force(0.0)  //  for implementing MPM
   {
 //   std::cout << xNumPoints << " " << yNumPoints << " " << zNumPoints << std::endl;
    Body::setProperties(density, young, poission);
    setMinPoint(minPoint);
    setMaxPoint(maxPoint);
    setNumElementsX(xNumElements);
    setNumElementsY(yNumElements);
    setNumElementsZ(zNumElements);
    setBodyForce(bodyForce);
    setTractionSurface(tracSurface);
    setTractionForce(surForce);
    setPointSurface(none);
    setFixedDirichletBoundary(boundary);
    setBoundaryDisplacement(boundary_position);
    createTractionForce();
//    Body::makeGlobalMatrices();
  }



 Box::Box(const Vector3D& minPoint, const Vector3D& maxPoint,
          const int& xNumElements, const int& yNumElements, const int& zNumElements,
          const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
          const double& density, const double& poission, const double& young,
          const double& bodyForce, const Vector3D& pointForce,
          const double& surForce, const ForceSurface& tracSurface, const ForceSurface& pointSurface,
          const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition)

   {
//    std::cout << xNumPoints << " " << yNumPoints << " " << zNumPoints << std::endl;
    Body::setProperties(density, young, poission);
    setMinPoint(minPoint);
    setMaxPoint(maxPoint);
    setNumElementsX(xNumElements);
    setNumElementsY(yNumElements);
    setNumElementsZ(zNumElements);
    setNumPointsX(xNumPoints); // MPM
    setNumPointsY(yNumPoints); // MPM
    setNumPointsZ(zNumPoints); // MPM
    setBodyForce(bodyForce);
    setTractionSurface(tracSurface);
    setTractionForce(surForce);
    setPointForce(pointForce);
    setPointSurface(pointSurface);
    setFixedDirichletBoundary(boundary);
    setBoundaryDisplacement(boundaryPosition);
    createTractionForce();
    createMaterialPointsArray();  // MPM 
    findPointsMassAndVolume();  // MPM 
//    Body::makeGlobalMatrices();
  }


 Box::Box(const Vector3D& minPoint, const Vector3D& maxPoint,
          const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
          const double& density, const double& poission, const double& young, const double& yield,
          const double& bodyForce, const double& surForce, const ForceSurface& tracSurface,
          const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition)
     
   {
//    std::cout << xNumPoints << " " << yNumPoints << " " << zNumPoints << std::endl;
    Body::setProperties(density, young, poission, yield);
    setMinPoint(minPoint);
    setMaxPoint(maxPoint);
    setNumElementsX(0.0);
    setNumElementsY(0.0);
    setNumElementsZ(0.0);
    setNumPointsX(xNumPoints); // MPM
    setNumPointsY(yNumPoints); // MPM
    setNumPointsZ(zNumPoints); // MPM
    setBodyForce(bodyForce);
    setTractionSurface(tracSurface);
    setTractionForce(surForce);
//    setPointForce(pointForce);
//    setPointSurface(pointSurface);
    setFixedDirichletBoundary(boundary);
    setBoundaryDisplacement(boundaryPosition);
//    createTractionForce();
//    std::cout << getNumPointsX() << " " << getNumPointsY() << " " << getNumPointsZ() << std::endl;
    createMaterialPointsArray();  // MPM & Peridynamics
    findPointsMassAndVolume();  // MPM & Peridynamics
//    Body::makeGlobalMatrices();
  }



 Box::Box(BoxSP box)
  {
    Body::setProperties(box->getDensity(), box->getPoission(), box->getYoung());
    setMinPoint(box->getMinPoint());
    setMaxPoint(box->getMaxPoint());
    setNumElementsX(box->getNumElementsX());
    setNumElementsY(box->getNumElementsY());
    setNumElementsZ(box->getNumElementsZ());
    setNumPointsX(box->getNumPointsX()); // MPM
    setNumPointsY(box->getNumPointsY()); // MPM
    setNumPointsZ(box->getNumPointsZ()); // MPM
  }

/* Box::Box(Box box)
  {
    Body::setProperties(box.getDensity(), box.getPoission(), box.getYoung());
    setMinPoint(box.getMinPoint());
    setMaxPoint(box.getMaxPoint());
    setNumElementsX(box.getNumElementsX());
    setNumElementsY(box.getNumElementsY());
    setNumElementsZ(box.getNumElementsZ());
  }*/

 Box::~Box()
 {
 }



 int
 Box::nodeID(int i, int j, int k)
   {
     int nx = d_num_elements_x;
     int ny = d_num_elements_y;
//     int nz = d_num_elements_z;
     return k*(ny+1)*(nx+1)+j*(nx+1)+i+1;
   }


 int
 Box::elementID(int i, int j, int k)
   {
     int nx = d_num_elements_x;
     int ny = d_num_elements_y;
//     int nz = d_num_elements_z;
     return k*ny*nx+j*nx+i+1;
   }


 void
 Box::createNodesArray()
   {
     Body::clearMass();
     std::vector <NodeP> nodesArray;
     nodesArray.clear();

     int nx = d_num_elements_x;
     int ny = d_num_elements_y;
     int nz = d_num_elements_z;

//     std::cout << nx << " " << ny << " " << nz << std::endl;

     bool surfaceOnX;
     bool surfaceOnY;
     bool surfaceOnZ;

     double dx = (d_max_point.x()-d_min_point.x())/nx;
     double dy = (d_max_point.y()-d_min_point.y())/ny;
     double dz = (d_max_point.z()-d_min_point.z())/nz;

     for (int k=0; k < nz+1; k++)
     {
       for (int j=0; j < ny+1; j++)
       {
         for (int i=0; i < nx+1; i++)
         {
           NodeP node(new Node(nodeID(i, j, k),
                      d_min_point.x()+i*dx,
                      d_min_point.y()+j*dy,
                      d_min_point.z()+k*dz, false));
           surfaceOnX = ((node->x() == d_min_point.x()) || (node->x() == d_max_point.x()));
           surfaceOnY = ((node->y() == d_min_point.y()) || (node->y() == d_max_point.y()));
           surfaceOnZ = ((node->z() == d_min_point.z()) || (node->z() == d_max_point.z()));
           node->onSurface(surfaceOnX || surfaceOnY || surfaceOnZ);
           node->densityNode(Body::getDensity());
//           std::cout << "NodeID= " << node->getID() << " Velocity= " << node->velocity() << std::endl;
           nodesArray.emplace_back(node);
         }
       }
     }

//     std::cout << nodesArray.size() << std::endl;

     Body::setNodes(nodesArray);

   }

 

 void
 Box::createElementsArray()
   {
     Body::clearStiff();
     std::vector <SymmetricMaterialTrilinearElementSP> elementsArray;
     elementsArray.clear();

     int nx = d_num_elements_x;
     int ny = d_num_elements_y;
     int nz = d_num_elements_z;

     int n1, n2, n3, n4;
     int n5, n6, n7, n8;

     double density = Body::getDensity();
     double poission = Body::getPoission();
     double young = Body::getYoung();

     for (int k=0; k < nz; k++)
     {
       for (int j=0; j < ny; j++)
       {
         for (int i=0; i < nx; i++)
         {
           n1 = nodeID(i, j, k);
           n2 = nodeID(i+1, j, k);
           n3 = nodeID(i+1, j+1, k);//nodeID(i, j+1, k);
           n4 = nodeID(i, j+1, k);//nodeID(i+1, j+1, k);
           n5 = nodeID(i, j+1, k+1);//nodeID(i, j, k+1);
           n6 = nodeID(i, j, k+1);//nodeID(i+1, j, k+1);
           n7 = nodeID(i+1, j, k+1);//nodeID(i, j+1, k+1);
           n8 = nodeID(i+1, j+1, k+1);
           SymmetricMaterialTrilinearElementSP element(new SymmetricMaterialTrilinearElement
                                                      (elementID(i, j, k), density, young, poission, Body::getNodes(), 
                                                      n1, n2, n3, n4, n5, n6, n7, n8));
                                                        
           elementsArray.push_back(element);
         }
       }
     }
    Body::setElements(elementsArray);
  }


 void
 Box::createTractionForce()
   {
     createNodesArray();

     unsigned int size = Body::getNodes().size();
     ForceSurface  trac_surf = getTractionSurface();
     Vector3D force_vec(0.0);
     double force = getTractionForce();
     for (unsigned int i = 0; i < size; i++)
       {     
          NodeP cur_node = Body::getNodes()[i];
          switch (trac_surf)
          {
          case right:
            if ((cur_node->onSurface()) && (cur_node->x() == getMaxPoint().x()))
            {
//              force_vec.x(-1*force);
              force_vec.x(force);
              cur_node->surfaceForce(force_vec);
            }
           break;

          case left:
            if ((cur_node->onSurface()) && (cur_node->x() == getMinPoint().x()))
            {
//              force_vec.x(force);
              force_vec.x(-1*force);
              cur_node->surfaceForce(force_vec);
            }
           break;

          case up:
            if ((cur_node->onSurface()) && (cur_node->y() == getMaxPoint().y()))
            {
//              force_vec.y(-1*force);
              force_vec.y(force);
              cur_node->surfaceForce(force_vec);
            }
           break;

          case down:
            if ((cur_node->onSurface()) && (cur_node->y() == getMinPoint().y()))
            {
//              force_vec.y(force);
              force_vec.y(-1*force);
              cur_node->surfaceForce(force_vec);
            }
           break;

          case front:
            if ((cur_node->onSurface()) && (cur_node->z() == getMaxPoint().z()))
            {
//              force_vec.z(-1*force);
              force_vec.z(force);
              cur_node->surfaceForce(force_vec);
            }
           break;

          case back:
            if ((cur_node->onSurface()) && (cur_node->z() == getMinPoint().z()))
            {
//              force_vec.z(force);
              force_vec.z(-1*force);
              cur_node->surfaceForce(force_vec);
            }
          break;

          case upDownY:
            if ((cur_node->onSurface()) && (cur_node->y() == getMinPoint().y()))
            {
//              force_vec.z(force);
              force_vec.z(-1*force);
              cur_node->surfaceForce(force_vec);
            }
            else if ((cur_node->onSurface()) && (cur_node->y() == getMaxPoint().y()))
            {
              force_vec.z(force);
//              force_vec.z(-1*force);
              cur_node->surfaceForce(force_vec);
            }
          break;


          case none:
              cur_node->surfaceForce(force_vec);
          break;

          default:
            std::cout <<
          " The surface can be right, left, up, down, front, back, or none which means we do not have any traction force."
            << std::endl;
            break;
          }
       }
     createElementsArray();                                 
    }


void 
Box::createFixedDirichletBoundary()
  {
    
  }


/*********************************
 *
 * Creating material points for 
 *      implementing MPM
 *
 ********************************/


int
Box::pointID(int i, int j, int k)
   {
     int nx = d_num_points_x;
     int ny = d_num_points_y;
//     int nz = d_num_elements_z;
//     return k*(ny+1)*(nx+1)+j*(nx+1)+i+1;
     return k*ny*nx+j*nx+i+1;
   }



void
Box::createMaterialPointsArray()
{
     std::vector <MaterialPointP> pointsArray;
     pointsArray.clear();

     int nx = d_num_points_x;
     int ny = d_num_points_y;
     int nz = d_num_points_z;

//     std::cout << nx << " " << ny << " " << nz << std::endl;

     double dx = (d_max_point.x()-d_min_point.x())/(nx-1);
     double dy = (d_max_point.y()-d_min_point.y())/(ny-1);
     double dz = (d_max_point.z()-d_min_point.z())/(nz-1);

     for (int k=0; k < nz; k++)
     {
       for (int j=0; j < ny; j++)
       {
         for (int i=0; i < nx; i++)
         {
           MaterialPointP point(new MaterialPoint(pointID(i, j, k),
                                    d_min_point.x()+i*dx,
                                    d_min_point.y()+j*dy,
                                    d_min_point.z()+k*dz));
//           surfaceOnX = ((node->x() == d_min_point.x()) || (node->x() == d_max_point.x()));
//           surfaceOnY = ((node->y() == d_min_point.y()) || (node->y() == d_max_point.y()));
//           surfaceOnZ = ((node->z() == d_min_point.z()) || (node->z() == d_max_point.z()));
//           node->onSurface(surfaceOnX || surfaceOnY || surfaceOnZ);
//           node->densityNode(Body::getDensity());
           pointsArray.emplace_back(point);
         }
       }
     }
     Body::setMaterialPoints(pointsArray);

//    std::cout << Body::getMaterialPoints().size() << std::endl;

//     printPointsPosition();
}


void
Box::printPointsPosition()
{
  std::vector <MaterialPointP> pointsArray = getMaterialPoints();
  int size = pointsArray.size(); 
  for (int i = 0; i < size; i++)
  {
    MaterialPointP cur_point = pointsArray[i];
    std::cout << "Position of point number " << cur_point->getID() << " is: ("
              << cur_point->getPosOld().x() << ", "
              << cur_point->getPosOld().y() << ", "
              << cur_point->getPosOld().z() << ")" << std::endl;
  }

}


void 
Box::findPointsMassAndVolume()
{
  double sum = 0.0;
  int condition = 8;
  double density = Body::getDensity();
  double volume = 0;
  volume = (d_max_point.x() - d_min_point.x())*(d_max_point.y() - d_min_point.y())*(d_max_point.z() - d_min_point.z());
  double unit_mass = volume*density/(8*(d_num_points_x-1)*(d_num_points_y-1)*(d_num_points_z-1));
  std::vector <MaterialPointP> pointsArray = getMaterialPoints();
  int size = pointsArray.size(); 
  for (int i = 0; i < size; i++)
  {
    MaterialPointP cur_point = pointsArray[i]; 
    if ((cur_point->getPosOld().x() == d_max_point.x()) || (cur_point->getPosOld().x() == d_min_point.x()))
        condition /= 2;        
    if ((cur_point->getPosOld().y() == d_max_point.y()) || (cur_point->getPosOld().y() == d_min_point.y()))
        condition /= 2;        
    if ((cur_point->getPosOld().z() == d_max_point.z()) || (cur_point->getPosOld().z() == d_min_point.z()))
        condition /= 2;
    cur_point->setMass(condition*unit_mass);
    cur_point->setInitialVolume(condition*unit_mass/density);

    sum += condition*unit_mass;

//    std::cout << cur_point->getID() << "   " << "Condition: " << condition << "  Mass: " << cur_point->getMass() 
//              << std::endl;

//    std::cout << cur_point->getID() << "   " << "Condition: " << condition << "  Volume: " 
//              << cur_point->getInitialVolume() << std::endl;


    condition = 8;
              
  }

//  std::cout << "Total mass is: " << sum << std::endl;

}


void
Box::findPointsArea()
{
   





}



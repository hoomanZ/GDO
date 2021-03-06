#ifndef FINITEELEMENT_BOX_H
#define FINITEELEMENT_BOX_H


#include <Body.h>
#include <GeometryMath/Vector3D.h>
#include <BoxSP.h>
#include <string>


namespace FiniteElement {

  class Box : public Body {
    public:

      enum ForceSurface
      {
        right, //xMax
        left,  //xMin
        up,    //yMax
        down,  //yMin
        front, //zMax
        back,  //zMin
        upDownX, //xMaxANDxMin 
        upDownY, //yMaxANDyMin 
        upDownZ, //zMaxANDzMin 
        tear, // Velocity in two opposite directions
        none   //no traction force
      };

      enum FixedDirichletBoundary
      {
        xMax,
        xMin,
        yMax,
        yMin,
        zMax,
        zMin,
        noFixedBoundary
      };

      Box();
      Box(Vector3D minPoint, Vector3D maxPoint, int xNumElements, int yNumElements, int zNumElements);
      Box(Vector3D minPoint, Vector3D maxPoint, int xNumElements, int yNumElements, int zNumElements,
                                                double density, double poission, double young);
      Box(const Vector3D& minPoint, const Vector3D& maxPoint, 
          const int& xNumElements, const int& yNumElements, const int& zNumElements,
          const double& density, const double& poission, const double& young, const double& bodyForce, 
          const double& surForce, const ForceSurface& surface,
          const FixedDirichletBoundary& boundary, const Vector3D& boundary_position);


      Box(const Vector3D& minPoint, const Vector3D& maxPoint,
          const int& xNumElements, const int& yNumElements, const int& zNumElements,
          const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
          const double& density, const double& poission, const double& young,
          const double& bodyForce, const double& surForce, const ForceSurface& surface,
          const FixedDirichletBoundary& boundary, const Vector3D& boundary_position);

//********************************************** Constructor for MPM **************************************** 

      Box(const Vector3D& minPoint, const Vector3D& maxPoint,
          const int& xNumElements, const int& yNumElements, const int& zNumElements,
          const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
          const double& density, const double& poission, const double& young,
          const double& bodyForce, const Vector3D& pointForce,
          const double& surForce, const ForceSurface& tracSurface, const ForceSurface& pointSurface,
          const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition);

//***********************************************************************************************************


//********************************************** Constructor for Peridynamics **************************************** 

      Box(const Vector3D& minPoint, const Vector3D& maxPoint,
          const int& xNumPoints, const int& yNumPoints, const int& zNumPoints,
          const double& density, const double& poission, const double& young,
          const double& yield, const double& bodyForce, 
          const double& surForce, const ForceSurface& tracSurface,
          const FixedDirichletBoundary& boundary, const Vector3D& boundaryPosition);

//***********************************************************************************************************


      Box(BoxSP box);
      
//      Box(Box box);      
      ~Box();



      inline void setMinPoint(const Vector3D& minPoint) {d_min_point = minPoint;}
      inline Vector3D getMinPoint() const {return d_min_point;}

      inline void setMaxPoint(const Vector3D& maxPoint) {d_max_point = maxPoint;}
      inline Vector3D getMaxPoint() const {return d_max_point;}

      inline void setNumElementsX(const int& numElementsX) {d_num_elements_x = numElementsX;}
      inline int getNumElementsX() const {return d_num_elements_x;}

      inline void setNumElementsY(const int& numElementsY) {d_num_elements_y = numElementsY;}
      inline int getNumElementsY() const {return d_num_elements_y;}

      inline void setNumElementsZ(const int& numElementsZ) {d_num_elements_z = numElementsZ;}
      inline int getNumElementsZ() const {return d_num_elements_z;}

      inline void setTractionForce(const double& force) {d_surface_force = force;}
      inline double getTractionForce() const {return d_surface_force;}

      inline void setPointForce(const Vector3D& force) {d_point_force = force;}
      inline Vector3D getPointForce() const {return d_point_force;}


      inline void setTractionSurface(const ForceSurface& surface) {d_surface = surface;}
      inline ForceSurface getTractionSurface() const {return d_surface;}

      inline void setPointSurface(const ForceSurface& surface) {d_point_surface = surface;}
      inline ForceSurface getPointSurface() const {return d_point_surface;}


      inline void setFixedDirichletBoundary(const FixedDirichletBoundary& boundary) {d_boundary = boundary;}
      inline FixedDirichletBoundary getFixedDirichletBoundary() const {return d_boundary;}


      inline void setBodyForce(const double& force) {d_body_force = force;}
      inline double getBodyForce() const {return d_body_force;}

      inline void setBoundaryDisplacement(const Vector3D& dis) {d_boundary_displacement = dis;}
      inline Vector3D getBoundaryDisplacement() const {return d_boundary_displacement;}
                   
      int nodeID(int i, int j, int k);
      int elementID(int i, int j, int k);

      void createNodesArray();
      void createElementsArray();
      void createTractionForce();
      void createFixedDirichletBoundary();

      void initialMatrix(std::vector <std::vector <double> > Matrix, int size);



//*******  funnctions for MPM  *********
      
      inline void setNumPointsX(const int& points) {d_num_points_x = points;}
      inline int getNumPointsX() const {return d_num_points_x;}

      inline void setNumPointsY(const int& points) {d_num_points_y = points;}
      inline int getNumPointsY() const {return d_num_points_y;}

      inline void setNumPointsZ(const int& points) {d_num_points_z = points;}
      inline int getNumPointsZ() const {return d_num_points_z;}

      int pointID(int i, int j, int k);
      void createMaterialPointsArray();

      void findPointsMassAndVolume();
      void findPointsArea();  


      void printPointsPosition();




   private:
     Vector3D d_min_point;
     Vector3D d_max_point;

     int d_num_elements_x;
     int d_num_elements_y;
     int d_num_elements_z;

     double d_surface_force;
     double d_body_force;
     Vector3D d_boundary_displacement;
     Vector3D d_point_force;

     ForceSurface d_surface;
     ForceSurface d_point_surface;
     FixedDirichletBoundary d_boundary;


//*******  variables for MPM  *********

     int d_num_points_x;
     int d_num_points_y;
     int d_num_points_z;
     
     


  }; //end of class
} //end of namespace
#endif

#ifndef FINITEELEMENT_COMPLICATEDGEOMETRY_H
#define FINITEELEMENT_COMPLICATEDGEOMETRY_H


#include <Body.h>
#include <GeometryMath/Vector3D.h>
#include <BoxSP.h>
#include <Box.h>
#include <string>
//#include <vector>


#include <MaterialPointPArray.h>
#include <NodePArray.h>
#include <SymmetricMaterialTrilinearElementSP.h>

namespace FiniteElement {

  class ComplicatedGeometry : public Body {
    public:
      ComplicatedGeometry();
      ComplicatedGeometry(std::string nodeFileName, std::string elementFileName);
      ComplicatedGeometry(double density, double young, double poission,
                          Vector3D tractionForce, Vector3D boundaryDisplacement);                  
      ComplicatedGeometry(std::string nodeFileName, std::string elementFileName, 
                          double density, double young, double poission,
                          Vector3D tractionForce, Vector3D boundaryDisplacement);
      ComplicatedGeometry(std::string nodeFileName, std::string elementFileName,
                          std::string tractionForceNodesFile, std::string dirichletBoundaryNodesFile,
                          double density, double young, double poission,
                          Vector3D tractionForce, Vector3D boundaryDisplacement);
      ComplicatedGeometry(std::string nodeFileName, std::string elementFileName, 
                          double density, double young, double poission,
                          Vector3D tractionForce, Vector3D boundaryDisplacement,
                          int numSai1, int numSai2, int numSai3);
      ComplicatedGeometry(std::string nodeFileName, std::string elementFileName,
                          std::string tractionForceNodesFile, std::string dirichletBoundaryNodesFile, 
                          double density, double young, double poission,
                          Vector3D tractionForce, Vector3D boundaryDisplacement,
                          int numSai1, int numSai2, int numSai3);      
      
      ComplicatedGeometry(std::string nodeFileName, std::string elementFileName,
                          std::string tractionForceNodesFile, std::string dirichletBoundaryNodesFile,
                          double density, double young, double poission,
                          Vector3D tractionForce, Vector3D boundaryDisplacement,
                          int numSai1, int numSai2, int numSai3,
                          double xNumElem, double yNumElem, double zNumElem);

      ~ComplicatedGeometry();      


      inline void setNodesFileName(const std::string& nodeName) {d_nodes_file_name = nodeName;}
      inline std::string getNodesFileName() const {return d_nodes_file_name;}

      inline void setElementsFileName(const std::string& elementName) {d_elements_file_name = elementName;}
      inline std::string getElementsFileName() const {return d_elements_file_name;}

      inline void setMaterialPointsArray(const MaterialPointPArray& pointsArray) {d_points_array = pointsArray;}
      inline MaterialPointPArray getMaterialPointsArray() const {return d_points_array;}

      inline void setNumPointsSai1(const int& numSai1) {d_num_points_sai1 = numSai1;}
      inline int getNumPointsSai1() const {return d_num_points_sai1;}

      inline void setNumPointsSai2(const int& numSai2) {d_num_points_sai2 = numSai2;}
      inline int getNumPointsSai2() const {return d_num_points_sai2;}

      inline void setNumPointsSai3(const int& numSai3) {d_num_points_sai3 = numSai3;}
      inline int getNumPointsSai3() const {return d_num_points_sai3;}

      inline void setRefinedNodesArray(const std::vector <NodeP>& nodesArray) {d_refined_nodes_array = nodesArray;}
      inline std::vector <NodeP> getRefinedNodesArray() const {return d_refined_nodes_array;}

      inline void setRefinedElementsArray(const std::vector <SymmetricMaterialTrilinearElementSP>& elementsArray) {d_refined_elements_array = elementsArray;}
      inline std::vector <SymmetricMaterialTrilinearElementSP> getRefinedElementsArray() const {return d_refined_elements_array;}



      inline void setPointForcePoints(const std::vector <int>& points) {d_point_force_points = points;}
      inline std::vector <int> getPointForcePoints() const {return d_point_force_points;}

      inline void setTractionForceNodes(const std::vector <int>& nodes) {d_traction_force_nodes = nodes;}
      inline std::vector <int> getTractionForceNodes() const {return d_traction_force_nodes;}

      inline void setDirichletBoundaryNodes(const std::vector <int>& nodes) {d_dirichlet_boundary_nodes = nodes;}
      inline std::vector <int> getDirichletBoundaryNodes() const {return d_dirichlet_boundary_nodes;}

      inline void setPointForce(const Vector3D& force) {d_point_force = force;}
      inline Vector3D getPointForce() const {return d_point_force;}


      inline void setTractionForce(const Vector3D& force) {d_traction_force = force;}
      inline Vector3D getTractionForce() const {return d_traction_force;}

      inline void setBoundaryDisplacement(const Vector3D& dis) {d_boundary_displacement = dis;}
      inline Vector3D getBoundaryDisplacement() const {return d_boundary_displacement;}

      inline void setBox(const BoxSP& box) {d_box = box;}
      inline BoxSP getBox() const {return d_box;}


      void createNodesArray();
      void createElementsArray();
      void createPointsArray();
      void createRefinedElements();
      bool pointRepeated(const Vector3D& position, const MaterialPointPArray& pointArray);
      bool nodeRepeated(const Vector3D& position, const std::vector <NodeP>& refinedNodes, int& index);
      int findNodeIndex(const SymmetricMaterialTrilinearElementSP& element,
                        const double& sai1, const double& sai2, const double& sai3);
      Vector3D findGlobalPosition(const SymmetricMaterialTrilinearElementSP& element,
                                  const double& sai1, const double& sai2, const double& sai3);

      void calculateNodesVolumeOfElement(SymmetricMaterialTrilinearElementSP& element);

      double nodeVolume(Vector3D v1, Vector3D v2, Vector3D v3, Vector3D v4,
                        Vector3D v5, Vector3D v6, Vector3D v7, Vector3D v8);

      double nodeArea(Vector3D v1, Vector3D v2, Vector3D v3, Vector3D v4);


      void calculateNodesVolume();
      void pointsFile();
      void nodesFile(std::string fileName);

      void boundaryConditionNodesID(const std::string& nodesFileName,
                                    std::vector <int>& nodesIDArray);

      bool isNodeNumberInArray(const NodeP& node, const std::vector <int>& nodesNumber);
      bool isPointNumberInArray(const MaterialPointP& point, const std::vector <int>& pointsNumber);

      void createTractionForce();
      void createRefinedTractionForce();
      void createPointForce();
 

      void calculateNodesAreaOfElement(SymmetricMaterialTrilinearElementSP& element);
      void calculateNodesArea();

      void modifyNodesID();//const std::string& nodesFileName);
      void modifyPointsID(std::vector <int>& pointsID);
      double dist(NodeP node1, NodeP node2);
      double dist(NodeP node, MaterialPointP point);

   private:
     std::string d_nodes_file_name;
     std::string d_elements_file_name;

     MaterialPointPArray d_points_array;
     int d_num_points_sai1;
     int d_num_points_sai2;     
     int d_num_points_sai3;

     std::vector <NodeP> d_refined_nodes_array;
     std::vector <SymmetricMaterialTrilinearElementSP>  d_refined_elements_array;


     std::vector <int>  d_traction_force_nodes;
     std::vector <int>  d_dirichlet_boundary_nodes;

     std::vector <int>  d_point_force_points;

     Vector3D d_point_force;
     Vector3D d_traction_force;
     Vector3D d_boundary_displacement;

     BoxSP d_box;   // for MPM

  }; // end of class
} // end of namespace

#endif

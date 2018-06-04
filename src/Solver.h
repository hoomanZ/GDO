#ifndef FINITEELEMENT_SOLVER_H
#define FINITEELEMENT_SOLVER_H

#include <Body.h>
#include <BodySP.h>
#include <BoxSP.h>
#include <Box.h>
#include <ComplicatedGeometry.h>
#include <ComplicatedGeometrySP.h>

//#include <eigen/Eigen/Dense>
#include <GeometryMath/Types.h>


namespace FiniteElement {
  class Solver {
    public:

      enum geometryType {Boxes, Arbitrary};
     
      enum solverType {Static, Explicit, Implicit};
      Solver();

      Solver(BoxSP box);
      Solver(solverType type, double time, double delT, BoxSP box);

      Solver(ComplicatedGeometrySP geometry);
      Solver(solverType type, double time, double delT, ComplicatedGeometrySP geometry);

      
      ~Solver();

      inline void setSolverType(const solverType& type) {d_solver_type = type;}
      inline double getSolverType() const {return d_solver_type;}

      inline void setGeometryType(const geometryType& type) {d_geometry_type = type;}
      inline double getGeometryType() const {return d_geometry_type;}


      inline void setTime(const double& time) {d_time = time;}
      inline double getTime() const {return d_time;}
      
      inline void setDelT(const double& delT) {d_del_t = delT;}
      inline double getDelT() const {return d_del_t;}

      inline void setBox(const BoxSP& box) {d_box = box;}
      inline BoxSP getBox() const {return d_box;}

      inline void setGeometry(const ComplicatedGeometrySP& geometry) {d_geometry = geometry;}
      inline ComplicatedGeometrySP getGeometry() const {return d_geometry;}


      void EIBoxSolvers();
      void EIGeometrySolvers();  
      void newDisplacementExplicit(const std::vector<std::vector<double> >& stiffMatrix,
                           const std::vector<std::vector<double> >& massLumpedMatrix,
                           const std::vector<double>& extForce, const std::vector<double>& oldDisp,
                           const std::vector<double>& disp, std::vector<double>& newDisp);
      void newDisplacementImplicit(const std::vector<std::vector<double> >& stiffMatrix,
                           const std::vector<std::vector<double> >& massMatrix,
                           const std::vector<double>& extForce, const std::vector<double>& oldDisp,
                           const std::vector<double>& disp, std::vector<double>& newDisp);
      void newDisplacementStatic(const std::vector<std::vector<double> >& stiffMatrix,
                                 const std::vector<double>& extForce,
                                 std::vector<double>& newDisp);


      void MUMPS(const std::vector<std::vector<double> >& mat, 
                 const std::vector<double>& rhs,
                 std::vector<double>& newDisp);

      std::vector<double> massDotVec(const std::vector<std::vector<double> >& matrix, 
                                     const std::vector<double>& vec);

      void implicitSolver(); 
      void applyBoundaryCondition(std::vector<std::vector<double> >& stiffMatrix);
      void applyBoundaryConditionGeometry(std::vector<std::vector<double> >& stiffMatrix);
      bool condition(const NodeP& node, const Box::FixedDirichletBoundary& boundary,
                     const Vector3D& max, const Vector3D& min); 

      void nodesFile(const int& iter_num, const std::string& currentFolder);

      void updateNodesPos(const std::vector<double>& disp);  




    private:
      double d_time;
      double d_del_t;
      BoxSP d_box;
      ComplicatedGeometrySP d_geometry;
//      BodySP d_box;
      solverType d_solver_type;
      geometryType d_geometry_type;


  }; //end of class
} //end of namespace
#endif

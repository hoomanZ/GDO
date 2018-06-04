#ifndef FINITEELEMENT_SOLVER_H
#define FINITEELEMENT_SOLVER_H

#include <Body.h>
#include <BoxSP.h>
#include <Box.h>
//#include <GeometryMath/Vector3D.h>


namespace FiniteElement {
  class Solver {
    public:
     
      enum solverType {Explicit, Implicit};
      Solver();
      Solver(solverType type, double time, double delT, BoxSP box);
      ~Solver();

      inline void setSolverType(const solverType& type) {d_solver_type = type;}
      inline double getSolverType() const {return d_solver_type;}

      inline void setTime(const double& time) {d_time = time;}
      inline double getTime() const {return d_time;}
      
      inline void setDelT(const double& delT) {d_del_t = delT;}
      inline double getDelT() const {return d_del_t;}

      inline void setBox(const BoxSP& box) {d_box = box;}
      inline BoxSP getBox() const {return d_box;}

      void explicitSolver(); 
      void newDisplacement(const std::vector<std::vector<double> >& stiffMatrix,
                           const std::vector<std::vector<double> >& massLumpedMatrix,
                           const std::vector<double>& extForce, const std::vector<double>& oldDisp,
                           const std::vector<double>& disp, std::vector<double>& newDisp);

      void implicitSolver();    




    private:
      double d_time;
      double d_del_t;
      BoxSP d_box;
      solverType d_solver_type;



  }; //end of class
} //end of namespace
#endif

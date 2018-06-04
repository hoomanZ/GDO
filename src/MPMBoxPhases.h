#ifndef FINITEELEMENT_MPMBOXPHASES_H
#define FINITEELEMENT_MPMBOXPHASES_H


#include <MPMBox.h>
#include <MPMBoxP.h>
#include <BoxSP.h>
#include <Box.h>

#include <GeometryMath/Vector3D.h>
#include <GeometryMath/Point3D.h>
#include <GeometryMath/Matrix3D.h>

#include <string>

namespace FiniteElement {

  class MPMBoxPhases{
    public:

      enum solverType {Lumped, Consistent, LumpedUsingMomentum, ConsistentUsingMomentum};
      enum updateStressType {USL, MUSL, USF};

      MPMBoxPhases();
      MPMBoxPhases(MPMBoxP box, solverType solver, double time, double totalTime);
      MPMBoxPhases(MPMBoxP box, double parameter, double time, double totalTime);
      MPMBoxPhases(MPMBoxP box, solverType solver, double parameter, double time, double totalTime);
      MPMBoxPhases(MPMBoxP box, solverType solver, updateStressType update, double parameter, double time, double totalTime);

      ~MPMBoxPhases();

      inline void setSolverType(const solverType& type) {d_solver_type = type;}
      inline solverType getSolverType() const {return d_solver_type;}

      inline void setUpdateStressType(const updateStressType& type) {d_update_stress_type = type;}
      inline updateStressType getUpdateStressType() const {return d_update_stress_type;}

      inline void setDeltaT(const double& time) {d_del_t = time;}
      inline double getDeltaT() const {return d_del_t;}

      inline void setTotalT(const double& time) {d_time = time;}
      inline double getTotalT() const {return d_time;}


      inline void setMPMBox(const MPMBoxP& box) {d_mpm_box = box;}
      inline MPMBoxP getMPMBox() const {return d_mpm_box;}

      inline void setNodesConsistentMass(const std::vector<std::vector<double> >& mass) {d_consistent_mass = mass;}
      inline std::vector<std::vector<double> > getNodesConsistentMass() const {return d_consistent_mass;}

      inline void setNodesInverseMass(const std::vector<std::vector<double> >& mass) {d_inverse_mass = mass;}
      inline std::vector<std::vector<double> > getNodesInverseMass() const {return d_inverse_mass;}
 
      
      inline void setNodesLumpedMass(const std::vector<double>& mass) {d_lumped_mass = mass;}
      inline std::vector<double> getNodesLumpedMass() const {return d_lumped_mass;}

      inline void setNodesVelocity(const std::vector<Vector3D>& vel) {d_velocity = vel;}
      inline std::vector<Vector3D> getNodesVelocity() const {return d_velocity;}

      inline void setNodesInternalForce(const std::vector<Vector3D>& force) {d_internal_force = force;}
      inline std::vector<Vector3D> getNodesInternalForce() const {return d_internal_force;}

      inline void setNodesExternalForce(const std::vector<Vector3D>& force) {d_external_force = force;}
      inline std::vector<Vector3D> getNodesExternalForce() const {return d_external_force;}


      inline void setImplicitParameter(const double& param) {d_implicit_parameter = param;}
      inline double getImplicitParameter() const {return d_implicit_parameter;}


      BoxSP getGridBox();
      std::vector<NodeP> getGridNodes();
      std::vector<SymmetricMaterialTrilinearElementSP> getGridElements();
      std::vector<MaterialPointP> getGeometryPoints();

//    ******************************  Material Point Method Initialization Phase  ********************************

      void initializationPhaseLumped();
      void initializationPhaseConsistent();
      void initializationPhaseLumpedUsingMomentum();
      void initializationPhaseConsistentUsingMomentum();

      void createNodesConsistentMass();
      void createNodesInverseMass();
      void createNodesLumpedMass();
      void createNodesVelocityLumped();
      void createNodesVelocityConsistent();
      void createNodesMomentum();
      void createNodesInternalForce();
      void createNodesExternalForce();


      void createNodesIncrementForce();      // MPM Implicit

//    ************************************************************************************************************





//    ********************************  Material Point Method Lagrangian Phase  **********************************

      void lagrangianPhaseLumped();
      void lagrangianPhaseConsistent();
      void lagrangianPhaseLumpedUsingMomentum();
      void lagrangianPhaseConsistentUsingMomentum();

      void findNodesLumpedAcceleration();
      void findNodesConsistentAcceleration();
      void updateNodesVelocity();
      void updateNodesMomentum();
      void updateNodesMomentumImplicit();
      void updatePointsVelocityAndPosition();
      void updatePointsVelocityAndPositionConsistentUsingMomentum();
      void updatePointsVelocityAndPositionLumpedUsingMomentum();
      void updateNodesVelocityLumpedUsingMomentum();
      void updateNodesVelocityConsistentUsingMomentum();

      void findPointsStrainIncrementAndUpdateStrain();
      void findPointsStressIncrementAndUpdateStress();    // constitutive law


      void updateFunctions(MaterialPointP& cur_point, const Vector3D& accSum, const Vector3D& velSum);


      void findNodesVelocityImplicitMethod();     //  MPM Impilict method
      int delta(int i, int j); 
      void updateAcceleration();

      void findNodesTotalForce();       // MPM Implicit Method Using Momentum

//    ************************************************************************************************************





//    *********************************  Material Point Method Convective Phase  *********************************

      void convectivePhase();

      void updatePointsProperties();
      void redefineBackgroundGrid();

      double findMax(const std::vector<double>& vector);
      double findMin(const std::vector<double>& vector);


//    ************************************************************************************************************

      void makeOutputFolder(int& status, int& ret, const std::string& direcName);
      void mpmSolver();
      void mpmApplyBoundaryCondition();
      void mpmApplyBoundaryConditionGeometry();      
      bool nodesBoundaryCondition(NodeP node);
      void mpmApplyPointsForce();
      void mpmApplyPointsForceGeometry();

      void findBoundaryNodesOfPoint(const MaterialPointP& point);
      void findBoundaryNodesOfPoints();

      void updateBoundaryConditionAndPointsForce();

      void mpmPrintOld(int iter);
      void mpmPrintNew();

      void MUMPS(const std::vector<std::vector<double> >& mat, 
                 const std::vector<double>& rhs,
                 std::vector<double>& newDisp);

      void initialMatrix(std::vector <std::vector <double> >& Matrix, int size);
      void initialVector(std::vector <double>& Vector, int size);
      void inverse(const std::vector <std::vector <double> >& mat, std::vector <std::vector <double> >& invMatrix, int size);


      void nodesFile(int& iter_num, const std::string& currentFolder);
      void elementsFile();
      void pointsFile(int& iter_num, const std::string& currentFolder);

   private:

      double d_del_t;
      double d_time;

      MPMBoxP d_mpm_box;
      std::vector<std::vector<double> >  d_consistent_mass;   // for nodes
      std::vector<std::vector<double> >  d_inverse_mass;
      std::vector<double> d_lumped_mass;  // for nodes
      std::vector<Vector3D> d_velocity;  // for nodes
      std::vector<Vector3D> d_internal_force; // for nodes
      std::vector<Vector3D> d_external_force; // for nodes

      solverType d_solver_type;
      updateStressType d_update_stress_type;

      double d_implicit_parameter;

      const double MTOL = 0.001;


  };
}
#endif

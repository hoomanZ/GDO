#ifndef FINITEELEMENT_BODY_H
#define FINITEELEMENT_BODY_H

#include <SymmetricMaterialTrilinearElement.h>
#include <Node.h>
#include <vector>

#include <MaterialPointP.h>  //for MPM


#include <PeriMaterialPointP.h>  //for Peridynamics
#include <StateMaterialPointP.h>
#include <StateMaterialPoint.h>

#include <GeometryMath/Matrix3D.h>

#include <StateBondP.h>



namespace FiniteElement {
  class Body {
    public:
      Body();

      ~Body();


      inline void setDensity(const double& den) {d_density = den;}
      inline double getDensity() const {return d_density;}

      inline void setPoission(const double& poission) {d_poission_ratio = poission;}
      inline double getPoission() const {return d_poission_ratio;}

      inline void setYoung(const double& young) {d_young_modulus = young;}
      inline double getYoung() const {return d_young_modulus;}


      inline void setElements(const std::vector<SymmetricMaterialTrilinearElementSP>& elementsArray) {d_elements_array = elementsArray;}
      inline std::vector<SymmetricMaterialTrilinearElementSP> getElements() const {return d_elements_array;}


      inline void setNodes(const std::vector<NodeP>& nodesArray) {d_nodes_array = nodesArray;}
      inline std::vector<NodeP> getNodes() const {return d_nodes_array;}

      inline void setGlobalStiffness(const std::vector<std::vector<double> >& globalStiffness)
                                                                 {d_global_stiffness = globalStiffness;}
      inline std::vector<std::vector<double> > getGlobalStiffness() const {return d_global_stiffness;}

      inline void setGlobalMass(const std::vector<std::vector<double> >& globalMass) {d_global_mass = globalMass;}
      inline std::vector<std::vector<double> > getGlobalMass() const {return d_global_mass;}

      inline void setExternalForce(const std::vector<double>& globalExternalForce) {d_external_force = globalExternalForce;}
      inline std::vector<double> getExternalForce() const {return d_external_force;}


      inline void setOldGlobDisp(const std::vector<double>& oldDisp) {d_glob_old_disp = oldDisp;}
      inline std::vector<double> getOldGlobDisp() const {return d_glob_old_disp;}

      inline void setGlobDisp(const std::vector<double>& disp) {d_glob_disp = disp;}
      inline std::vector<double> getGlobDisp() const {return d_glob_disp;}

      inline void setNewGlobDisp(const std::vector<double>& newDisp) {d_glob_new_disp = newDisp;}
      inline std::vector<double> getNewGlobDisp() const {return d_glob_new_disp;}

      inline void setGlobalLumpedMass(const std::vector<std::vector<double> >& globalLumpedMass) 
                                                                    {d_global_lumped_mass = globalLumpedMass;}
      inline std::vector<std::vector<double> > getGlobalLumpedMass() const {return d_global_lumped_mass;}


      void clearStiff();
      void clearMass();
      void clearLumpedMass();
      void clearExternalForce();
      void clearOldGlobDisp();
      void clearGlobDisp();
      void clearNewGlobDisp();

      void clearNodesArray();
      void clearElementsArray();


      void initialMatrix(std::vector <std::vector <double> >& Matrix, int size);
      void initialVector(std::vector <double>& Vector, int size);


      void setProperties(double density, double young, double poission);
      void setProperties(double density, double young, double poission, double yield);
      bool isNodeInElement(NodeP node, SymmetricMaterialTrilinearElementSP element);
      int globalIndex(SymmetricMaterialTrilinearElementSP element, int i);      
      void makeGlobalMatrices();




/********************************************************************************************
* 
*  Functions needed for implementing Material Point Method (MPM)
*
*********************************************************************************************/

      
      inline void setMaterialPoints(const std::vector<MaterialPointP>& pointsArray) {d_points_array = pointsArray;}
      inline std::vector<MaterialPointP> getMaterialPoints() const {return d_points_array;}

       
/********************************************************************************************
*
* Functions needed for implementing Bond-based Peridynamics
*
*********************************************************************************************/


      inline void setMaterialPeriPoints(const std::vector<PeriMaterialPointP>& pointsArray)
                                       {d_peri_points_array = pointsArray;}
      inline std::vector<PeriMaterialPointP> getMaterialPeriPoints() const {return d_peri_points_array;}



      inline void setStatePoints(const std::vector<StateMaterialPointP>& pointsArray)
                                       {d_state_points_array = pointsArray;}
      inline std::vector<StateMaterialPointP> getStatePoints() const {return d_state_points_array;}




      void setPeriPoints();


      void createPeriBondsOfEachPeriPoints(const double& horizon);

      void calculateAccelerationOfEachPeriPoints();
      void calculateDisplacementAndVelocityOfPointsUsingVelocityVerlet(const double& del_t);

//********************************************************************************************



/********************************************************************************************
*
* Functions needed for implementing State-based Peridynamics
*
*********************************************************************************************/


      inline void setYieldStrength(const double& yield) {d_yield_strength = yield;}
      inline double getYieldStrength() const {return d_yield_strength;}


      void setStatePoints();
      void setStatePoints(const Matrix3D& initialVelocityGradient);
      void createStateBondsOfEachStatePoints(const double& horizon);
      void createStateBondsOfEachStatePoints(const double& delX, const double& horizon);
      void calculateAccelerationOfEachStatePoints(const double& delT);
      void calculateAccelerationOfEachStatePoints(const double& delT, const double& springConstant);
      StateBondP twinBond(const StateBondP& bond);


      void accelerationOfEachStatePoints(const double& delT);




//********************************************************************************************

    private:

      double d_density;
      double d_poission_ratio;
      double d_young_modulus;

      std::vector<SymmetricMaterialTrilinearElementSP> d_elements_array;
      std::vector<NodeP> d_nodes_array;

      std::vector<std::vector<double> > d_global_stiffness;
      std::vector<std::vector<double> > d_global_mass;
      std::vector<std::vector<double> > d_global_lumped_mass;
 
      std::vector<double> d_external_force;
      std::vector<double> d_glob_old_disp;
      std::vector<double> d_glob_disp;
      std::vector<double> d_glob_new_disp;     
     

/********************************************************************************************
*
* Variables needed for implementing Material Point Method (MPM)
*
*********************************************************************************************/
      
      std::vector<MaterialPointP> d_points_array;
      


/********************************************************************************************
*
* Variables needed for implementing Bond-based Peridynamics
*
*********************************************************************************************/

      std::vector<PeriMaterialPointP> d_peri_points_array;


/********************************************************************************************
*
* Variables needed for implementing State-based Peridynamics
*
*********************************************************************************************/


      std::vector<StateMaterialPointP> d_state_points_array;
      
      double d_yield_strength;


  }; // end of class.
} // end of namespace
#endif


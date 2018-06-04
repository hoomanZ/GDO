#ifndef FINITEELEMENT_SYMMETRIC_MATERIAL_TRILINEAR_ELEMENT_H
#define FINITEELEMENT_SYMMETRIC_MATERIAL_TRILINEAR_ELEMENT_H

#include <Tensor.h>
#include <NodeP.h>
#include <SymmetricMaterial.h>
#include <SymmetricMaterialSP.h>
#include <SymmetricMaterialTrilinearElementSP.h>
#include <TrilinearVolumeElementSP.h>
#include <TrilinearVolumeElement.h>
#include <GeometryMath/Types.h>

#include <vector>

namespace FiniteElement {
  class SymmetricMaterialTrilinearElement : public SymmetricMaterial {
    public:
      SymmetricMaterialTrilinearElement();
      SymmetricMaterialTrilinearElement(const int& id, const double& density, const double& young, const double& poission,
                                        const NodeP& node1, const NodeP& node2, const NodeP& node3, const NodeP& node4,
                                        const NodeP& node5, const NodeP& node6, const NodeP& node7, const NodeP& node8);
      SymmetricMaterialTrilinearElement(const int& id, const double& density, const double& young, const double& poission,
                                        std::vector <NodeP> nodes,
                                        int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
                                        
      SymmetricMaterialTrilinearElement(const int& id, const double& density, const double& young, const double& poission,
                                        const TrilinearVolumeElementSP& triElement); 
      SymmetricMaterialTrilinearElement(const SymmetricMaterialTrilinearElement& mat);
      SymmetricMaterialTrilinearElement(const SymmetricMaterialTrilinearElementSP& mat);
      ~SymmetricMaterialTrilinearElement();

      inline void triElement(const TrilinearVolumeElementSP& t) {d_tri_element = t;}
      inline TrilinearVolumeElementSP triElement() const {return d_tri_element;}

      inline void stiffMat(const Matrix24& mat) {d_stiffness_matrix= mat;}
      inline Matrix24 stiffMat() const {return d_stiffness_matrix;}

      inline void massMat(const Matrix24& mat) {d_mass_matrix= mat;}
      inline Matrix24 massMat() const {return d_mass_matrix;}

      inline void lumpedMassMat(const Matrix24& mat) {d_lumped_mass_matrix = mat;}
      inline Matrix24 lumpedMassMat() const {return d_lumped_mass_matrix;}

      inline void surfaceVec(const Vector24& vec) {d_surface_force_vector = vec;}
      inline Vector24 surfaceVec() const {return d_surface_force_vector;}

      inline void bodyVec(const Vector24& vec) {d_body_force_vector = vec;}
      inline Vector24 bodyVec() const {return d_body_force_vector;}

      inline void oldDispVec(const Vector24& vec) {d_old_disp_vector = vec;}
      inline Vector24 oldDispVec() const {return d_old_disp_vector;}

      inline void dispVec(const Vector24& vec) {d_disp_vector = vec;}
      inline Vector24 dispVec() const {return d_disp_vector;}

      inline void newDispVec(const Vector24& vec) {d_new_disp_vector = vec;}
      inline Vector24 newDispVec() const {return d_new_disp_vector;}



//      inline void extForceVec(const Vector24& forceVec) {d_body_force_vector = vec;}
      inline Vector24 extForceVec() const {return d_body_force_vector + d_surface_force_vector;}


      inline void setRefinedNodes(const std::vector <std::vector <double> >& nodes) {d_refined_nodes = nodes;}
      inline std::vector <std::vector <double> > getRefinedNodes() const {return d_refined_nodes;}




      void operator = (const SymmetricMaterialTrilinearElement& material);
      void operator = (const SymmetricMaterialTrilinearElementSP& material);



      double twoPointsIntegral(double (SymmetricMaterialTrilinearElement::*function)(double xi1, double xi2, double xi3),
                                                                                   SymmetricMaterialTrilinearElement& a);


      double threePointsIntegral(double (SymmetricMaterialTrilinearElement::*function)(double xi1, double xi2, double xi3),
                                                                                   SymmetricMaterialTrilinearElement& a);


      double fourPointsIntegral(double (SymmetricMaterialTrilinearElement::*function)(double xi1, double xi2, double xi3),
                                                                                   SymmetricMaterialTrilinearElement& a);


      double fivePointsIntegral(double (SymmetricMaterialTrilinearElement::*function)(double xi1, double xi2, double xi3),
                                                                                   SymmetricMaterialTrilinearElement& a);


      double sixPointsIntegral(double (SymmetricMaterialTrilinearElement::*function)(double xi1, double xi2, double xi3),
                                                                                   SymmetricMaterialTrilinearElement& a);


      double symStiffnessIntegrand(unsigned int a, unsigned int b,
                                   unsigned int i, unsigned int k,
                                   double xi1, double xi2, double xi3);
      double symMassIntegrand(unsigned int a, unsigned int b,
                              double xi1, double xi2, double xi3);
      
      double symSurfaceIntegrand(unsigned int b, unsigned int i,
                                 unsigned int ni, unsigned int nj,
                                 unsigned int nk, unsigned int nl,
                                 double xi1, double xi2);
                                 
      double symBodyIntegrand(unsigned int b, unsigned int i,
                              double xi1, double xi2, double xi3);

      void stiffnessMatrix(unsigned int number);
      void massMatrix(unsigned int number);
      void surfaceForceVector(unsigned int number);
      void bodyForceVector(unsigned int number);
      

      void  createVectorsOfPointsAndWeights(unsigned int num, std::vector<double>& gaussianVec,
                                                              std::vector<double>& quadratureVec);

      
      double examFunc(double a, double b, double c) {return a*a+2*b*b+3*c*c;}

      void printExamFunc();
      void createDispVectors(); 












    private:
      TrilinearVolumeElementSP d_tri_element;
      
      std::array <double, 2>  d_two_gaussian_points;
      std::array <double, 2>  d_two_quadrature_weights;

      std::array <double, 3>  d_three_gaussian_points;
      std::array <double, 3>  d_three_quadrature_weights;

      std::array <double, 4>  d_four_gaussian_points;
      std::array <double, 4>  d_four_quadrature_weights;

      std::array <double, 5>  d_five_gaussian_points;
      std::array <double, 5>  d_five_quadrature_weights;

      std::array <double, 6>  d_six_gaussian_points;
      std::array <double, 6>  d_six_quadrature_weights;


//      Matrix24 d_stiffness_matrix;
//      Matrix24 d_mass_matrix;
//      Matrix24 d_lumped_mass_matrix;

//      Vector24 d_surface_force_vector;
//      Vector24 d_body_force_vector;
//      Vector24 d_old_disp_vector;
//      Vector24 d_disp_vector;
//      Vector24 d_new_disp_vector;



      Matrix24 d_stiffness_matrix = Matrix24::Zero();
      Matrix24 d_mass_matrix = Matrix24::Zero();
      Matrix24 d_lumped_mass_matrix = Matrix24::Zero();

      Vector24 d_surface_force_vector = Vector24::Zero();
      Vector24 d_body_force_vector = Vector24::Zero();
      Vector24 d_old_disp_vector = Vector24::Zero();
      Vector24 d_disp_vector = Vector24::Zero();
      Vector24 d_new_disp_vector = Vector24::Zero();



      std::vector <std::vector <double> >  d_refined_nodes;  // Making material points using elements

 
  };//end of class
} //end of namespace
#endif

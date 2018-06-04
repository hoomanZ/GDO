#ifndef FINITEELEMENT_SYMMETRIC_MATERIAL_H
#define FINITEELEMENT_SYMMETRIC_MATERIAL_H

#include <Tensor.h>
#include <SymmetricMaterialSP.h>

namespace FiniteElement {
  class SymmetricMaterial {
    public:
      SymmetricMaterial();
      SymmetricMaterial(int id, double density, double young, double poission);
      SymmetricMaterial(const SymmetricMaterial& mat);
      SymmetricMaterial(const SymmetricMaterialSP& mat);
      ~SymmetricMaterial();

      void operator = (const SymmetricMaterial& material);
      void operator = (const SymmetricMaterialSP& material);

      inline void id(const int& ID) {d_id = ID;}
      inline int id() const {return d_id;}

      inline void density(const double& den) {d_density = den;}
      inline double density() const {return d_density;}

      inline void youngModulus(const double& young) {d_young_modulus = young;}
      inline double youngModulus() const {return d_young_modulus;}

      inline void poissionRatio(const double& poission) {d_poission_ratio = poission;}
      inline double poissionRatio() const {return d_poission_ratio;}

      inline void firstLame(const double& lambda) {d_lambda = lambda;}
      inline double firstLame() const {return d_lambda;}

      inline void shearModulus(const double& mu) {d_mu = mu;}
      inline double shearModulus() const {return d_mu;}


      inline void tensor(const Tensor& t) {d_c = t;}
      inline Tensor tensor() const {return d_c;}
   


      int deltaFunction(int i, int j) {return (i != j) ? 0 : 1;}








    private:
      int d_id;
      
      double d_density;
      double d_young_modulus;
      double d_poission_ratio;
      double d_lambda;  //Lame's first parameter
      double d_mu; //shear modulus (Lame's second parameter)

      Tensor d_c;
      





  }; //end of class
} //end of namespace
#endif

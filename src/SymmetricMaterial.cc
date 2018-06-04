#include <SymmetricMaterial.h>

#include <iostream>

using namespace FiniteElement;

SymmetricMaterial::SymmetricMaterial() 
  : d_id(0), d_density(0.0),
    d_young_modulus(0.0), d_poission_ratio(0.0)
{
  d_lambda = (d_young_modulus*d_poission_ratio)/((1+d_poission_ratio)*(1-2*d_poission_ratio));
  d_mu = d_young_modulus/(2*(1+d_poission_ratio));
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          d_c(i, j, k, l)=d_lambda*deltaFunction(i,j)*deltaFunction(k,l)+
                          d_mu*(deltaFunction(i,k)*deltaFunction(j,l)+deltaFunction(i,l)*deltaFunction(j,k));
        }
      }
    }
  }
}

SymmetricMaterial::SymmetricMaterial(int id, double density, double young, double poission)
  : d_id(id), d_density(density),
    d_young_modulus(young), d_poission_ratio(poission)
{
  d_lambda = (d_young_modulus*d_poission_ratio)/((1+d_poission_ratio)*(1-2*d_poission_ratio));
  d_mu = d_young_modulus/(2*(1+d_poission_ratio));
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          d_c(i, j, k, l)=d_lambda*deltaFunction(i,j)*deltaFunction(k,l)+
                          d_mu*(deltaFunction(i,k)*deltaFunction(j,l)+deltaFunction(i,l)*deltaFunction(j,k));
        }
      }
    }
  }
}

SymmetricMaterial::SymmetricMaterial(const SymmetricMaterial& mat) 
  : d_id(mat.d_id), d_density(mat.d_density),
    d_young_modulus(mat.d_young_modulus), d_poission_ratio(mat.d_poission_ratio)
{
  d_lambda = (d_young_modulus*d_poission_ratio)/((1+d_poission_ratio)*(1-2*d_poission_ratio));
  d_mu = d_young_modulus/(2*(1+d_poission_ratio));
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          d_c(i, j, k, l)=d_lambda*deltaFunction(i,j)*deltaFunction(k,l)+
                          d_mu*(deltaFunction(i,k)*deltaFunction(j,l)+deltaFunction(i,l)*deltaFunction(j,k));
        }
      }
    }
  }
}

SymmetricMaterial::SymmetricMaterial(const SymmetricMaterialSP& mat) 
  : d_id(mat->d_id), d_density(mat->d_density),
    d_young_modulus(mat->d_young_modulus), d_poission_ratio(mat->d_poission_ratio)
{
  d_lambda = (d_young_modulus*d_poission_ratio)/((1+d_poission_ratio)*(1-2*d_poission_ratio));
  d_mu = d_young_modulus/(2*(1+d_poission_ratio));
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
        for (int l=0;l<3;l++){
          d_c(i, j, k, l)=d_lambda*deltaFunction(i,j)*deltaFunction(k,l)+
                          d_mu*(deltaFunction(i,k)*deltaFunction(j,l)+deltaFunction(i,l)*deltaFunction(j,k));
        }
      }
    }
  }
}

SymmetricMaterial::~SymmetricMaterial()
{
}

void
SymmetricMaterial::operator = (const SymmetricMaterial& material)
{
  id(material.id());
  density(material.density());
  youngModulus(material.youngModulus());
  poissionRatio(material.poissionRatio());
  firstLame(material.firstLame());
  shearModulus(material.shearModulus());
  tensor(material.tensor());
}

void
SymmetricMaterial::operator = (const SymmetricMaterialSP& material)
{
  id(material->id());
  density(material->density());
  youngModulus(material->youngModulus());
  poissionRatio(material->poissionRatio());
  firstLame(material->firstLame());
  shearModulus(material->shearModulus());
  tensor(material->tensor());
}













   

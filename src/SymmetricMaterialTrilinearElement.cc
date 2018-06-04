#include <SymmetricMaterialTrilinearElement.h>
#include <SymmetricMaterialTrilinearElementSP.h>

#include <iostream>

using namespace FiniteElement;
//using namespace SymmetricMaterial;

SymmetricMaterialTrilinearElement::SymmetricMaterialTrilinearElement() 
  : SymmetricMaterial::SymmetricMaterial(), d_tri_element(new TrilinearVolumeElement())
{
  std::fill_n(d_two_gaussian_points.begin(), 1, -0.577350269189626);
  std::fill_n(d_two_gaussian_points.begin()+1, 1, 0.577350269189626);
  std::fill_n(d_two_quadrature_weights.begin(), 1, 1.0);
  std::fill_n(d_two_quadrature_weights.begin()+1, 1, 1.0);


  std::fill_n(d_three_gaussian_points.begin(), 1, -0.774596669241483);
  std::fill_n(d_three_gaussian_points.begin()+1, 1, 0.0);
  std::fill_n(d_three_gaussian_points.begin()+2, 1, 0.774596669241483);
  std::fill_n(d_three_quadrature_weights.begin(), 1, 0.555555555555556);
  std::fill_n(d_three_quadrature_weights.begin()+1, 1, 0.888888888888889);
  std::fill_n(d_three_quadrature_weights.begin()+2, 1, 0.555555555555556);


  std::fill_n(d_four_gaussian_points.begin(), 1, -0.861136311594053);
  std::fill_n(d_four_gaussian_points.begin()+1, 1, -0.339981043584856);
  std::fill_n(d_four_gaussian_points.begin()+2, 1, 0.339981043584856);
  std::fill_n(d_four_gaussian_points.begin()+3, 1, 0.861136311594053);
  std::fill_n(d_four_quadrature_weights.begin(), 1, 0.347854845137454);
  std::fill_n(d_four_quadrature_weights.begin()+1, 1, 0.652145154862546);
  std::fill_n(d_four_quadrature_weights.begin()+2, 1, 0.652145154862546);
  std::fill_n(d_four_quadrature_weights.begin()+3, 1, 0.347854845137454);


  std::fill_n(d_five_gaussian_points.begin(), 1, -0.906179845938664);
  std::fill_n(d_five_gaussian_points.begin()+1, 1, -0.538469310105683);
  std::fill_n(d_five_gaussian_points.begin()+2, 1, 0.0);
  std::fill_n(d_five_gaussian_points.begin()+3, 1, 0.538469310105683);
  std::fill_n(d_five_gaussian_points.begin()+4, 1, 0.906179845938664);
  std::fill_n(d_five_quadrature_weights.begin(), 1, 0.236926885056189);
  std::fill_n(d_five_quadrature_weights.begin()+1, 1, 0.478628670499366);
  std::fill_n(d_five_quadrature_weights.begin()+2, 1, 0.568888888888889);
  std::fill_n(d_five_quadrature_weights.begin()+3, 1, 0.478628670499366);
  std::fill_n(d_five_quadrature_weights.begin()+4, 1, 0.236926885056189);


  std::fill_n(d_six_gaussian_points.begin(), 1, -0.932469514203152);
  std::fill_n(d_six_gaussian_points.begin()+1, 1, -0.661209386466265);
  std::fill_n(d_six_gaussian_points.begin()+2, 1, -0.238619186083197);
  std::fill_n(d_six_gaussian_points.begin()+3, 1, 0.238619186083197);
  std::fill_n(d_six_gaussian_points.begin()+4, 1, 0.661209386466265);
  std::fill_n(d_six_gaussian_points.begin()+5, 1, 0.932469514203152);
  std::fill_n(d_six_quadrature_weights.begin(), 1, 0.171324492379170);
  std::fill_n(d_six_quadrature_weights.begin()+1, 1, 0.360761573048139);
  std::fill_n(d_six_quadrature_weights.begin()+2, 1, 0.467913934572691);
  std::fill_n(d_six_quadrature_weights.begin()+3, 1, 0.467913934572691);
  std::fill_n(d_six_quadrature_weights.begin()+4, 1, 0.360761573048139);
  std::fill_n(d_six_quadrature_weights.begin()+5, 1, 0.171324492379170);

////  d_stiffness_matrix.zeros();
////  d_mass_matrix.zeros();
////  d_surface_force_vector.zeros();
////  d_body_force_vector.zeros();
////  d_lumped_mass_matrix.zeros();
////  d_old_disp_vector.zeros();
////  d_disp_vector.zeros();
////  d_new_disp_vector.zeros();
}

SymmetricMaterialTrilinearElement::SymmetricMaterialTrilinearElement(const int& id, const double& density, 
                                                                     const double& young, const double& poission,
                                                                     const NodeP& node1, const NodeP& node2,
                                                                     const NodeP& node3, const NodeP& node4,
                                                                     const NodeP& node5, const NodeP& node6,
                                                                     const NodeP& node7, const NodeP& node8)
  : SymmetricMaterial::SymmetricMaterial(id, density, young, poission), 
    d_tri_element(new TrilinearVolumeElement(node1, node2, node3, node4,
                                             node5, node6, node7, node8))
{
  std::fill_n(d_two_gaussian_points.begin(), 1, -0.577350269189626);
  std::fill_n(d_two_gaussian_points.begin()+1, 1, 0.577350269189626);
  std::fill_n(d_two_quadrature_weights.begin(), 1, 1.0);
  std::fill_n(d_two_quadrature_weights.begin()+1, 1, 1.0);


  std::fill_n(d_three_gaussian_points.begin(), 1, -0.774596669241483);
  std::fill_n(d_three_gaussian_points.begin()+1, 1, 0.0);
  std::fill_n(d_three_gaussian_points.begin()+2, 1, 0.774596669241483);
  std::fill_n(d_three_quadrature_weights.begin(), 1, 0.555555555555556);
  std::fill_n(d_three_quadrature_weights.begin()+1, 1, 0.888888888888889);
  std::fill_n(d_three_quadrature_weights.begin()+2, 1, 0.555555555555556);


  std::fill_n(d_four_gaussian_points.begin(), 1, -0.861136311594053);
  std::fill_n(d_four_gaussian_points.begin()+1, 1, -0.339981043584856);
  std::fill_n(d_four_gaussian_points.begin()+2, 1, 0.339981043584856);
  std::fill_n(d_four_gaussian_points.begin()+3, 1, 0.861136311594053);
  std::fill_n(d_four_quadrature_weights.begin(), 1, 0.347854845137454);
  std::fill_n(d_four_quadrature_weights.begin()+1, 1, 0.652145154862546);
  std::fill_n(d_four_quadrature_weights.begin()+2, 1, 0.652145154862546);
  std::fill_n(d_four_quadrature_weights.begin()+3, 1, 0.347854845137454);


  std::fill_n(d_five_gaussian_points.begin(), 1, -0.906179845938664);
  std::fill_n(d_five_gaussian_points.begin()+1, 1, -0.538469310105683);
  std::fill_n(d_five_gaussian_points.begin()+2, 1, 0.0);
  std::fill_n(d_five_gaussian_points.begin()+3, 1, 0.538469310105683);
  std::fill_n(d_five_gaussian_points.begin()+4, 1, 0.906179845938664);
  std::fill_n(d_five_quadrature_weights.begin(), 1, 0.236926885056189);
  std::fill_n(d_five_quadrature_weights.begin()+1, 1, 0.478628670499366);
  std::fill_n(d_five_quadrature_weights.begin()+2, 1, 0.568888888888889);
  std::fill_n(d_five_quadrature_weights.begin()+3, 1, 0.478628670499366);
  std::fill_n(d_five_quadrature_weights.begin()+4, 1, 0.236926885056189);


  std::fill_n(d_six_gaussian_points.begin(), 1, -0.932469514203152);
  std::fill_n(d_six_gaussian_points.begin()+1, 1, -0.661209386466265);
  std::fill_n(d_six_gaussian_points.begin()+2, 1, -0.238619186083197);
  std::fill_n(d_six_gaussian_points.begin()+3, 1, 0.238619186083197);
  std::fill_n(d_six_gaussian_points.begin()+4, 1, 0.661209386466265);
  std::fill_n(d_six_gaussian_points.begin()+5, 1, 0.932469514203152);
  std::fill_n(d_six_quadrature_weights.begin(), 1, 0.171324492379170);
  std::fill_n(d_six_quadrature_weights.begin()+1, 1, 0.360761573048139);
  std::fill_n(d_six_quadrature_weights.begin()+2, 1, 0.467913934572691);
  std::fill_n(d_six_quadrature_weights.begin()+3, 1, 0.467913934572691);
  std::fill_n(d_six_quadrature_weights.begin()+4, 1, 0.360761573048139);
  std::fill_n(d_six_quadrature_weights.begin()+5, 1, 0.171324492379170);

////  d_stiffness_matrix.zeros();
////  d_mass_matrix.zeros();
////  d_surface_force_vector.zeros();
////  d_body_force_vector.zeros();
////  d_lumped_mass_matrix.zeros();
////  d_old_disp_vector.zeros();
////  d_disp_vector.zeros();
////  d_new_disp_vector.zeros();
}


SymmetricMaterialTrilinearElement::SymmetricMaterialTrilinearElement(const int& id, const double& density,
                                                                     const double& young, const double& poission,
                                                                     std::vector <NodeP> nodes,
                                                                     int n1, int n2, int n3, int n4, 
                                                                     int n5, int n6, int n7, int n8)
  : SymmetricMaterial::SymmetricMaterial(id, density, young, poission), 
    d_tri_element(new TrilinearVolumeElement(nodes, n1, n2, n3, n4, n5, n6, n7, n8))
{
  std::fill_n(d_two_gaussian_points.begin(), 1, -0.577350269189626);
  std::fill_n(d_two_gaussian_points.begin()+1, 1, 0.577350269189626);
  std::fill_n(d_two_quadrature_weights.begin(), 1, 1.0);
  std::fill_n(d_two_quadrature_weights.begin()+1, 1, 1.0);


  std::fill_n(d_three_gaussian_points.begin(), 1, -0.774596669241483);
  std::fill_n(d_three_gaussian_points.begin()+1, 1, 0.0);
  std::fill_n(d_three_gaussian_points.begin()+2, 1, 0.774596669241483);
  std::fill_n(d_three_quadrature_weights.begin(), 1, 0.555555555555556);
  std::fill_n(d_three_quadrature_weights.begin()+1, 1, 0.888888888888889);
  std::fill_n(d_three_quadrature_weights.begin()+2, 1, 0.555555555555556);


  std::fill_n(d_four_gaussian_points.begin(), 1, -0.861136311594053);
  std::fill_n(d_four_gaussian_points.begin()+1, 1, -0.339981043584856);
  std::fill_n(d_four_gaussian_points.begin()+2, 1, 0.339981043584856);
  std::fill_n(d_four_gaussian_points.begin()+3, 1, 0.861136311594053);
  std::fill_n(d_four_quadrature_weights.begin(), 1, 0.347854845137454);
  std::fill_n(d_four_quadrature_weights.begin()+1, 1, 0.652145154862546);
  std::fill_n(d_four_quadrature_weights.begin()+2, 1, 0.652145154862546);
  std::fill_n(d_four_quadrature_weights.begin()+3, 1, 0.347854845137454);


  std::fill_n(d_five_gaussian_points.begin(), 1, -0.906179845938664);
  std::fill_n(d_five_gaussian_points.begin()+1, 1, -0.538469310105683);
  std::fill_n(d_five_gaussian_points.begin()+2, 1, 0.0);
  std::fill_n(d_five_gaussian_points.begin()+3, 1, 0.538469310105683);
  std::fill_n(d_five_gaussian_points.begin()+4, 1, 0.906179845938664);
  std::fill_n(d_five_quadrature_weights.begin(), 1, 0.236926885056189);
  std::fill_n(d_five_quadrature_weights.begin()+1, 1, 0.478628670499366);
  std::fill_n(d_five_quadrature_weights.begin()+2, 1, 0.568888888888889);
  std::fill_n(d_five_quadrature_weights.begin()+3, 1, 0.478628670499366);
  std::fill_n(d_five_quadrature_weights.begin()+4, 1, 0.236926885056189);


  std::fill_n(d_six_gaussian_points.begin(), 1, -0.932469514203152);
  std::fill_n(d_six_gaussian_points.begin()+1, 1, -0.661209386466265);
  std::fill_n(d_six_gaussian_points.begin()+2, 1, -0.238619186083197);
  std::fill_n(d_six_gaussian_points.begin()+3, 1, 0.238619186083197);
  std::fill_n(d_six_gaussian_points.begin()+4, 1, 0.661209386466265);
  std::fill_n(d_six_gaussian_points.begin()+5, 1, 0.932469514203152);
  std::fill_n(d_six_quadrature_weights.begin(), 1, 0.171324492379170);
  std::fill_n(d_six_quadrature_weights.begin()+1, 1, 0.360761573048139);
  std::fill_n(d_six_quadrature_weights.begin()+2, 1, 0.467913934572691);
  std::fill_n(d_six_quadrature_weights.begin()+3, 1, 0.467913934572691);
  std::fill_n(d_six_quadrature_weights.begin()+4, 1, 0.360761573048139);
  std::fill_n(d_six_quadrature_weights.begin()+5, 1, 0.171324492379170);

////  d_stiffness_matrix.zeros();
////  d_mass_matrix.zeros();
////  d_surface_force_vector.zeros();
////  d_body_force_vector.zeros();
////  d_lumped_mass_matrix.zeros();
////  d_old_disp_vector.zeros();
////  d_disp_vector.zeros();
////  d_new_disp_vector.zeros();
}
                                             


SymmetricMaterialTrilinearElement::SymmetricMaterialTrilinearElement(const int& id, const double& density, 
                                                                     const double& young, const double& poission,
                                                                     const TrilinearVolumeElementSP& triElement)
  : SymmetricMaterial::SymmetricMaterial(id, density, young, poission)
{
  d_tri_element = triElement;
  std::fill_n(d_two_gaussian_points.begin(), 1, -0.577350269189626);
  std::fill_n(d_two_gaussian_points.begin()+1, 1, 0.577350269189626);
  std::fill_n(d_two_quadrature_weights.begin(), 1, 1.0);
  std::fill_n(d_two_quadrature_weights.begin()+1, 1, 1.0);


  std::fill_n(d_three_gaussian_points.begin(), 1, -0.774596669241483);
  std::fill_n(d_three_gaussian_points.begin()+1, 1, 0.0);
  std::fill_n(d_three_gaussian_points.begin()+2, 1, 0.774596669241483);
  std::fill_n(d_three_quadrature_weights.begin(), 1, 0.555555555555556);
  std::fill_n(d_three_quadrature_weights.begin()+1, 1, 0.888888888888889);
  std::fill_n(d_three_quadrature_weights.begin()+2, 1, 0.555555555555556);


  std::fill_n(d_four_gaussian_points.begin(), 1, -0.861136311594053);
  std::fill_n(d_four_gaussian_points.begin()+1, 1, -0.339981043584856);
  std::fill_n(d_four_gaussian_points.begin()+2, 1, 0.339981043584856);
  std::fill_n(d_four_gaussian_points.begin()+3, 1, 0.861136311594053);
  std::fill_n(d_four_quadrature_weights.begin(), 1, 0.347854845137454);
  std::fill_n(d_four_quadrature_weights.begin()+1, 1, 0.652145154862546);
  std::fill_n(d_four_quadrature_weights.begin()+2, 1, 0.652145154862546);
  std::fill_n(d_four_quadrature_weights.begin()+3, 1, 0.347854845137454);


  std::fill_n(d_five_gaussian_points.begin(), 1, -0.906179845938664);
  std::fill_n(d_five_gaussian_points.begin()+1, 1, -0.538469310105683);
  std::fill_n(d_five_gaussian_points.begin()+2, 1, 0.0);
  std::fill_n(d_five_gaussian_points.begin()+3, 1, 0.538469310105683);
  std::fill_n(d_five_gaussian_points.begin()+4, 1, 0.906179845938664);
  std::fill_n(d_five_quadrature_weights.begin(), 1, 0.236926885056189);
  std::fill_n(d_five_quadrature_weights.begin()+1, 1, 0.478628670499366);
  std::fill_n(d_five_quadrature_weights.begin()+2, 1, 0.568888888888889);
  std::fill_n(d_five_quadrature_weights.begin()+3, 1, 0.478628670499366);
  std::fill_n(d_five_quadrature_weights.begin()+4, 1, 0.236926885056189);


  std::fill_n(d_six_gaussian_points.begin(), 1, -0.932469514203152);
  std::fill_n(d_six_gaussian_points.begin()+1, 1, -0.661209386466265);
  std::fill_n(d_six_gaussian_points.begin()+2, 1, -0.238619186083197);
  std::fill_n(d_six_gaussian_points.begin()+3, 1, 0.238619186083197);
  std::fill_n(d_six_gaussian_points.begin()+4, 1, 0.661209386466265);
  std::fill_n(d_six_gaussian_points.begin()+5, 1, 0.932469514203152);
  std::fill_n(d_six_quadrature_weights.begin(), 1, 0.171324492379170);
  std::fill_n(d_six_quadrature_weights.begin()+1, 1, 0.360761573048139);
  std::fill_n(d_six_quadrature_weights.begin()+2, 1, 0.467913934572691);
  std::fill_n(d_six_quadrature_weights.begin()+3, 1, 0.467913934572691);
  std::fill_n(d_six_quadrature_weights.begin()+4, 1, 0.360761573048139);
  std::fill_n(d_six_quadrature_weights.begin()+5, 1, 0.171324492379170);

////  d_stiffness_matrix.zeros();
////  d_mass_matrix.zeros();
////  d_surface_force_vector.zeros();
////  d_body_force_vector.zeros();
////  d_lumped_mass_matrix.zeros();
////  d_old_disp_vector.zeros();
////  d_disp_vector.zeros();
////  d_new_disp_vector.zeros();
}


SymmetricMaterialTrilinearElement::SymmetricMaterialTrilinearElement(const SymmetricMaterialTrilinearElement& mat)
  :SymmetricMaterial::SymmetricMaterial(mat)
{
  triElement(mat.triElement());
  stiffMat(mat.stiffMat());
  massMat(mat.massMat());
  lumpedMassMat(mat.lumpedMassMat());
  surfaceVec(mat.surfaceVec());
  bodyVec(mat.bodyVec());
  oldDispVec(mat.oldDispVec());
  dispVec(mat.dispVec());
  newDispVec(mat.newDispVec());

  std::fill_n(d_two_gaussian_points.begin(), 1, -0.577350269189626);
  std::fill_n(d_two_gaussian_points.begin()+1, 1, 0.577350269189626);
  std::fill_n(d_two_quadrature_weights.begin(), 1, 1.0);
  std::fill_n(d_two_quadrature_weights.begin()+1, 1, 1.0);


  std::fill_n(d_three_gaussian_points.begin(), 1, -0.774596669241483);
  std::fill_n(d_three_gaussian_points.begin()+1, 1, 0.0);
  std::fill_n(d_three_gaussian_points.begin()+2, 1, 0.774596669241483);
  std::fill_n(d_three_quadrature_weights.begin(), 1, 0.555555555555556);
  std::fill_n(d_three_quadrature_weights.begin()+1, 1, 0.888888888888889);
  std::fill_n(d_three_quadrature_weights.begin()+2, 1, 0.555555555555556);


  std::fill_n(d_four_gaussian_points.begin(), 1, -0.861136311594053);
  std::fill_n(d_four_gaussian_points.begin()+1, 1, -0.339981043584856);
  std::fill_n(d_four_gaussian_points.begin()+2, 1, 0.339981043584856);
  std::fill_n(d_four_gaussian_points.begin()+3, 1, 0.861136311594053);
  std::fill_n(d_four_quadrature_weights.begin(), 1, 0.347854845137454);
  std::fill_n(d_four_quadrature_weights.begin()+1, 1, 0.652145154862546);
  std::fill_n(d_four_quadrature_weights.begin()+2, 1, 0.652145154862546);
  std::fill_n(d_four_quadrature_weights.begin()+3, 1, 0.347854845137454);


  std::fill_n(d_five_gaussian_points.begin(), 1, -0.906179845938664);
  std::fill_n(d_five_gaussian_points.begin()+1, 1, -0.538469310105683);
  std::fill_n(d_five_gaussian_points.begin()+2, 1, 0.0);
  std::fill_n(d_five_gaussian_points.begin()+3, 1, 0.538469310105683);
  std::fill_n(d_five_gaussian_points.begin()+4, 1, 0.906179845938664);
  std::fill_n(d_five_quadrature_weights.begin(), 1, 0.236926885056189);
  std::fill_n(d_five_quadrature_weights.begin()+1, 1, 0.478628670499366);
  std::fill_n(d_five_quadrature_weights.begin()+2, 1, 0.568888888888889);
  std::fill_n(d_five_quadrature_weights.begin()+3, 1, 0.478628670499366);
  std::fill_n(d_five_quadrature_weights.begin()+4, 1, 0.236926885056189);


  std::fill_n(d_six_gaussian_points.begin(), 1, -0.932469514203152);
  std::fill_n(d_six_gaussian_points.begin()+1, 1, -0.661209386466265);
  std::fill_n(d_six_gaussian_points.begin()+2, 1, -0.238619186083197);
  std::fill_n(d_six_gaussian_points.begin()+3, 1, 0.238619186083197);
  std::fill_n(d_six_gaussian_points.begin()+4, 1, 0.661209386466265);
  std::fill_n(d_six_gaussian_points.begin()+5, 1, 0.932469514203152);
  std::fill_n(d_six_quadrature_weights.begin(), 1, 0.171324492379170);
  std::fill_n(d_six_quadrature_weights.begin()+1, 1, 0.360761573048139);
  std::fill_n(d_six_quadrature_weights.begin()+2, 1, 0.467913934572691);
  std::fill_n(d_six_quadrature_weights.begin()+3, 1, 0.467913934572691);
  std::fill_n(d_six_quadrature_weights.begin()+4, 1, 0.360761573048139);
  std::fill_n(d_six_quadrature_weights.begin()+5, 1, 0.171324492379170);
}  


SymmetricMaterialTrilinearElement::SymmetricMaterialTrilinearElement(const SymmetricMaterialTrilinearElementSP& mat)
  :SymmetricMaterial::SymmetricMaterial(mat)
{
  triElement(mat->triElement());
  stiffMat(mat->stiffMat());
  massMat(mat->massMat());
  lumpedMassMat(mat->lumpedMassMat());
  surfaceVec(mat->surfaceVec());
  bodyVec(mat->bodyVec());
  oldDispVec(mat->oldDispVec());
  dispVec(mat->dispVec());
  newDispVec(mat->newDispVec());


  std::fill_n(d_two_gaussian_points.begin(), 1, -0.577350269189626);
  std::fill_n(d_two_gaussian_points.begin()+1, 1, 0.577350269189626);
  std::fill_n(d_two_quadrature_weights.begin(), 1, 1.0);
  std::fill_n(d_two_quadrature_weights.begin()+1, 1, 1.0);


  std::fill_n(d_three_gaussian_points.begin(), 1, -0.774596669241483);
  std::fill_n(d_three_gaussian_points.begin()+1, 1, 0.0);
  std::fill_n(d_three_gaussian_points.begin()+2, 1, 0.774596669241483);
  std::fill_n(d_three_quadrature_weights.begin(), 1, 0.555555555555556);
  std::fill_n(d_three_quadrature_weights.begin()+1, 1, 0.888888888888889);
  std::fill_n(d_three_quadrature_weights.begin()+2, 1, 0.555555555555556);


  std::fill_n(d_four_gaussian_points.begin(), 1, -0.861136311594053);
  std::fill_n(d_four_gaussian_points.begin()+1, 1, -0.339981043584856);
  std::fill_n(d_four_gaussian_points.begin()+2, 1, 0.339981043584856);
  std::fill_n(d_four_gaussian_points.begin()+3, 1, 0.861136311594053);
  std::fill_n(d_four_quadrature_weights.begin(), 1, 0.347854845137454);
  std::fill_n(d_four_quadrature_weights.begin()+1, 1, 0.652145154862546);
  std::fill_n(d_four_quadrature_weights.begin()+2, 1, 0.652145154862546);
  std::fill_n(d_four_quadrature_weights.begin()+3, 1, 0.347854845137454);


  std::fill_n(d_five_gaussian_points.begin(), 1, -0.906179845938664);
  std::fill_n(d_five_gaussian_points.begin()+1, 1, -0.538469310105683);
  std::fill_n(d_five_gaussian_points.begin()+2, 1, 0.0);
  std::fill_n(d_five_gaussian_points.begin()+3, 1, 0.538469310105683);
  std::fill_n(d_five_gaussian_points.begin()+4, 1, 0.906179845938664);
  std::fill_n(d_five_quadrature_weights.begin(), 1, 0.236926885056189);
  std::fill_n(d_five_quadrature_weights.begin()+1, 1, 0.478628670499366);
  std::fill_n(d_five_quadrature_weights.begin()+2, 1, 0.568888888888889);
  std::fill_n(d_five_quadrature_weights.begin()+3, 1, 0.478628670499366);
  std::fill_n(d_five_quadrature_weights.begin()+4, 1, 0.236926885056189);


  std::fill_n(d_six_gaussian_points.begin(), 1, -0.932469514203152);
  std::fill_n(d_six_gaussian_points.begin()+1, 1, -0.661209386466265);
  std::fill_n(d_six_gaussian_points.begin()+2, 1, -0.238619186083197);
  std::fill_n(d_six_gaussian_points.begin()+3, 1, 0.238619186083197);
  std::fill_n(d_six_gaussian_points.begin()+4, 1, 0.661209386466265);
  std::fill_n(d_six_gaussian_points.begin()+5, 1, 0.932469514203152);
  std::fill_n(d_six_quadrature_weights.begin(), 1, 0.171324492379170);
  std::fill_n(d_six_quadrature_weights.begin()+1, 1, 0.360761573048139);
  std::fill_n(d_six_quadrature_weights.begin()+2, 1, 0.467913934572691);
  std::fill_n(d_six_quadrature_weights.begin()+3, 1, 0.467913934572691);
  std::fill_n(d_six_quadrature_weights.begin()+4, 1, 0.360761573048139);
  std::fill_n(d_six_quadrature_weights.begin()+5, 1, 0.171324492379170);
}  


SymmetricMaterialTrilinearElement::~SymmetricMaterialTrilinearElement()
{
}


void 
SymmetricMaterialTrilinearElement::operator = (const SymmetricMaterialTrilinearElement& material)
{
  //SymmetricMaterial::operator = (material);
  this->operator = (material);
  triElement(material.triElement());
  stiffMat(material.stiffMat());
  massMat(material.massMat());
  lumpedMassMat(material.lumpedMassMat());
  surfaceVec(material.surfaceVec());
  bodyVec(material.bodyVec());
  oldDispVec(material.oldDispVec());
  dispVec(material.dispVec());
  newDispVec(material.newDispVec());
}


void 
SymmetricMaterialTrilinearElement::operator = (const SymmetricMaterialTrilinearElementSP& material)
{
  SymmetricMaterial::operator = (material);
  triElement(material->triElement());
  stiffMat(material->stiffMat());
  massMat(material->massMat());
  lumpedMassMat(material->lumpedMassMat());
  surfaceVec(material->surfaceVec());
  bodyVec(material->bodyVec());
  oldDispVec(material->oldDispVec());
  dispVec(material->dispVec());
  newDispVec(material->newDispVec());
}


double 
SymmetricMaterialTrilinearElement::twoPointsIntegral(
                        double (SymmetricMaterialTrilinearElement::*function)(double xi1, double xi2, double xi3),
                                                                            SymmetricMaterialTrilinearElement& a)
{
/*  double kesi1 = -0.5773502691;
  double kesi2 = 0.5773502691;
  return (a.*function)(kesi1, kesi1, kesi1) + (a.*function)(kesi1, kesi1, kesi2)+ 
         (a.*function)(kesi1, kesi2, kesi1) + (a.*function)(kesi1, kesi2, kesi2)+ 
         (a.*function)(kesi2, kesi1, kesi1) + (a.*function)(kesi2, kesi1, kesi2)+ 
         (a.*function)(kesi2, kesi2, kesi1) + (a.*function)(kesi2, kesi2, kesi2); */
 double sum = 0.0;
 double kesi1 = 0.0;
 double kesi2 = 0.0;
 double kesi3 = 0.0;
 double w1 = 0.0;
 double w2 = 0.0;
 double w3 = 0.0;
 int size = d_two_gaussian_points.size();
 for(int ii = 0; ii < size; ii++) {
   kesi1 = d_two_gaussian_points[ii];
   w1 = d_two_quadrature_weights[ii];
   for(int jj = 0; jj < size; jj++) {
      kesi2 = d_two_gaussian_points[jj];
      w2 = d_two_quadrature_weights[jj];
      for(int kk = 0; kk < size; kk++) {
          kesi3 = d_two_gaussian_points[kk];
          w3 = d_two_quadrature_weights[kk];
          sum += w1*w2*w3*(a.*function)(kesi1, kesi2, kesi3);
      }
   }
 }
return sum;
}

         
double
SymmetricMaterialTrilinearElement::threePointsIntegral(
                        double (SymmetricMaterialTrilinearElement::*function)(double xi1, double xi2, double xi3),
                                                                            SymmetricMaterialTrilinearElement& a)
{
 double sum = 0.0;
 double kesi1 = 0.0;
 double kesi2 = 0.0;
 double kesi3 = 0.0;
 double w1 = 0.0;
 double w2 = 0.0;
 double w3 = 0.0;
 int size = d_three_gaussian_points.size();
 for(int ii = 0; ii < size; ii++) {
   kesi1 = d_three_gaussian_points[ii];
   w1 = d_three_quadrature_weights[ii];
   for(int jj = 0; jj < size; jj++) {
      kesi2 = d_three_gaussian_points[jj];
      w2 = d_three_quadrature_weights[jj];
      for(int kk = 0; kk < size; kk++) {
          kesi3 = d_three_gaussian_points[kk];
          w3 = d_three_quadrature_weights[kk];
          sum += w1*w2*w3*(a.*function)(kesi1, kesi2, kesi3);
      }
   }
 }
return sum;
}


double
SymmetricMaterialTrilinearElement::fourPointsIntegral(
                        double (SymmetricMaterialTrilinearElement::*function)(double xi1, double xi2, double xi3),
                                                                            SymmetricMaterialTrilinearElement& a)
{
 double sum = 0.0;
 double kesi1 = 0.0;
 double kesi2 = 0.0;
 double kesi3 = 0.0;
 double w1 = 0.0;
 double w2 = 0.0;
 double w3 = 0.0;
 int size = d_four_gaussian_points.size();
 for(int ii = 0; ii < size; ii++) {
   kesi1 = d_four_gaussian_points[ii];
   w1 = d_four_quadrature_weights[ii];
   for(int jj = 0; jj < size; jj++) {
      kesi2 = d_four_gaussian_points[jj];
      w2 = d_four_quadrature_weights[jj];
      for(int kk = 0; kk < size; kk++) {
          kesi3 = d_four_gaussian_points[kk];
          w3 = d_four_quadrature_weights[kk];
          sum += w1*w2*w3*(a.*function)(kesi1, kesi2, kesi3);
      }
   }
 }
return sum;
}


double
SymmetricMaterialTrilinearElement::fivePointsIntegral(
                        double (SymmetricMaterialTrilinearElement::*function)(double xi1, double xi2, double xi3),
                                                                            SymmetricMaterialTrilinearElement& a)
{
 double sum = 0.0;
 double kesi1 = 0.0;
 double kesi2 = 0.0;
 double kesi3 = 0.0;
 double w1 = 0.0;
 double w2 = 0.0;
 double w3 = 0.0;
 int size = d_five_gaussian_points.size();
 for(int ii = 0; ii < size; ii++) {
   kesi1 = d_five_gaussian_points[ii];
   w1 = d_five_quadrature_weights[ii];
   for(int jj = 0; jj < size; jj++) {
      kesi2 = d_five_gaussian_points[jj];
      w2 = d_five_quadrature_weights[jj];
      for(int kk = 0; kk < size; kk++) {
          kesi3 = d_five_gaussian_points[kk];
          w3 = d_five_quadrature_weights[kk];
          sum += w1*w2*w3*(a.*function)(kesi1, kesi2, kesi3);
      }
   }
 }
return sum;
}


double
SymmetricMaterialTrilinearElement::sixPointsIntegral(
                        double (SymmetricMaterialTrilinearElement::*function)(double xi1, double xi2, double xi3),
                                                                            SymmetricMaterialTrilinearElement& a)
{
 double sum = 0.0;
 double kesi1 = 0.0;
 double kesi2 = 0.0;
 double kesi3 = 0.0;
 double w1 = 0.0;
 double w2 = 0.0;
 double w3 = 0.0;
 int size = d_six_gaussian_points.size();
 for(int ii = 0; ii < size; ii++) {
   kesi1 = d_six_gaussian_points[ii];
   w1 = d_six_quadrature_weights[ii];
   for(int jj = 0; jj < size; jj++) {
      kesi2 = d_six_gaussian_points[jj];
      w2 = d_six_quadrature_weights[jj];
      for(int kk = 0; kk < size; kk++) {
          kesi3 = d_six_gaussian_points[kk];
          w3 = d_six_quadrature_weights[kk];
          sum += w1*w2*w3*(a.*function)(kesi1, kesi2, kesi3);
      }
   }
 }
return sum;
}


void
SymmetricMaterialTrilinearElement::printExamFunc() {
 SymmetricMaterialTrilinearElement aa;
  std::cout << aa.sixPointsIntegral(&SymmetricMaterialTrilinearElement::examFunc, aa) << std::endl;
}


double
SymmetricMaterialTrilinearElement::symStiffnessIntegrand(unsigned int a, unsigned int b,
                                             unsigned int i, unsigned int k,
                                             double xi1, double xi2, double xi3) 
{
  double sum = 0.0;
  for (unsigned int j = 1; j < 4; j++) {
    for (unsigned int l = 1; l < 4; l++) {
      sum += this->tensor()(i-1, j-1, k-1, l-1)
                               *d_tri_element->symDerivativeTestFunctionToPosition(a, l, xi1, xi2, xi3)
                               *d_tri_element->symDerivativeTestFunctionToPosition(b, j, xi1, xi2, xi3);
    }
  }
  sum *= d_tri_element->symJacobianMatrix(xi1, xi2, xi3).Determinant();
  return sum;
}

double
SymmetricMaterialTrilinearElement::symMassIntegrand(unsigned int a, unsigned int b,
                                                    double xi1, double xi2, double xi3) 
{
  double sum = 0.0;
//  for (unsigned int j = 1; j < 4; j++) {
//    for (unsigned int l = 1; l < 4; l++) {
//      std::cout << d_tri_element->symInterpolatedValue(d_tri_element->densityArray(), xi1, xi2, xi3) << std::endl;
//      sum += d_tri_element->symInterpolatedValue(d_tri_element->densityArray(), xi1, xi2, xi3)
      sum = d_tri_element->symInterpolatedValue(d_tri_element->densityArray(), xi1, xi2, xi3)
                               *d_tri_element->symTestFunction(a, xi1, xi2, xi3)
                               *d_tri_element->symTestFunction(b, xi1, xi2, xi3);
//    }
//  }
  sum *= d_tri_element->symJacobianMatrix(xi1, xi2, xi3).Determinant();
  return sum;
}



double
SymmetricMaterialTrilinearElement::symSurfaceIntegrand(unsigned int b, unsigned int i,
                                                       unsigned int ni, unsigned int nj,
                                                       unsigned int nk, unsigned int nl,
                                                       double xi1, double xi2)
{
  double jacob = 0.0;
  double interpolatedForce = 0.0;     
  std::array<Vector3D, 8>  arrayOfForces = d_tri_element->surfaceForceArray(ni, nj, nk, nl);
  if ((ni == 5) && (nj == 6) && (nk == 7) && (nl == 8)) {
//    std::cout << "Xi3Max" << std::endl;
    jacob = d_tri_element->symNormVecXi3Max(xi1, xi2).length();
    interpolatedForce = d_tri_element->symInterpolatedVector(arrayOfForces, xi1, xi2, 1)[i-1];
    return interpolatedForce*(d_tri_element->symTestFunction(b, xi1, xi2, 1))*jacob;
  } else if ((ni == 1) && (nj == 2 ) && (nk == 3) && (nl == 4)) {
//    std::cout << "Xi3Min" << std::endl;
    jacob = d_tri_element->symNormVecXi3Min(xi1, xi2).length();
    interpolatedForce = d_tri_element->symInterpolatedVector(arrayOfForces, xi1, xi2, -1)[i-1];
    return interpolatedForce*(d_tri_element->symTestFunction(b, xi1, xi2, -1))*jacob;
  } else if ((ni == 3) && (nj == 4) && (nk == 5) && (nl == 8)) {
//    std::cout << "Xi2Max" << std::endl;
    jacob = d_tri_element->symNormVecXi2Max(xi1, xi2).length();
    interpolatedForce = d_tri_element->symInterpolatedVector(arrayOfForces, xi1, 1, xi2)[i-1];
    return interpolatedForce*(d_tri_element->symTestFunction(b, xi1, 1, xi2))*jacob;
  } else if ((ni == 1) && (nj == 2) && (nk == 6) && (nl == 7)) {
//    std::cout << "Xi2Min" << std::endl;
    jacob = d_tri_element->symNormVecXi2Min(xi1, xi2).length();
    interpolatedForce = d_tri_element->symInterpolatedVector(arrayOfForces, xi1, -1, xi2)[i-1];
    return interpolatedForce*(d_tri_element->symTestFunction(b, xi1, -1, xi2))*jacob;
  } else if ((ni == 2) && (nj == 3) && (nk == 7) && (nl == 8)) {
//    std::cout << "Xi1Max" << std::endl;
    jacob = d_tri_element->symNormVecXi1Max(xi1, xi2).length();
//      std::cout << "(" << b << "," << i << "," << xi1 << "," << xi2 << ")  ";
    interpolatedForce = d_tri_element->symInterpolatedVector(arrayOfForces, 1, xi1, xi2)[i-1];
//      std::cout << "interpolatedForce= " << interpolatedForce << "  ";
//      std::cout << "TestFunction= " << d_tri_element->symTestFunction(b, 1, xi1, xi2) << std::endl;   
    return interpolatedForce*(d_tri_element->symTestFunction(b, 1, xi1, xi2))*jacob;
  } else if ((ni == 1) && (nj == 4) && (nk == 5) && (nl == 6)) {
//    std::cout << "Xi1Min" << std::endl;
    jacob = d_tri_element->symNormVecXi1Max(xi1, xi2).length();
    interpolatedForce = d_tri_element->symInterpolatedVector(arrayOfForces, -1, xi1, xi2)[i-1];
    return interpolatedForce*(d_tri_element->symTestFunction(b, -1, xi1, xi2))*jacob; 
  } else {
           std::cout << " The indices i, j, k, and l present the numbers of local nodes in the element." << std::endl;
           std::cout << " The possible amounts of i, j, k, and l are:" << std::endl;
           std::cout << " (i, j, k, l) = (5, 6, 7, 8) for plane xi3 = 1" << std::endl;
           std::cout << " (i, j, k, l) = (1, 2, 3, 4) for plane xi3 = -1" << std::endl;
           std::cout << " (i, j, k, l) = (3, 4, 5, 8) for plane xi2 = 1" << std::endl;
           std::cout << " (i, j, k, l) = (1, 2, 6, 7) for plane xi2 = -1" << std::endl;
           std::cout << " (i, j, k, l) = (2, 3, 7, 8) for plane xi1 = 1" << std::endl;
           std::cout << " (i, j, k, l) = (1, 4, 5, 6) for plane xi1 = -1" << std::endl;
           return -1;
  }
}                                                                    


double
SymmetricMaterialTrilinearElement::symBodyIntegrand(unsigned int b, unsigned int i,
                                                       double xi1, double xi2, double xi3) 
{
    std::array<Vector3D, 8>  arrayOfBodyForces = d_tri_element->bodyForceArray();
   double jacob = 0.0;
   double interpolatedForce = 0.0;
   interpolatedForce = d_tri_element->symInterpolatedVector(arrayOfBodyForces, xi1, xi2, xi3)[i-1];
   jacob = d_tri_element->symJacobianMatrix(xi1, xi2, xi3).Determinant();
   return interpolatedForce*d_tri_element->symTestFunction(b, xi1, xi2, xi3)*jacob;
}    




void
SymmetricMaterialTrilinearElement::createVectorsOfPointsAndWeights(unsigned int num, 
                                                                   std::vector<double>& gaussianVec, 
                                                                   std::vector<double>& quadratureVec)
{
   gaussianVec.clear();
   quadratureVec.clear();
   switch (num)
   {
   case 2:
     {
       for (unsigned int ii = 0; ii < d_two_gaussian_points.size(); ii++) {
         gaussianVec.push_back(d_two_gaussian_points[ii]);
         quadratureVec.push_back(d_two_quadrature_weights[ii]);
       }
       break;
     }
     case 3:
     {
       for (unsigned int ii = 0; ii < d_three_gaussian_points.size(); ii++) {
         gaussianVec.push_back(d_three_gaussian_points[ii]);
         quadratureVec.push_back(d_three_quadrature_weights[ii]);
       }
         break;
     }
     case 4:
     {
       for (unsigned int ii = 0; ii < d_four_gaussian_points.size(); ii++) {
         gaussianVec.push_back(d_four_gaussian_points[ii]);
         quadratureVec.push_back(d_four_quadrature_weights[ii]);
       }
       break;
     }
     case 5:
     {
       for (unsigned int ii = 0; ii < d_five_gaussian_points.size(); ii++) {
         gaussianVec.push_back(d_five_gaussian_points[ii]);
         quadratureVec.push_back(d_five_quadrature_weights[ii]);
       }
       break;
     }
     case 6:
     {
       for (unsigned int ii = 0; ii < d_six_gaussian_points.size(); ii++) {
         gaussianVec.push_back(d_six_gaussian_points[ii]);
         quadratureVec.push_back(d_six_quadrature_weights[ii]);
       }
       break;
     }
     default:
     {
       std::cout << " The index num shows the number of gaussian points which is from 2 to 6.";
       break;
     }
   } //end of switch
} //end of function


void
SymmetricMaterialTrilinearElement::stiffnessMatrix(unsigned int number)
{ 
  double kesi1 = 0.0;
  double kesi2 = 0.0;
  double kesi3 = 0.0;
  double w1 = 0.0;
  double w2 = 0.0;
  double w3 = 0.0;
  std::vector<double> gaussianVector;
  std::vector<double> quadratureVector;
  createVectorsOfPointsAndWeights(number, gaussianVector, quadratureVector);
  double vectorSize = gaussianVector.size(); 
  double integral = 0.0;
  unsigned int index1 = 0;
  unsigned int index2 = 0;
  for (unsigned int a = 1; a < 9; a++) {
    for (unsigned int b = 1; b < 9; b++) {
      for (unsigned int i = 1; i < 4; i++) {
        for (unsigned int k = 1; k < 4; k++) {
          integral = 0.0;
          index1 = 0;
          index2 = 0;
          for (unsigned int numVec1 = 0; numVec1 < vectorSize; numVec1++) {
            kesi1 = gaussianVector[numVec1];
            w1 = quadratureVector[numVec1];
            for (unsigned int numVec2 = 0; numVec2 < vectorSize; numVec2++) {
              kesi2 = gaussianVector[numVec2];
              w2 = quadratureVector[numVec2];
              for (unsigned int numVec3 = 0; numVec3 < vectorSize; numVec3++) {
                kesi3 = gaussianVector[numVec3];
                w3 = quadratureVector[numVec3];
                integral += w1*w2*w3*symStiffnessIntegrand(a, b, i, k, kesi1, kesi2, kesi3);
              }  //numVec3
            }  //numVec2
          }  //numVec1
          index1 = 3*(b-1)+i;
          index2 = 3*(a-1)+k;
          d_stiffness_matrix(index1-1, index2-1) = integral;
        }  //k
      }  //i
    }  //b
  }  //a
   
//  if (id() == 45)
//  {
//    for (int ii = 0; ii < 3; ii++)
//      for (int jj = 0; jj < 3; jj++)
//        for (int kk = 0; kk < 3; kk++)
//          for (int ll = 0; ll < 3; ll++)
//            std::cout << "(" << ii << "," << jj << ","
//                             << kk << "," << ll << ")= " << tensor()(ii, jj, kk, ll) << std::endl;


//    for (int row = 0; row < 24; row++)
//    {
//      std::cout << std::endl;
//      for (int column = 0; column < 24; column++)
//      {
//        std::cout << d_stiffness_matrix(row, column) << " "; 
//      }
//    }  
//  }

  
}  //end of function



void
SymmetricMaterialTrilinearElement::massMatrix(unsigned int number)
{ 
  double kesi1 = 0.0;
  double kesi2 = 0.0;
  double kesi3 = 0.0;
  double w1 = 0.0;
  double w2 = 0.0;
  double w3 = 0.0;
  std::vector<double> gaussianVector;
  std::vector<double> quadratureVector;
  createVectorsOfPointsAndWeights(number, gaussianVector, quadratureVector);
  double vectorSize = gaussianVector.size(); 
  double integral = 0.0;
  unsigned int index1 = 0;
  unsigned int index2 = 0;
  for (unsigned int a = 1; a < 9; a++) {
    for (unsigned int b = 1; b < 9; b++) {
          integral = 0.0;
          index1 = 0;
          index2 = 0;
          for (unsigned int numVec1 = 0; numVec1 < vectorSize; numVec1++) {
            kesi1 = gaussianVector[numVec1];
            w1 = quadratureVector[numVec1];
            for (unsigned int numVec2 = 0; numVec2 < vectorSize; numVec2++) {
              kesi2 = gaussianVector[numVec2];
              w2 = quadratureVector[numVec2];
              for (unsigned int numVec3 = 0; numVec3 < vectorSize; numVec3++) {
                kesi3 = gaussianVector[numVec3];
                w3 = quadratureVector[numVec3];
//                if ((a ==1) && (b == 1))
//                   std::cout << std::endl << symMassIntegrand(a, b, kesi1, kesi2, kesi3) << std::endl;
                integral += w1*w2*w3*symMassIntegrand(a, b, kesi1, kesi2, kesi3);
              }  //numVec3
            }  //numVec2
          }  //numVec1
          index1 = 3*b-2;
          index2 = 3*a-2;
          d_mass_matrix(index1-1, index2-1) = integral;
          d_mass_matrix(index1, index2) = integral;
          d_mass_matrix(index1+1, index2+1) = integral;
    }  //b
  }  //a


  if (id() == 45)
  {
//    for (int ii = 0; ii < 3; ii++)
//      for (int jj = 0; jj < 3; jj++)
//        for (int kk = 0; kk < 3; kk++)
//          for (int ll = 0; ll < 3; ll++)
//            std::cout << "(" << ii << "," << jj << ","
//                             << kk << "," << ll << ")= " << tensor()(ii, jj, kk, ll) << std::endl;


    for (int row = 0; row < 24; row++)
    {
      std::cout << std::endl;
      for (int column = 0; column < 24; column++)
      {
        std::cout << d_mass_matrix(row, column) << " "; 
      }
    }  
  }

  
  //  calculate lumped mass matrix
  for (unsigned int ii = 0; ii < 24; ii++) {
    for (unsigned int jj = 0; jj < 24; jj++) { 
      d_lumped_mass_matrix(ii, ii) += d_mass_matrix(ii, jj);
    }
  } 
}  //end of function



void
SymmetricMaterialTrilinearElement::surfaceForceVector(unsigned int number)
{
  double kesi1 = 0.0;
  double kesi2 = 0.0;
  double w1 = 0.0;
  double w2 = 0.0;
  std::vector<double> gaussianVector;
  std::vector<double> quadratureVector;
  createVectorsOfPointsAndWeights(number, gaussianVector, quadratureVector);
  double vectorSize = gaussianVector.size(); 
  double integral = 0.0;
  unsigned int index = 0;
  int n1 = 1; int n2 = 1;
  int n3 = 1; int n4 = 1;
  int n5 = 1; int n6 = 1;
  for (unsigned int a = 1; a < 9; a++) {
    for (unsigned int i = 1; i < 4; i++) {
          integral = 0.0;
          index = 0;
          for (unsigned int numVec1 = 0; numVec1 < vectorSize; numVec1++) {
            kesi1 = gaussianVector[numVec1];
            w1 = quadratureVector[numVec1];
            for (unsigned int numVec2 = 0; numVec2 < vectorSize; numVec2++) {
              kesi2 = gaussianVector[numVec2];
              w2 = quadratureVector[numVec2];
              if (d_tri_element->face(5, 6, 7, 8)) { // std::cout << id() << " xi3= 1" << std::endl;
                 integral += w1*w2*symSurfaceIntegrand(a, i, 5, 6, 7, 8, kesi1, kesi2);
//                 integral = w1*w2*symSurfaceIntegrand(a, i, 5, 6, 7, 8, kesi1, kesi2);
                if (/*(this->id() == 65)) && */(numVec1 == vectorSize - 1) && (numVec2 == vectorSize - 1))
                {
                 std::cout << "ElementID= " << id()  
                 /*std::cout*/ << " NodeP" << a << 
                 " Number n1= " << n1 << id() << " xi3= 1 " << "integral" << "(" << i << ")= " << integral
                 << std::endl;
                }
                 n1 = n1 + 1;
              }
              if (d_tri_element->face(1, 2, 3, 4)) { // std::cout << id() << " xi3= -1" << std::endl;
                 integral += w1*w2*symSurfaceIntegrand(a, i, 1, 2, 3, 4, kesi1, kesi2);
                 std::cout << "Number n2= " << n2 << " " << id() << " xi3= -1 " << "integral= " << integral << std::endl;
                 n2 = n2 + 1;
              } 
              if (d_tri_element->face(3, 4, 5, 8)) {  std::cout << id() << " xi2= 1" << std::endl;
                 integral += w1*w2*symSurfaceIntegrand(a, i, 3, 4, 5, 8, kesi1, kesi2);
                 std::cout << "Number n3= " << n3 << " " << id() << " xi2= 1 " << "integral= " << integral << std::endl;
                 n3 = n3 + 1;
              } 
              if (d_tri_element->face(1, 2, 6, 7)) {  std::cout << id() << " xi2= -1" << std::endl;
                 integral += w1*w2*symSurfaceIntegrand(a, i, 1, 2, 6, 7, kesi1, kesi2);
                 std::cout << "Number n4= " << n4 << " " << id() << " xi2= -1 " << "integral= " << integral << std::endl;
                 n4 = n4 + 1;
              } 
              if (d_tri_element->face(2, 3, 7, 8)) {  std::cout << id() << " xi1= 1" << std::endl;
                 integral += w1*w2*symSurfaceIntegrand(a, i, 2, 3, 7, 8, kesi1, kesi2);
                 std::cout << "Number n5= " << n5 << " " << id() << " xi1= 1 " << "integral= " << integral << std::endl;
                 n5 = n5 + 1;
              } 
              if (d_tri_element->face(1, 4, 5, 6)) {  std::cout << id() << " xi1= -1" << std::endl;
                 integral += w1*w2*symSurfaceIntegrand(a, i, 1, 4, 5, 6, kesi1, kesi2);
                 std::cout << "Number n6= " << n6 << " " << id() << " xi1= -1 " << "integral= " << integral << std::endl;
                 n6 = n6 + 1;
              }       
            }  //numVec2
          }  //numVec1
          index = 3*a+i-4;
          d_surface_force_vector(index) = integral;
    }  //i
  }  //a 
  
//  int ID = 45;
//  if (id() == ID)
//  {
//    int size = d_surface_force_vector.size();
//    for (int k = 0; k < size; k++)
//    {
//      std::cout << "s(" << k << ")= " << d_surface_force_vector[k] << " ";
//    }
//  }
  
} //end of function


void
SymmetricMaterialTrilinearElement::bodyForceVector(unsigned int number)
{
  double kesi1 = 0.0;
  double kesi2 = 0.0;
  double kesi3 = 0.0;
  double w1 = 0.0;
  double w2 = 0.0;
  double w3 = 0.0;
  std::vector<double> gaussianVector;
  std::vector<double> quadratureVector;
  createVectorsOfPointsAndWeights(number, gaussianVector, quadratureVector);
  double vectorSize = gaussianVector.size(); 
  double integral = 0.0;
  unsigned int index = 0;
  for (unsigned int a = 1; a < 9; a++) {
    for (unsigned int i = 1; i < 4; i++) {
          integral = 0.0;
          index = 0;
          for (unsigned int numVec1 = 0; numVec1 < vectorSize; numVec1++) {
            kesi1 = gaussianVector[numVec1];
            w1 = quadratureVector[numVec1];
            for (unsigned int numVec2 = 0; numVec2 < vectorSize; numVec2++) {
              kesi2 = gaussianVector[numVec2];
              w2 = quadratureVector[numVec2];
              for (unsigned int numVec3 = 0; numVec3 < vectorSize; numVec3++) {
                kesi3 = gaussianVector[numVec3];
                w3 = quadratureVector[numVec3];
                integral += w1*w2*w3*symBodyIntegrand(a, i, kesi1, kesi2, kesi3);
              }  //numVec3
            }  //numVec2
          }  //numVec1
          index = 3*a+i-4;
          d_body_force_vector(index) = integral;
    }  //i
  }  //a 
} //end of function
               

void
SymmetricMaterialTrilinearElement::createDispVectors()
{
  for (int i = 0; i < 8; i++)
  {
    NodeP cur_node = d_tri_element->getArray()[i];
    for (int j = 0; j < 3; j++)
    {
      d_old_disp_vector[3*i+j] = cur_node->oldDisplacement()[j];
      d_disp_vector[3*i+j] = cur_node->displacement()[j];
      d_new_disp_vector[3*i+j] = cur_node->newDisplacement()[j];
    }
  }
}




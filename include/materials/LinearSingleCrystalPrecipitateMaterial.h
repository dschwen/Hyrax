/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  12 January 2012
*
*************************************************************************/


#ifndef LINEARSINGLECRYSTALPRECIPITATEMATERIAL_H
#define LINEARSINGLECRYSTALPRECIPITATEMATERIAL_H

#include "TensorElasticMaterial.h"
#include "ElasticityTensorR4.h"
#include "RankTwoTensorTonks.h"

/**
 * LinearSingleCrystalPrecipitateMaterial handles anisotropic, single-crystal material elastic
 * constants. It handles a single
 * crystal of matrix with an arbitrary number of orientation variants of a coherent precipitate.
 */

//Forward declaration
class LinearSingleCrystalPrecipitateMaterial;

template<>
InputParameters validParams<LinearSingleCrystalPrecipitateMaterial>();

class LinearSingleCrystalPrecipitateMaterial : public TensorElasticMaterial
{
public:
  LinearSingleCrystalPrecipitateMaterial(const std:: string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();

  /**
   * Computes the overall material elasticity tensor for all space.
   */
  virtual void computeQpElasticityTensor();

  /**
   * Computes the overall material eigenstrain tensor for all space.
   */
  virtual void computeQpEigenstrain();

  /**
   * Computes the elastic strain for all space.
   */
  virtual void computeQpElasticStrain();

  /**
   * Computes the elastic stress for all space.
   */
  virtual void computeQpElasticStress();

private:
  // vectors to get the input values
  //std::vector<Real> _Cijkl_matrix_vector;
  // std::vector<Real> _Cijkl_precipitate_vector;
  std::vector<Real> _eigenstrain_vector;

  // number of orientation variants for the precipitate in a single matrix crystal
  unsigned int _n_variants;

  // Indicate if 21 (true) or 9 (false) values given for input of elastic stiffness
  //bool _all_21;

  // Individual material information
  //SymmAnisotropicElasticityTensor _Cijkl_matrix;
  //ElasticityTensorR4 _Cijkl_precipitate;
  RankTwoTensorTonks _eigenstrain;

  //std::vector<ElasticityTensorR4> _Cijkl_precipitates_rotated;
  std::vector<RankTwoTensorTonks> _eigenstrains_rotated;

  MaterialProperty<RankTwoTensorTonks> & _local_strain;
  MaterialProperty<RankTwoTensorTonks> & _misfit_strain;
  MaterialProperty<std::vector<RankTwoTensorTonks> > & _eigenstrains_rotated_MP;
  MaterialProperty<ElasticityTensorR4> & _Cijkl_matrix_MP;
  //MaterialProperty<std::vector<ElasticityTensorR4> > & _Cijkl_precipitates_rotated_MP;

  // derivatives of the local elasticity tensor and the misfit strain with respect
  // to order parameter
  //MaterialProperty<std::vector<ElasticityTensorR4> > & _d_elasticity_tensor;
  MaterialProperty<std::vector<RankTwoTensorTonks> > & _d_eigenstrains_rotated_MP;

  // Vector of references to the coupled order parameters
  std::vector<VariableValue *> _coupled_variables;
};

#endif //LINEARSINGLECRYSTALPRECIPITATEMATERIAL_H

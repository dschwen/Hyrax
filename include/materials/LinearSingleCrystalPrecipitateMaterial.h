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

//#include "LinearElasticMaterial.h"
#include "TensorMechanicsMaterial.h"
#include "ElasticityTensorR4.h"
#include "RankTwoTensor.h"

/**
 * LinearSingleCrystalPrecipitateMaterial handles anisotropic, single-crystal
 * material elastic constants. It handles a single crystal of matrix with an
 * arbitrary number of orientation variants of a coherent precipitate.
 */

//Forward declaration
class LinearSingleCrystalPrecipitateMaterial;

template<>
InputParameters validParams<LinearSingleCrystalPrecipitateMaterial>();

class LinearSingleCrystalPrecipitateMaterial : public TensorMechanicsMaterial //LinearElasticMaterial
{
public:
  LinearSingleCrystalPrecipitateMaterial(const std:: string & name, InputParameters parameters);

protected:
  virtual void computeProperties();

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

  virtual void computeQpStrain();

  virtual void computeQpStress();

private:
  // vectors to get the input values
  std::vector<Real> _Cijkl_precipitate_vector;
  std::vector<Real> _eigenstrain_vector;

  // Individual material information
  ElasticityTensorR4 _Cijkl_precipitate;
  RankTwoTensor _eigenstrain;

  // number of orientation variants for the precipitate in a single matrix crystal
  unsigned int _n_variants;

  //vectors to hold each individual rotated tensor variant
  //std::vector<ElasticityTensorR4> _Cijkl_precipitates_rotated;
  std::vector<RankTwoTensor> _eigenstrains_rotated;

  // turn all these tensors into materials properties
  MaterialProperty<RankTwoTensor> & _local_strain; //total strain
  MaterialProperty<RankTwoTensor> & _misfit_strain; //sum of all the active eigenstrains in space
  MaterialProperty<std::vector<RankTwoTensor> > & _eigenstrains_MP; //whether each variant is active or not in space
  MaterialProperty<ElasticityTensorR4> & _Cijkl_MP; //holds the matrix Cijkl everywhere
  MaterialProperty<ElasticityTensorR4> & _Cijkl_precipitates_MP; //holds the precipitates Cijkl everywhere

  // derivatives of the local elasticity tensor and the misfit strain with respect
  // to order parameter
  //MaterialProperty<std::vector<ElasticityTensorR4> > & _d_elasticity_tensor;
  MaterialProperty<std::vector<RankTwoTensor> > & _d_eigenstrains_MP;//derivatives of the eigenstrains operating in space

  // These have been added to allow the auxkernel to access the elasticity tensor and
  // rotated eigenstrain values for precipitates
   MaterialProperty<std::vector<RankTwoTensor> > & _precipitate_eigenstrain;

  // Vector of references to the coupled order parameters
  std::vector<VariableValue *> _coupled_variables;
};

#endif //LINEARSINGLECRYSTALPRECIPITATEMATERIAL_H

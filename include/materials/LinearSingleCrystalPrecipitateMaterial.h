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

#include "SolidMechanicsMaterial.h"
#include "SymmTensor.h"
#include "SymmAnisotropicElasticityTensor.h"

/**
 * LinearSingleCrystalPrecipitateMaterial handles anisotropic, single-crystal material elastic
 * constants.  It's designed to work with SymmSingleCrystalMaterial.  It handles a single
 * crystal of matrix with an arbitrary number of orientation variants of a coherent precipitate.
 * The rotation of orientation variants is currently set up around the C-axis (2D rotation).
 */

//Forward declaration
class LinearSingleCrystalPrecipitateMaterial;

template<>
InputParameters validParams<LinearSingleCrystalPrecipitateMaterial>();

class LinearSingleCrystalPrecipitateMaterial : public SolidMechanicsMaterial
{
public:
  LinearSingleCrystalPrecipitateMaterial(const std:: string & name, InputParameters parameters);
  virtual ~LinearSingleCrystalPrecipitateMaterial();

protected:
  virtual void computeQpProperties();
  
  virtual void computeQpElasticityTensor();
  virtual void computeQpElasticStrain();
  virtual void computeQpElasticStress();

private:
  // vectors to get the input values
  std::vector<Real> _Cijkl_matrix_vector;
  std::vector<Real> _Cijkl_precipitate_vector;
  std::vector<Real> _eigenstrain_vector;
  
  // number of orientation variants for the precipitate in a single matrix crystal
  int _n_variants;
  
  // Individual material information
  SymmAnisotropicElasticityTensor _Cijkl_matrix;
  std::vector<SymmAnisotropicElasticityTensor > _Cijkl_precipitates_rotated;
  std::vector<SymmTensor > _eigenstrains_rotated;

  MaterialProperty<SymmTensor> & _local_strain;
  MaterialProperty<SymmTensor> & _misfit_strain;
//  MaterialProperty<std::vector<SymmTensor > > & _eigenstrains_rotated_MP;
  // MaterialProperty<SymmAnisotropicElasticityTensor > & _Cijkl_matrix_MP;
  // MaterialProperty<std::vector<SymmAnisotropicElasticityTensor > > & _Cijkl_precipitates_rotated_MP;

  // Vector of references to the coupled order parameters
  std::vector<VariableValue *> _coupled_variables;
};

#endif //LINEARSINGLECRYSTALPRECIPITATEMATERIAL_H

/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  13 February 2014
*
*************************************************************************/

#ifndef PRECIPITATEMATRIXMISFITMATERIAL_H
#define PRECIPITATEMATRIXMISFITMATERIAL_H

#include "LinearSingleCrystalPrecipitateMaterial.h"

//forward declaration
class PrecipitateMatrixMisfitMaterial;

template<>
InputParameters validParams<PrecipitateMatrixMisfitMaterial>();

class PrecipitateMatrixMisfitMaterial : public LinearSingleCrystalPrecipitateMaterial
{
public:
  PrecipitateMatrixMisfitMaterial(const std::string & name, InputParameters parameters);

protected:
  virtual void computeProperties();

  virtual void computeQpElasticityTensor();
  virtual void computeQpEigenstrain();
  virtual void computeQpPrecipitateEigenstrain();
  virtual void computeQpMatrixEigenstrain();

  virtual void computeQpElasticStrain();
  virtual void computeQpMisfitStrain();
  //virtual void computeQpElasticStress();
  //virtual void computeQpStrain();
  //virtual void computeQpStress();

private:

  std::vector<Real> _eigenstrain_matrix_vector;
  RankTwoTensor _eigenstrain_matrix;

  MaterialProperty<RankTwoTensor> & _eigenstrain_matrix_MP;
  MaterialProperty<std::vector<RankTwoTensor> > & _dn_eigenstrain_matrix_MP;
  MaterialProperty<RankTwoTensor> & _dc_eigenstrain_matrix_MP;
  MaterialProperty<std::vector<ElasticityTensorR4> > & _dn_elasticity_tensor;

  MaterialProperty<std::vector<RankTwoTensor> > & _dn_misfit_strain;
  MaterialProperty<RankTwoTensor> & _dc_misfit_strain;

  VariableValue & _solute;
};

#endif //PRECIPITATEMATRIXMISFITMATERIAL_H


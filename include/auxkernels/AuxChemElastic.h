/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  6 February 2013
*
*************************************************************************/

#ifndef AUXCHEMELASTIC_H
#define AUXCHEMELASTIC_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "ElasticityTensorR4.h"


//forward declaration
class AuxChemElastic;

template<>
InputParameters validParams<AuxChemElastic>();

class AuxChemElastic : public AuxKernel
{
public:
  AuxChemElastic(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  virtual Real computeEnergy(Real & conserved, Real & nonconserved, bool matrix);
  virtual Real computeDifferential(Real & coupled_conserved,Real & nonconserved);

  virtual Real computeFchem(Real & conserved, Real & nonconserved);
  virtual Real computeSelfElasticEnergy(bool matrix);
  virtual Real computeInteractionElasticEnergy(bool matrix);

  virtual Real computeDfchemDcons(Real & coupled_conserved, Real & coupled_nonconserved);
  virtual Real computeDselfDcons();
  virtual Real computeDintDcons();

  virtual Real computeDfchemDnoncons(Real & coupled_conserved, Real & coupled_nonconserved);
  virtual Real computeDselfDnoncons();
  virtual Real computeDintDnoncons();

  VariableValue & _coupled_cons;
  VariableValue & _coupled_noncons;

  Real _precip_conserved;
  Real _precip_nonconserved;

  MaterialProperty<Real> & _a1;
  MaterialProperty<Real> & _a2;
  MaterialProperty<Real> & _a3;
  MaterialProperty<Real> & _a4;
  MaterialProperty<Real> & _c1;
  MaterialProperty<Real> & _c2;

  // unsigned int _n_variants;
  unsigned int _noncons_var_num;

  MaterialProperty<std::vector<RankTwoTensor> > & _eigenstrains_rotated_MP;
  MaterialProperty<ElasticityTensorR4> & _elasticity_tensor;

  MaterialProperty<std::vector<RankTwoTensor> > & _precipitate_eigenstrain_rotated;
  MaterialProperty<ElasticityTensorR4> & _precipitate_elasticity;

  MaterialProperty<RankTwoTensor> & _local_strain;

  MaterialProperty<std::vector<RankTwoTensor> > & _d_eigenstrains_rotated_MP;

private:
};

#endif //AUXCHEMELASTIC_H

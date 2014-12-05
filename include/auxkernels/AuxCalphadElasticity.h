/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  22 November 2014
*
*************************************************************************/

#ifndef AUXCALPHADELASTICITY_H
#define AUXCALPHADELASTICITY_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "ElasticityTensorR4.h"

#include "CalphadAB1CD1.h"
#include "CalphadAB1CD2.h"

#include "ElementIntegralVariablePostprocessor.h"

//forward declaration
class AuxCalphadElasticity;

template<>
InputParameters validParams<AuxCalphadElasticity>();

class AuxCalphadElasticity : public AuxKernel
{
public:
  AuxCalphadElasticity(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  virtual void computeHeaviside();
  virtual void computeBarrier();

  virtual void computeDHeaviside();
  virtual void computeDBarrier();

  virtual Real computeChemMatrixEnergy();
  virtual Real computeChemPrecipEnergy();
  virtual Real computeDifferential();
  virtual Real computeElasticEnergy();

  Real _precip_cons;
  Real _precip_noncons;
  Real _self_energy;

  MaterialProperty<Real> & _G_alpha;
  MaterialProperty<Real> & _G_delta;
  MaterialProperty<Real> & _G_alpha_precip;
  MaterialProperty<Real> & _G_delta_precip;
  MaterialProperty<Real> & _dG_alpha;
  MaterialProperty<Real> & _dG_delta;

  MaterialProperty<RankTwoTensor> & _stress;
  MaterialProperty<std::vector<RankTwoTensor> > & _precipitate_eigenstrain;

  MaterialProperty<Real> & _Omega;
  MaterialProperty<Real> & _W;

  unsigned int _OP_number;
  unsigned int _n_OP_vars;
  std::vector<VariableValue *> _OP;

  VariableValue & _X;
  Real _scaling_factor;
  bool _use_elastic_energy;

  Real _H;
  Real _g;
  Real _dH_dOP;
  Real _dg_dOP;

private:
};

#endif //AUXCALPHADELASTICITY_H

/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  27 December 2013
*
*************************************************************************/

#ifndef AUXCALPHADENERGY_H
#define AUXCALPHADENERGY_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "ElasticityTensorR4.h"

#include "CalphadAB1CD1.h"
#include "CalphadAB1CD2.h"

//forward declaration
class AuxCalphadEnergy;

template<>
InputParameters validParams<AuxCalphadEnergy>();

class AuxCalphadEnergy : public AuxKernel
{
public:
  AuxCalphadEnergy(const std::string & name, InputParameters parameters);

protected:
   virtual Real computeValue();

  //virtual Real computeEnergy(Real & conserved, Real & nonconserved, bool matrix);
  virtual Real computeDifferential();
  virtual Real computeMatrixEnergy();
  virtual Real computePrecipEnergy();

  virtual void computeHeaviside();
  virtual void computeBarrier();

  virtual void computeDHeaviside();
  virtual void computeDBarrier();

  Real _precip_cons;
  Real _precip_noncons;

  MaterialProperty<Real> & _G_alpha;
  MaterialProperty<Real> & _G_delta;
  MaterialProperty<Real> & _G_alpha_precip;
  MaterialProperty<Real> & _G_delta_precip;
  MaterialProperty<Real> & _dG_alpha;
  MaterialProperty<Real> & _dG_delta;

  MaterialProperty<std::vector<RankTwoTensor> > & _dn_misfit_strain;
  MaterialProperty<RankTwoTensor> & _dc_misfit_strain;
  MaterialProperty<RankTwoTensor> & _elastic_strain;
  MaterialProperty<RankTwoTensor> & _local_strain;

  MaterialProperty<ElasticityTensorR4> & _elasticity_tensor;
  MaterialProperty<ElasticityTensorR4> & _Cijkl_MP;  //matrix
  MaterialProperty<ElasticityTensorR4> & _Cijkl_precipitate_MP; //precipitate


  MaterialProperty<std::vector<RankTwoTensor> > & _precipitate_eigenstrain;
  MaterialProperty<RankTwoTensor> _matrix_eigenstrain;

  MaterialProperty<Real> & _Omega;
  MaterialProperty<Real> & _W;

  unsigned int _OP_number;
  unsigned int _n_OP_vars;
  std::vector<VariableValue *> _OP;

  VariableValue & _X;

  Real _H;
  Real _g;
  Real _dH_dOP;
  Real _dg_dOP;

private:
};

#endif //AUXCALPHADENERGY_H

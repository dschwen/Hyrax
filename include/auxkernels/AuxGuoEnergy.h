/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  27 December 2013
*
*************************************************************************/

#ifndef AUXGUOENERGY_H
#define AUXGUOENERGY_H

#include "AuxChemElastic.h"
#include "RankTwoTensor.h"
#include "ElasticityTensorR4.h"

//forward declaration
class AuxGuoEnergy;

template<>
InputParameters validParams<AuxGuoEnergy>();

class AuxGuoEnergy : public AuxChemElastic
{
public:
  AuxGuoEnergy(const std::string & name, InputParameters parameters);

protected:
  // virtual Real computeValue();

  //virtual Real computeEnergy(Real & conserved, Real & nonconserved, bool matrix);
  //virtual Real computeDifferential(Real & coupled_conserved,Real & nonconserved);

  virtual Real computeFchem(Real & conserved, Real & nonconserved);
  //virtual Real computeSelfElasticEnergy(bool matrix);
  //virtual Real computeInteractionElasticEnergy(bool matrix);

  virtual Real computeDfchemDcons(Real & coupled_conserved, Real & coupled_nonconserved);
  virtual Real computeDselfDcons();
  virtual Real computeDintDcons();

  virtual Real computeDfchemDnoncons(Real & coupled_conserved, Real & coupled_nonconserved);
  //virtual Real computeDselfDnoncons();
  //virtual Real computeDintDnoncons();

  /*
  VariableValue & _coupled_cons;
  VariableValue & _coupled_noncons;

  Real _precip_conserved;
  Real _precip_nonconserved;
  */


  // unsigned int _n_variants;
  /*
  unsigned int _noncons_var_num;

  MaterialProperty<std::vector<RankTwoTensor> > & _eigenstrains_rotated_MP;
  MaterialProperty<ElasticityTensorR4> & _elasticity_tensor;

  MaterialProperty<std::vector<RankTwoTensor> > & _precipitate_eigenstrain_rotated;
  MaterialProperty<ElasticityTensorR4> & _precipitate_elasticity;

  MaterialProperty<RankTwoTensor> & _local_strain;

  MaterialProperty<std::vector<RankTwoTensor> > & _d_eigenstrains_rotated_MP;
  */

private:
  MaterialProperty<Real> & _a1;
  MaterialProperty<Real> & _a2;
  MaterialProperty<Real> & _a3;
  MaterialProperty<Real> & _a4;
  MaterialProperty<Real> & _c1;
  MaterialProperty<Real> & _c2;
};

#endif //AUXGUOENERGY_H

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

#include "AuxChemElastic.h"
#include "RankTwoTensor.h"
#include "ElasticityTensorR4.h"

#include "CalphadAB1CD1.h"
#include "CalphadAB1CD2.h"

//forward declaration
class AuxCalphadEnergy;

template<>
InputParameters validParams<AuxCalphadEnergy>();

class AuxCalphadEnergy : public AuxChemElastic
{
public:
  AuxCalphadEnergy(const std::string & name, InputParameters parameters);

protected:
   virtual Real computeValue();

  //virtual Real computeEnergy(Real & conserved, Real & nonconserved, bool matrix);
  virtual Real computeDifferential(Real & coupled_conserved,Real & nonconserved);

  virtual Real computeFchem(Real & conserved, Real & nonconserved);
  //virtual Real computeSelfElasticEnergy(bool matrix);
  //virtual Real computeInteractionElasticEnergy(bool matrix);

  virtual Real computeDfchemDcons(Real & coupled_conserved, Real & coupled_nonconserved);
  virtual Real computeDselfDcons();
  virtual Real computeDintDcons();

  virtual Real computeDfchemDnoncons(Real & coupled_conserved, Real & coupled_nonconserved);
  //virtual Real computeDselfDnoncons();
  //virtual Real computeDintDnoncons();

  virtual Real computeHeaviside(Real & nonconserved);
  virtual Real computeBarrier(Real & nonconserved);

  virtual Real computeDHeaviside(Real & nonconserved);
  virtual Real computeDBarrier(Real & nonconserved);

  /* std::vector<Real> _alpha_low_coeffs;
  std::vector<Real> _alpha_high_coeffs;
  std::vector<Real> _alpha_mix_coeffs;

  std::vector<Real> _delta_low_coeffs;
  std::vector<Real> _delta_high_coeffs;
  std::vector<Real> _delta_mix_coeffs;

  //vectors to hold the Rudlich-Kister polynomial coefficents
  std::vector<Real> _L0_coeffs;
    std::vector<Real> _L1_coeffs;
  */


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
  /*
  MaterialProperty<Real> & _a1;
  MaterialProperty<Real> & _a2;
  MaterialProperty<Real> & _a3;
  MaterialProperty<Real> & _a4;
  MaterialProperty<Real> & _c1;
  MaterialProperty<Real> & _c2;
  */

  Real _R;                                      //Universal gas constant

  //COUPLED VARIABLES
  VariableValue & _T;                           //coupled Temperature field

  MaterialProperty<Real> & _Omega;
  MaterialProperty<Real> & _W;

  Real _scaling_factor;

//  VariableValue & _T;

  CalphadAB1CD1 _alpha;
  CalphadAB1CD2 _delta;

  std::vector<Real> _hcp_Zr_coeffs;
  std::vector<Real> _hcp_ZrH_coeffs;

  std::vector<Real> _fcc_Zr_coeffs;
  std::vector<Real> _fcc_ZrH2_coeffs;

  std::vector<Real> _H2_coeffs;

  std::vector<Real> _L0_coeffs;
  std::vector<Real> _L1_coeffs;
};

#endif //AUXCALPHADENERGY_H

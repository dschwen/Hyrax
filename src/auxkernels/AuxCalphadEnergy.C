/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  27 December 2013
*
*************************************************************************/

#include "AuxCalphadEnergy.h"
#include <iostream>

template<>
InputParameters validParams<AuxCalphadEnergy>()
{
  InputParameters params = validParams<AuxChemElastic>();

/*
  params.addRequiredCoupledVar("coupled_conserved_var", "coupled conserved field variable");
  params.addRequiredCoupledVar("coupled_nonconserved_var", "coupled non-conserved field variable");
  params.addRequiredParam<Real>("precip_conserved", "value of the equilibrium 2nd phase conserved field variable");
  params.addRequiredParam<Real>("precip_nonconserved", "value of the equilibrium 2nd phase nonconserved field variable");
  params.addRequiredParam<int>("nonconserved_var_number", "the number (starting from 1) of the nonconserved variable");
*/

  params.addParam<Real>("scaling_factor", 1, "free energy scaling factor for nondimensionalization");

  //Default in J/mol-K
  params.addParam<Real>("gas_constant", 8.3144621, "Universal gas constant");

  params.addRequiredCoupledVar("coupled_temperature", "temperature to be used to calculating Gibbs energies");

  params.addRequiredParam<std::vector<Real> >("hcp_Zr_coeffs", " A, B, C, D, E Gibbs coeffs for HCP Zr");
  params.addRequiredParam<std::vector<Real> >("H2_coeffs", "A, B, C, D, E Gibbs coeffs for H2 gas");
  params.addRequiredParam<std::vector<Real> >("hcp_ZrH_coeffs", "A, B, C, D, E Gibbs coeffs for HCP ZrH");
  params.addRequiredParam<std::vector<Real> >("fcc_Zr_coeffs", " A, B, C, D, E Gibbs coeffs for FCC Zr");
  params.addRequiredParam<std::vector<Real> >("fcc_ZrH2_coeffs", "A, B, C, D, E Gibbs coeffs for FCC ZrH2");

  params.addRequiredParam<std::vector<Real> >("L0_coeffs", "Rudlich-Kister L0 polynomial coefficients");
  params.addRequiredParam<std::vector<Real> >("L1_coeffs", "Rudlich-Kister L1 polynomial coefficients");

  return params;
}

AuxCalphadEnergy::AuxCalphadEnergy(const std::string & name, InputParameters parameters) :
    AuxChemElastic(name, parameters),
    _R(getParam<Real>("gas_constant")),
    _T(coupledValue("coupled_temperature")),
    _Omega(getMaterialProperty<Real>("molar_volume")),
    _W(getMaterialProperty<Real>("well_height")),
    _scaling_factor(getParam<Real>("scaling_factor")),
    _alpha(),
    _delta(),
    _hcp_Zr_coeffs(getParam<std::vector<Real> >("hcp_Zr_coeffs")),
    _hcp_ZrH_coeffs(getParam<std::vector<Real> >("hcp_ZrH_coeffs")),
    _fcc_Zr_coeffs(getParam<std::vector<Real> >("fcc_Zr_coeffs")),
    _fcc_ZrH2_coeffs(getParam<std::vector<Real> >("fcc_ZrH2_coeffs")),
    _H2_coeffs(getParam<std::vector<Real> >("H2_coeffs")),
    _L0_coeffs(getParam<std::vector<Real> >("L0_coeffs")),
    _L1_coeffs(getParam<std::vector<Real> >("L1_coeffs"))
    /*
    _coupled_cons(coupledValue("coupled_conserved_var")),
    _coupled_noncons(coupledValue("coupled_nonconserved_var")),
    _precip_conserved(getParam<Real>("precip_conserved")),
    _precip_nonconserved(getParam<Real>("precip_nonconserved")),
*/

/*
    _noncons_var_num(getParam<int>("nonconserved_var_number")),
    _eigenstrains_rotated_MP(getMaterialProperty<std::vector<RankTwoTensor> >("eigenstrains_MP")),
    _elasticity_tensor(getMaterialProperty<ElasticityTensorR4>("elasticity_tensor")),
    _precipitate_eigenstrain_rotated(getMaterialProperty<std::vector<RankTwoTensor> >("precipitate_eigenstrain")),
    _precipitate_elasticity(getMaterialProperty<ElasticityTensorR4>("Cijkl_precipitates_MP")),
    _local_strain(getMaterialProperty<RankTwoTensor>("local_strain")),
    _d_eigenstrains_rotated_MP(getMaterialProperty<std::vector<RankTwoTensor> >("d_eigenstrains_MP"))
*/
{
  _alpha.parameterize(_R, _hcp_Zr_coeffs, _H2_coeffs, _hcp_ZrH_coeffs);
  _delta.parameterize(_R, _fcc_Zr_coeffs, _H2_coeffs, _fcc_ZrH2_coeffs,
                      _L0_coeffs, _L1_coeffs, _hcp_Zr_coeffs);
}

Real
AuxCalphadEnergy::computeValue()
{
  Real matrix_energy(0.0);
  Real precip_energy(0.0);
  Real differential(0.0);

  matrix_energy = computeEnergy(_coupled_cons[_qp], _coupled_noncons[_qp], true);
  precip_energy = computeEnergy(_precip_conserved, _precip_nonconserved, false);
  differential = computeDifferential(_coupled_cons[_qp], _coupled_noncons[_qp]);

  std::cout<<"differential = "<<differential<<std::endl;

  return (matrix_energy - precip_energy + differential);
}


/*
Real
AuxCalphadEnergy::computeEnergy(Real & conserved, Real & nonconserved, bool matrix)
{
  Real fchem(0.0);
  Real self_elastic_energy(0.0);
  Real interaction_elastic_energy(0.0);

  fchem = computeFchem(conserved, nonconserved);
  self_elastic_energy = computeSelfElasticEnergy(matrix);
  interaction_elastic_energy = computeInteractionElasticEnergy(matrix);

  return fchem + self_elastic_energy - interaction_elastic_energy;
}
*/


Real
AuxCalphadEnergy::computeDifferential(Real & coupled_conserved, Real & coupled_nonconserved)
{
  // partial derivative of f_chem with respect to conserved variable
  Real dfchem_dcons(0.0);
  // partial derivative of f_chem with respect to nonconserved variable
  Real dfchem_dnoncons(0.0);

  // partial derivatives of self-elastic energies
  Real dself_dcons(0.0);
  Real dself_dnoncons(0.0);

  // partial derivative of interaction elastic energy
  Real dint_dcons(0.0);
  Real dint_dnoncons(0.0);

  Real first_term(0.0);
  Real second_term(0.0);

  dfchem_dcons = computeDfchemDcons(coupled_conserved, coupled_nonconserved);
  dself_dcons = computeDselfDcons();
  dint_dcons = computeDintDcons();

  dfchem_dnoncons = computeDfchemDnoncons(coupled_conserved, coupled_nonconserved);
  dself_dnoncons = computeDselfDnoncons();
  dint_dnoncons = computeDintDnoncons();

  // first_term = (dfchem_dcons + dself_dcons + dint_dcons)*(coupled_conserved - _precip_conserved);
  // second_term = (dfchem_dnoncons + dself_dnoncons + dint_dnoncons)*(coupled_nonconserved - _precip_nonconserved);

  first_term = (dfchem_dcons + dself_dcons + dint_dcons)*(_precip_conserved - coupled_conserved);
  second_term = (dfchem_dnoncons + dself_dnoncons + dint_dnoncons)*(_precip_nonconserved - coupled_nonconserved);

  std::cout<<"dfchem_dcons = "<<dfchem_dcons<<std::endl;

  return first_term + second_term;
}


Real
AuxCalphadEnergy::computeFchem(Real & conserved, Real & nonconserved)
{
  Real heaviside = computeHeaviside(nonconserved);
  Real g = computeBarrier(nonconserved);

  Real G_alpha;
  Real G_delta;

  //set the cutoff values for the energy as done in the Calphad kernels
  // if(conserved < 0.001)
  //  G_alpha = computeGalpha(0.001);
  //else if (conserved > 0.499)
  //  G_alpha = computeGalpha(0.499);
  //else
  G_alpha = _alpha.computeGMix(conserved, _T[_qp]);

    //if (conserved < 0.100)
    //G_delta = computeGdelta(0.001);
    //else if (conserved > 0.655)
    //G_delta = computeGdelta(0.655);
    //else
  G_delta = _delta.computeGMix(conserved, _T[_qp]);

  return ( (1-heaviside)*G_alpha + heaviside*G_delta + _W[_qp]*g) / _Omega[_qp];


  //this needs to be the actual phase field free energy, including the order parameter...
  //remember, this energy gets computed for each order parameter separately.
  // f = ( 1-h(n) )*gmix_alpha(x,T) + h(n)*gmix_delta(x,T) + w*g(n)
}

Real
AuxCalphadEnergy::computeHeaviside(Real & nonconserved)
{
  return 3*std::pow(nonconserved, 2) - 2*std::pow(nonconserved, 3);
}

Real
AuxCalphadEnergy::computeBarrier(Real & nonconserved)
{
  return std::pow(nonconserved, 2) - 2*std::pow(nonconserved, 3) + std::pow(nonconserved, 4);
}

Real
AuxCalphadEnergy::computeDHeaviside(Real & nonconserved)
{
  return 6*nonconserved*(1 - nonconserved);
}

Real
AuxCalphadEnergy::computeDBarrier(Real & nonconserved)
{
  return 2*nonconserved*(1 - 3*nonconserved + 2*std::pow(nonconserved, 2) );
}

/*
Real
AuxCalphadEnergy::computeSelfElasticEnergy(bool matrix)
{
  RankTwoTensor eigenstrain;
  RankTwoTensor c;
  ElasticityTensorR4 elasticity;

  if(matrix)
  {
    eigenstrain = (_eigenstrains_rotated_MP[_qp])[_noncons_var_num-1];
    elasticity = _elasticity_tensor[_qp];
  }
  else
  {
    eigenstrain = (_precipitate_eigenstrain_rotated[_qp])[_noncons_var_num-1];
    elasticity = _precipitate_elasticity[_qp];
  }

  c = elasticity*eigenstrain;

  return 0.5*c.doubleContraction(eigenstrain);
}
*/

/*
Real
AuxCalphadEnergy::computeInteractionElasticEnergy(bool matrix)
{
  RankTwoTensor eigenstrain;
  RankTwoTensor c;
  ElasticityTensorR4 elasticity;

  if(matrix)
  {
    eigenstrain = (_eigenstrains_rotated_MP[_qp])[_noncons_var_num-1];
    elasticity = _elasticity_tensor[_qp];
  }
  else
  {
    eigenstrain = (_precipitate_eigenstrain_rotated[_qp])[_noncons_var_num-1];
    elasticity = _precipitate_elasticity[_qp];
  }

  c = elasticity*eigenstrain;

  return c.doubleContraction(_local_strain[_qp]);
}
*/

Real
AuxCalphadEnergy::computeDfchemDcons(Real & coupled_conserved, Real & coupled_nonconserved)
{
  Real heaviside = computeHeaviside(coupled_nonconserved);
  Real g = computeBarrier(coupled_nonconserved);

  Real dG_alpha;
  Real dG_delta;

  dG_alpha = _alpha.computeDGMixDc(coupled_conserved, _T[_qp]);

  dG_delta = _delta.computeDGMixDc(coupled_conserved, _T[_qp]);

  std::cout<<"dG_alpha = "<<dG_alpha<<std::endl;

  return ( (1-heaviside)*dG_alpha + heaviside*dG_delta + _W[_qp]*g) / _Omega[_qp];
}

Real
AuxCalphadEnergy::computeDselfDcons()
{
  //this is legitimately 0 for the current formulation

  return 0.0;
}

Real
AuxCalphadEnergy::computeDintDcons()
{
// this is legitimately 0 for the current formulation
  return 0.0;
}

Real
AuxCalphadEnergy::computeDfchemDnoncons(Real & coupled_conserved, Real & coupled_nonconserved)
{
  Real Dheaviside = computeDHeaviside(coupled_nonconserved);
  Real Dg = computeDBarrier(coupled_nonconserved);

  Real G_alpha;
  Real G_delta;

  G_alpha = _alpha.computeGMix(coupled_conserved, _T[_qp]);

  G_delta = _delta.computeGMix(coupled_conserved, _T[_qp]);

  return ( Dheaviside*(G_delta - G_alpha) + _W[_qp]*Dg) / _Omega[_qp];
}


/*
Real
AuxCalphadEnergy::computeDselfDnoncons()
{

  RankTwoTensor eigenstrain;
  RankTwoTensor d_eigenstrain;
  RankTwoTensor c;
  ElasticityTensorR4 elasticity;

  eigenstrain = (_eigenstrains_rotated_MP[_qp])[_noncons_var_num-1];
  d_eigenstrain =( _d_eigenstrains_rotated_MP[_qp])[_noncons_var_num-1];
  elasticity = _elasticity_tensor[_qp];

  c = elasticity*eigenstrain;

  return 2.0*c.doubleContraction(d_eigenstrain);
}
*/

 /*
Real
AuxCalphadEnergy::computeDintDnoncons()
{
  RankTwoTensor d_eigenstrain;
  RankTwoTensor c;
  ElasticityTensorR4 elasticity;

  d_eigenstrain = (_precipitate_eigenstrain_rotated[_qp])[_noncons_var_num-1];
  elasticity = _precipitate_elasticity[_qp];

  c = elasticity*d_eigenstrain;

  return -2.0*c.doubleContraction(_local_strain[_qp]);
}
*/

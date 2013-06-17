/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  13 June 2013
*
*************************************************************************/

#include "ZrHCalphad.h"

template<>
InputParameters validParams<ZrHCalphad>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredParam<std::vector<Real> >("G_hcp_Zr_vector", "A, B, C, D, E Gibbs coefficients: hcp Zr");
  params.addRequiredParam<std::vector<Real> >("G_fcc_Zr_vector", "A, B, C, D, E Gibbs coefficients: fcc Zr");
  params.addRequiredParam<std::vector<Real> >("G_hcp_ZrH_vector", "A, B Gibbs coefficients: hcp ZrH");
  params.addRequiredParam<std::vector<Real> >("G_fcc_ZrH2_vector", "A, B, C Gibbs coefficients: fcc ZrH2");
  params.addRequiredParam<std::vector<Real> >("G_H2_vector", "A, B, C, D, E Gibbs coefficients: H2 gas");
  params.addRequiredParam<std::vector<Real> >("L0_vector", "A, B Gibbs coefficients: L0 term in fcc ZrH2");
  params.addRequiredParam<std::vector<Real> >("L1_vector", "A, B Gibbs coefficients: L1 term in fcc ZrH2");
  params.addRequiredParam<Real>("H_diffusion_coeff", "Diffusion coefficient for H in hcp Zr (isotropic)");
  params.addRequiredParam<Real>("mobility_AC", "mobility coefficient for Allen-Cahn equation (isotropic)");
  params.addRequiredParam<Real>("kappa_CH", "CH gradient energy coefficient (isotropic)");
  params.addRequiredParam<Real>("kappa_AC", "AC gradient energy coefficient (isotropic)");

  return params;
}

ZrHCalphad::ZrHCalphad(const std::string & name, InputParameters parameters)
    : Material(name, parameters),
      _G_hcp_Zr_vector(getParam<std::vector<Real> >("G_hcp_Zr_vector")),
      _G_fcc_Zr_vector(getParam<std::vector<Real> >("G_fcc_Zr_vector")),
      _G_hcp_ZrH_vector(getParam<std::vector<Real> >("G_hcp_ZrH_vector")),
      _G_fcc_ZrH2_vector(getParam<std::vector<Real> >("G_fcc_ZrH2_vector")),
      _G_H2_vector(getParam<std::vector<Real> >("G_H2_vector")),
      _L0_vector(getParam<std::vector<Real> >("L0_vector")),
      _L1_vector(getParam<std::vector<Real> >("L1_vector")),
      _H_coeff(getParam<Real>("H_diffusion_coeff")),
      _kappa_CH(getParam<Real>("kappa_CH")),
      _mobility_AC(getParam<Real>("mobility_AC")),
      _kappa_AC(getParam<Real>("kappa_AC")),
      _G_hcp_Zr_a(declareProperty<Real>("G_hcp_Zr_a")),
      _G_hcp_Zr_b(declareProperty<Real>("G_hcp_Zr_b")),
      _G_hcp_Zr_c(declareProperty<Real>("G_hcp_Zr_c")),
      _G_hcp_Zr_d(declareProperty<Real>("G_hcp_Zr_d")),
      _G_hcp_Zr_e(declareProperty<Real>("G_hcp_Zr_e")),
      _G_fcc_Zr_a(declareProperty<Real>("G_fcc_Zr_a")),
      _G_fcc_Zr_b(declareProperty<Real>("G_fcc_Zr_b")),
      _G_fcc_Zr_c(declareProperty<Real>("G_fcc_Zr_c")),
      _G_fcc_Zr_d(declareProperty<Real>("G_fcc_Zr_d")),
      _G_fcc_Zr_e(declareProperty<Real>("G_fcc_Zr_e")),
      _G_hcp_ZrH_a(declareProperty<Real>("G_hcp_ZrH_a")),
      _G_hcp_ZrH_b(declareProperty<Real>("G_hcp_ZrH_b")),
      _G_fcc_ZrH2_a(declareProperty<Real>("G_fcc_ZrH2_a")),
      _G_fcc_ZrH2_b(declareProperty<Real>("G_fcc_ZrH2_b")),
      _G_fcc_ZrH2_c(declareProperty<Real>("G_fcc_ZrH2_c")),
      _G_H2_a(declareProperty<Real>("G_H2_a")),
      _G_H2_b(declareProperty<Real>("G_H2_b")),
      _G_H2_c(declareProperty<Real>("G_H2_c")),
      _G_H2_d(declareProperty<Real>("G_H2_d")),
      _G_H2_e(declareProperty<Real>("G_H2_e")),
      _L0_a(declareProperty<Real>("L0_a")),
      _L0_b(declareProperty<Real>("L0_b")),
      _L1_a(declareProperty<Real>("L1_a")),
      _L1_b(declareProperty<Real>("L1_b")),
      _M(declareProperty<Real>("M")),
      _grad_M(declareProperty<RealGradient>("grad_M")),
      _kappa_c(declareProperty<Real>("kappa_c")),
      _L(declareProperty<Real>("L")),
      _kappa_n(declareProperty<Real>("kappa_n"))
{
  if((_G_hcp_Zr_vector.size() != 5) ||
     (_G_fcc_Zr_vector.size() != 5) ||
     (_G_hcp_ZrH_vector.size() != 2) ||
     (_G_fcc_ZrH2_vector.size() != 3) ||
     (_G_H2_vector.size() != 5) ||
     (_L0_vector.size() != 2) ||
     (_L1_vector.size() != 2) )
    mooseError("Please supply the correct # of values for the Gibbs coefficients (ZrHCalphad).");
}

void
ZrHCalphad::computeProperties()
{
  for(unsigned int qp = 0; qp < _qrule->n_points(); qp++)
  {
    //THIS NEEDS TO BE CHANGED
    _M[qp] = _H_coeff;
    _grad_M[qp] = 0.0;

    _L[qp] = _mobility_AC;
    _kappa_c[qp] = _kappa_CH;
    _kappa_n[qp] = _kappa_AC;

    _G_hcp_Zr_a[qp] = _G_hcp_Zr_vector[0];
    _G_hcp_Zr_b[qp] = _G_hcp_Zr_vector[1];
    _G_hcp_Zr_c[qp] = _G_hcp_Zr_vector[2];
    _G_hcp_Zr_d[qp] = _G_hcp_Zr_vector[3];
    _G_hcp_Zr_e[qp] = _G_hcp_Zr_vector[4];

    _G_fcc_Zr_a[qp] = _G_fcc_Zr_vector[0];
    _G_fcc_Zr_b[qp] = _G_fcc_Zr_vector[1];
    _G_fcc_Zr_c[qp] = _G_fcc_Zr_vector[2];
    _G_fcc_Zr_d[qp] = _G_fcc_Zr_vector[3];
    _G_fcc_Zr_e[qp] = _G_fcc_Zr_vector[4];

    _G_hcp_ZrH_a[qp] = _G_hcp_ZrH_vector[0];
    _G_hcp_ZrH_b[qp] = _G_hcp_ZrH_vector[1];

   _G_fcc_ZrH2_a[qp] = _G_fcc_ZrH2_vector[0];
   _G_fcc_ZrH2_b[qp] = _G_fcc_ZrH2_vector[1];
   _G_fcc_ZrH2_c[qp] = _G_fcc_ZrH2_vector[2];

   _G_H2_a[qp] = _G_H2_vector[0];
   _G_H2_b[qp] = _G_H2_vector[1];
   _G_H2_c[qp] = _G_H2_vector[2];
   _G_H2_d[qp] = _G_H2_vector[3];
   _G_H2_e[qp] = _G_H2_vector[4];

   _L0_a[qp] = _L0_vector[0];
   _L0_b[qp] = _L0_vector[1];

   _L1_a[qp] = _L1_vector[0];
   _L1_b[qp] = _L1_vector[1];
  }
}



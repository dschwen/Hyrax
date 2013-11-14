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

  //params.addRequiredParam<Real>("H_Zr_D0", "Diffusion prefactor for H in hcp Zr (isotropic)");
  //params.addRequiredParam<Real>("H_ZrH2_D0", "Diffusion prefactor for H in fcc ZrH2 (isotropic)");
  //params.addRequiredParam<Real>("H_Zr_Q0", "Activation energy for H in hcp Zr");
  //params.addRequiredParam<Real>("H_ZrH2_Q0", "Activaton energy for H in fcc ZrH2");

  params.addRequiredParam<Real>("mobility_CH","mobility coefficient for Cahn-Hilliard equation (isotropic)");
  params.addRequiredParam<Real>("mobility_AC", "mobility coefficient for Allen-Cahn equation (isotropic)");
  params.addRequiredParam<Real>("kappa_CH", "CH gradient energy coefficient (isotropic)");
  params.addRequiredParam<Real>("kappa_AC", "AC gradient energy coefficient (isotropic)");
  params.addRequiredParam<Real>("well_height", "well height for double well g(n)");

  params.addRequiredParam<Real>("molar_volume", "molar volume of the material");

  params.addRequiredParam<Real>("thermal_diffusivity", "Thermal diffusivity of the material");
  params.addParam<Real>("dThermaldiffusivity_dT", 0, "derivitive of thermal diffusivity with temperature");
  params.addRequiredCoupledVar("coupled_temperature", "temperature to be used to calculating Gibbs energies");

  return params;
}

ZrHCalphad::ZrHCalphad(const std::string & name, InputParameters parameters)
    : Material(name, parameters),

      //_H_Zr_D0(getParam<Real>("H_Zr_D0")),
      //_H_ZrH2_D0(getParam<Real>("H_ZrH2_D0")),
      //_H_Zr_Q0(getParam<Real>("H_Zr_Q0")),
      //_H_ZrH2_Q0(getParam<Real>("H_ZrH2_Q0")),

      _mobility_CH(getParam<Real>("mobility_CH")),
      _mobility_AC(getParam<Real>("mobility_AC")),
      _kappa_CH(getParam<Real>("kappa_CH")),
      _kappa_AC(getParam<Real>("kappa_AC")),
      _well_height(getParam<Real>("well_height")),
      _molar_volume(getParam<Real>("molar_volume")),

      _thermal_diffusivity(getParam<Real>("thermal_diffusivity")),
      _dThermal_diffusivity_dT(getParam<Real>("dThermaldiffusivity_dT")),

      _M(declareProperty<Real>("M")),
      _grad_M(declareProperty<RealGradient>("grad_M")),
      _L(declareProperty<Real>("L")),

      _kappa_c(declareProperty<Real>("kappa_c")),
      _kappa_n(declareProperty<Real>("kappa_n")),
      _W(declareProperty<Real>("well_height")),
      _molar_vol(declareProperty<Real>("molar_volume")),

      _thermal_diff(declareProperty<Real>("thermal_diffusivity")),
      _dThermDiff_dT(declareProperty<Real>("dThermal_diffusivity_dT")),

      _temperature(coupledValue("coupled_temperature"))
{
}

void
ZrHCalphad::computeQpProperties()
{

//  _D_alpha = _H_Zr_D0*std::exp(-_H_Zr_Q0/(_R*_T));
//  _D_delta = _H_ZrH2_D0*std::exp(-_H_ZrH2_Q0/(_R*_T));
  // std::cout<<"in ZrHCalphad computeQpProperties"<<std::endl;

  _M[_qp] = _mobility_CH;
  _grad_M[_qp] = 0.0;

  _L[_qp] = _mobility_AC;

  _kappa_c[_qp] = _kappa_CH;
  _kappa_n[_qp] = _kappa_AC;

  _W[_qp] = _well_height;
  _molar_vol[_qp] = _molar_volume;

  _thermal_diff[_qp] = _thermal_diffusivity;
  _dThermDiff_dT[_qp] = _dThermal_diffusivity_dT;

//  std::cout<<"end ZrHCalphad computeQpProperties"<<std::endl;
}


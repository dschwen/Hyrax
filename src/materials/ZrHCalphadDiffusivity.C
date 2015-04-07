/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  13 February 2014
*
*************************************************************************/

#include "ZrHCalphadDiffusivity.h"

template<>
InputParameters validParams<ZrHCalphadDiffusivity>()
{
  InputParameters params = validParams<ZrHCalphad>();

  params.addRequiredParam<Real>("H_Zr_D0", "Diffusion prefactor for H in hcp Zr (isotropic)");
  params.addRequiredParam<Real>("H_ZrH2_D0", "Diffusion prefactor for H in fcc ZrH2 (isotropic)");
  params.addRequiredParam<Real>("H_Zr_Q0", "Activation energy for H in hcp Zr");
  params.addRequiredParam<Real>("H_ZrH2_Q0", "Activaton energy for H in fcc ZrH2");
  params.addParam<Real>("R", 8.3144, "gas constant");
  params.addParam<Real>("k", 1.38E-23, "Boltzmann constant");

  params.addRequiredParam<int>("n_OP_variables", "# of coupled OP variables, >=1");
  params.addRequiredCoupledVar("OP_variable_names", "Array of coupled OP variable names");
  params.addRequiredCoupledVar("concentration", "coupled concentration variable");
  params.addParam<Real>("CH_mobility_scaling", 1, "scaling factor to divide by to nondimensionalize mobility");
  params.addParam<Real>("Q_transport", 0, "heat of transport of H in hcp Zr");

  return params;
}

ZrHCalphadDiffusivity::ZrHCalphadDiffusivity(const std::string & name, InputParameters parameters)
    : ZrHCalphad(name, parameters),

      _H_Zr_D0(getParam<Real>("H_Zr_D0")),
      _H_ZrH2_D0(getParam<Real>("H_ZrH2_D0")),
      _H_Zr_Q0(getParam<Real>("H_Zr_Q0")),
      _H_ZrH2_Q0(getParam<Real>("H_ZrH2_Q0")),
      _R(getParam<Real>("R")),
      _k(getParam<Real>("k")),
      _mobility_CH_scaling(getParam<Real>("CH_mobility_scaling")),
      _d2Galpha_dc2(getMaterialProperty<Real>("d2GAB1CD1_dc2")),
      _d2Gdelta_dc2(getMaterialProperty<Real>("d2GAB1CD2_dc2")),
      _d2Gdelta_dc2_precip(getMaterialProperty<Real>("d2GAB1CD2_dc2_precip")),
      _D_alpha(declareProperty<Real>("D_alpha")),
      _D_delta(declareProperty<Real>("D_delta")),
      _n_OP_variables(getParam<int>("n_OP_variables")),
      _c(coupledValue("concentration")),
      _L1Q(declareProperty<Real>("L1Q")),
      _Q_transport(getParam<Real>("Q_transport")),
      _d2Galpha_dcdT(getMaterialProperty<Real>("d2GAB1CD1_dcdT"))
{
  // Create a vector of the coupled OP variables and gradients
  if(_n_OP_variables != coupledComponents("OP_variable_names"))
    mooseError("Please match the # of orientation variants to coupled OPs (CHCoupledCalphad)");

  _OP.resize(_n_OP_variables);

  for(unsigned int i=0; i< _n_OP_variables; i++)
    _OP[i] = &coupledValue("OP_variable_names", i);
}

void
ZrHCalphadDiffusivity::computeQpProperties()
{
  Real Heaviside = computeHeaviside();

  _D_alpha[_qp] = _H_Zr_D0*std::exp(-_H_Zr_Q0/(_R*_temperature[_qp]));
  _D_delta[_qp] = _H_ZrH2_D0*std::exp(-_H_ZrH2_Q0/(_R*_temperature[_qp]));

  //nondimensionalize the mobility here
  //using mobility calculated for interstitial dilute solutions
  Real solute = _c[_qp];
  if (solute < 0)
    solute = 0;

  Real OP = (*_OP[0])[_qp];
  if (OP < 0) OP = 0;
  if (OP > 1) OP = 1;

  Real heavi = 3*OP*OP - 2*OP*OP*OP;

  //_M[_qp] = ((solute*_D_alpha[_qp])/(_R*_temperature[_qp]))/_mobility_CH_scaling;

  //_M[_qp] = ((1-Heaviside)*(_D_alpha[_qp]/_d2Galpha_dc2[_qp]) + Heaviside*(_D_delta[_qp]/_d2Gdelta_dc2_precip[_qp]))/_mobility_CH_scaling;

// _M[_qp] = ((1-Heaviside)*(_D_alpha[_qp]/_d2Galpha_dc2[_qp]) + solute*(_D_delta[_qp]/_d2Gdelta_dc2_precip[_qp]))/_mobility_CH_scaling;

//  _M[_qp] = ((1-std::sqrt(OP))*(_D_alpha[_qp]/_d2Galpha_dc2[_qp]) + _D_delta[_qp]/_d2Gdelta_dc2_precip[_qp])/_mobility_CH_scaling;
//  _M[_qp] = ((1-OP)*(_D_alpha[_qp]/_d2Galpha_dc2[_qp]) + _D_delta[_qp]/_d2Gdelta_dc2_precip[_qp])/_mobility_CH_scaling;
//  _M[_qp] = ((1-OP)*(_D_alpha[_qp]/_d2Galpha_dc2[_qp]) + OP*_D_delta[_qp]/_d2Gdelta_dc2_precip[_qp])/_mobility_CH_scaling;
   _M[_qp] = ((1-heavi)*(_D_alpha[_qp]/_d2Galpha_dc2[_qp]) + heavi*_D_delta[_qp]/_d2Gdelta_dc2_precip[_qp])/_mobility_CH_scaling;




 if (_M[_qp] < 0)
   _M[_qp] = 0;

  //_console<<"earlier M = "<< (_D_alpha[_qp]/curvature)/_mobility_CH_scaling<<std::endl;
   // _console<<"curvature = "<<curvature<<std::endl;
//  _console<<"Mobility = "<<_M[_qp]<<std::endl;
//  _console<<"D_delta = "<<_D_delta[_qp]<<std::endl;
//  _console<<"D2gdelta_precip = "<<_d2Gdelta_dc2_precip[_qp]<<std::endl;
  //_console<<"mobility scaling = "<<_mobility_CH_scaling<<std::endl;

  _grad_M[_qp] = 0.0;

  _L[_qp] = _mobility_AC;

  _kappa_c[_qp] = _kappa_CH;
  _kappa_n[_qp] = _kappa_AC;

  _W[_qp] = _well_height;
  _molar_vol[_qp] = _molar_volume;

  _thermal_diff[_qp] = _thermal_diffusivity;
  _dThermDiff_dT[_qp] = _dThermal_diffusivity_dT;
}

Real
ZrHCalphadDiffusivity::computeHeaviside()
{
  Real heaviside_first(0);
  Real heaviside_second(0);

  Real OP;
  //may need to put some checking in here so that OP fixed between 0 and 1
  for(unsigned int i=0; i<_n_OP_variables; i++)
  {
    if ((*_OP[i])[_qp] < 0 )
      OP = 0;
    if ((*_OP[i])[_qp] > 1 )
      OP = 1;
    
    //heaviside_first += std::pow((*_OP[i])[_qp], 2);
    //heaviside_second += std::pow((*_OP[i])[_qp], 3);
    heaviside_first += std::pow(OP, 2);
    heaviside_second += std::pow(OP, 3);
  }

  return 3*heaviside_first - 2*heaviside_second;
}



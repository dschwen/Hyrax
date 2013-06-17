/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  13 June 2013
*
*************************************************************************/

#include "ACCoupledCalphad.h"
#include <cmath>

template<>
InputParameters validParams<ACCoupledCalphad>()
{
  InputParameters params = validParams<ACBulk>();
  params.addRequiredCoupledVar("coupled_CH_var", "the concentration to be coupled to the AC equation");
  params.addRequiredParam<int>("n_OP_vars", "# of coupled OP variables");
  params.addRequiredCoupledVar("OP_var_names", "Array of coupled OP variable names");
  params.addRequiredParam<int>("OP_number","# of the order parameter for this kernel, starting from 1");
  params.addRequiredParam<Real>("temperature", "Simulation temperature");
  params.addParam<Real>("gas_constant", 8.3144621, "Universal gas constant");
  //universal gas constant supplied here in J/mol-K
  params.addRequiredParam<Real>("well_height", "Free energy well height");

  return params;
}

ACCoupledCalphad::ACCoupledCalphad(const std::string & name, InputParameters parameters)
    : ACBulk(name, parameters),
      _G_hcp_Zr_a(getMaterialProperty<Real>("G_hcp_Zr_a")),
      _G_hcp_Zr_b(getMaterialProperty<Real>("G_hcp_Zr_b")),
      _G_hcp_Zr_c(getMaterialProperty<Real>("G_hcp_Zr_c")),
      _G_hcp_Zr_d(getMaterialProperty<Real>("G_hcp_Zr_d")),
      _G_hcp_Zr_e(getMaterialProperty<Real>("G_hcp_Zr_e")),
      _G_hcp_ZrH_a(getMaterialProperty<Real>("G_hcp_ZrH_a")),
      _G_hcp_ZrH_b(getMaterialProperty<Real>("G_hcp_ZrH_b")),
      _G_fcc_Zr_a(getMaterialProperty<Real>("G_fcc_Zr_a")),
      _G_fcc_Zr_b(getMaterialProperty<Real>("G_fcc_Zr_b")),
      _G_fcc_Zr_c(getMaterialProperty<Real>("G_fcc_Zr_c")),
      _G_fcc_Zr_d(getMaterialProperty<Real>("G_fcc_Zr_d")),
      _G_fcc_Zr_e(getMaterialProperty<Real>("G_fcc_Zr_e")),
      _G_fcc_ZrH2_a(getMaterialProperty<Real>("G_fcc_ZrH2_a")),
      _G_fcc_ZrH2_b(getMaterialProperty<Real>("G_fcc_ZrH2_b")),
      _G_fcc_ZrH2_c(getMaterialProperty<Real>("G_fcc_ZrH2_c")),
      _G_H2_a(getMaterialProperty<Real>("G_H2_a")),
      _G_H2_b(getMaterialProperty<Real>("G_H2_b")),
      _G_H2_c(getMaterialProperty<Real>("G_H2_c")),
      _G_H2_d(getMaterialProperty<Real>("G_H2_d")),
      _G_H2_e(getMaterialProperty<Real>("G_H2_e")),
      _L0_a(getMaterialProperty<Real>("L0_a")),
      _L0_b(getMaterialProperty<Real>("L0_b")),
      _L1_a(getMaterialProperty<Real>("L1_a")),
      _L1_b(getMaterialProperty<Real>("L1_b")),
      _T(getParam<Real>("temperature")),
      _R(getParam<Real>("gas_constant")),
      _w(getParam<Real>("well_height")),
      _coupled_CH_var(coupledValue("coupled_CH_var")),
      _n_OP_vars(getParam<int>("n_OP_vars")),
      _OP_number(getParam<int>("OP_number"))
{
  // Create a vector of the coupled OP variables and set = 0 the one that the kernel
  // is operating on
  if(_n_OP_vars != coupledComponents("OP_var_names"))
    mooseError("Please match the number of orientation variants to coupled OPs (ACCoupledCalphad).");

  _coupled_OP_vars.resize(_n_OP_vars);

  for(unsigned int i=0; i< _n_OP_vars; i++)
  {
    if(i == _OP_number-1)
      _coupled_OP_vars[i] = NULL;
    else
      _coupled_OP_vars[i] = &coupledValue("OP_var_names", i);
  }
}

Real
ACCoupledCalphad::computeDFDOP(PFFunctionType type)
{
  Real square_sum, quad_sum, square_mult;
  square_sum = quad_sum = 0.0;

  if (_n_OP_vars == 1)
    square_mult = 0.0;
  else
    square_mult = 1.0;

  Real dgdn, dHeavisidedn;

  // compute the coupled OP terms
  for(unsigned int i=0; i<_n_OP_vars; i++)
  {
    if(i != _OP_number-1)
    {
      square_sum += ((*_coupled_OP_vars[i])[_qp])*((*_coupled_OP_vars[i])[_qp]);

      quad_sum += ((*_coupled_OP_vars[i])[_qp])*((*_coupled_OP_vars[i])[_qp])*((*_coupled_OP_vars[i])[_qp])*((*_coupled_OP_vars[i])[_qp]);

      square_mult *= ((*_coupled_OP_vars[i])[_qp])*((*_coupled_OP_vars[i])[_qp]);
    }
  }

  dgdn = 2*_u[_qp] - 6*_u[_qp]*_u[_qp] + 4*_u[_qp]*_u[_qp]*_u[_qp]
    + 2*_u[_qp]*(square_sum + quad_sum) + 4*_u[_qp]*_u[_qp]*_u[_qp]*square_sum + 2*square_mult;

  dHeavisidedn = 6*_u[_qp]*(1 - _u[_qp]);

  switch (type)
  {
  case Residual:
    return (-1*computeGalphamix() + computeGdeltamix())*dHeavisidedn + _w*dgdn;

  case Jacobian:
    return 0;
    // return _phi[_j][_qp]*(6-6*_u[_qp])*(-1*computeGalphamix() + computeGdeltamix()) + _w*(2 - 6*_u[_qp] + 12*_u[_qp]*_u[_qp] + 2*(square_sum) + 2*(quad_sum) + 12*(square_sum) + 2*square_mult)*_phi[_j][_qp];
  }
  mooseError("Invalid type passed in");
}

Real
ACCoupledCalphad::computeGalphamix()
{
  Real reference, ideal;

  //do some checking so that equations are only calculated over the valid region
  //here, 0<atomic fraction H <= 0.5

  if(_coupled_CH_var[_qp] > 0.5 || _coupled_CH_var[_qp] < 0)
    return 0;

  reference = (1 - 2*_coupled_CH_var[_qp])*computeGhcpZr() + _coupled_CH_var[_qp]*computeGhcpZrH();

  ideal = _R*_T*( (1 - 2*_coupled_CH_var[_qp])*std::log( (2*_coupled_CH_var[_qp]-1)/(_coupled_CH_var[_qp] - 1))
                  +_coupled_CH_var[_qp]*std::log( _coupled_CH_var[_qp]/(1 - _coupled_CH_var[_qp])) );

  return  reference + ideal;
}


Real
ACCoupledCalphad::computeGdeltamix()
{
  Real reference, ideal, excess;
  Real L0, L1;

  //do some checking so that equations are only calculated over the valid region
  //here, 0<atomic fraction H <= 2/3

  if(_coupled_CH_var[_qp] > 2/3 || _coupled_CH_var[_qp] < 0)
    return 0;

  L0 = _L0_a[_qp] - _L0_b[_qp]*_T;
  L1 = _L1_a[_qp] + _L1_b[_qp]*_T;


  reference = (1/2)*( (2 - 3*_coupled_CH_var[_qp])*computeGfccZr() + _coupled_CH_var[_qp]*computeGfccZrH2() );

  ideal = _R*_T*( (1 - 2*_coupled_CH_var[_qp])*std::log( (2 - 3*_coupled_CH_var[_qp])/(2 - 2*_coupled_CH_var[_qp]))
                  +_coupled_CH_var[_qp]*std::log( _coupled_CH_var[_qp]/(2 - 2*_coupled_CH_var[_qp])) );

  excess = ((3*_coupled_CH_var[_qp] - 2)*_coupled_CH_var[_qp])/(4*std::pow(_coupled_CH_var[_qp] - 1, 2));
  excess *= ( (_coupled_CH_var[_qp] - 1)*L0 + (1 - 2*_coupled_CH_var[_qp])*L1 );

  return reference + ideal + excess;
}

Real
ACCoupledCalphad::computeGhcpZr()
{
  //Calculates molar Gibbs free energy for Zr in HCP form - molar standard enthalpy for Zr

  return _G_hcp_Zr_a[_qp] + _G_hcp_Zr_b[_qp]*_T + _G_hcp_Zr_c[_qp]*_T*std::log(_T) + _G_hcp_Zr_d[_qp]*_T*_T
    + _G_hcp_Zr_e[_qp]/_T;
}

Real
ACCoupledCalphad::computeGhcpZrH()
{
  //Calculates molar Gibbs free energy for ZrH compound - molar standard enthalpies

  return _G_hcp_ZrH_a[_qp] + _G_hcp_ZrH_b[_qp]*_T + computeGhcpZr() + 0.5*computeGH2();
}

Real
ACCoupledCalphad::computeGfccZr()
{
  //Calculates molar Gibbs free energy for fcc Zr - molar standard enthalpy for Zr

  return  _G_fcc_Zr_a[_qp] + _G_fcc_Zr_b[_qp]*_T + _G_fcc_Zr_c[_qp]*_T*std::log(_T) + _G_fcc_Zr_d[_qp]*_T*_T
    + _G_fcc_Zr_e[_qp]/_T;
}

Real
ACCoupledCalphad::computeGfccZrH2()
{
  //Calculates molar Gibbs free energy for fcc ZrH2 - molar standard enthalpies

  return _G_fcc_ZrH2_a[_qp] + _G_fcc_ZrH2_b[_qp]*_T + _G_fcc_ZrH2_c[_qp]*_T*std::log(_T) + computeGhcpZr()
    + computeGH2();
}

Real
ACCoupledCalphad::computeGH2()
{
  //Calculates molar Gibbs free energy for H2 gas - molar standard enthalpy for H2

  return _G_H2_a[_qp] + _G_H2_b[_qp]*_T + _G_H2_c[_qp]*_T*std::log(_T) + _G_H2_d[_qp]*_T*_T + _G_H2_e[_qp]/_T;
}

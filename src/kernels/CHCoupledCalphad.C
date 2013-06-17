/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*
*  11 June 2013
*
*************************************************************************/

#include "CHCoupledCalphad.h"
#include <cmath>

/**
 * CHCoupledCalphad handles the conserved order parameter (atomic fraction)
 * evolved using the Cahn-Hilliard equation coupled to an arbitrary number of
 * non-conserved structural order parameters for the Zr-delta hydride system.
 * It uses free energy equations taken from A.T. Dinsdale 15 (1991) CALPHAD 317,
 * N. Dupin et al, J. Nuc. Mater. 275 (1999) 287 and J.D. Cox et al., CODATA Key
 * Values for Thermodynamics, Hemisphere Publishing Corporation, 1989.
 */

template<>
InputParameters validParams<CHCoupledCalphad>()
{
  InputParameters params = validParams<CHBulk>();
  params.addRequiredParam<int>("n_OP_variables", "# of coupled OP variables, >=1");
  params.addRequiredCoupledVar("OP_variable_names", "Array of coupled OP variable names");
  params.addRequiredParam<Real>("temperature", "Simulation temperature");
  params.addParam<Real>("gas_constant", 8.3144621, "Universal gas constant");
  //universal gas constant supplied here in J/mol-K

  return params;
}

CHCoupledCalphad::CHCoupledCalphad(const std::string & name, InputParameters parameters)
    : CHBulk(name, parameters),
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
      _n_OP_variables(getParam<int>("n_OP_variables"))
{
  // Create a vector of the coupled OP variables and gradients
  if(_n_OP_variables != coupledComponents("OP_variable_names"))
    mooseError("Please match the # of orientation variants to coupled OPs (CHCoupledCalphad)");

  _coupled_OP_vars.resize(_n_OP_variables);
  _coupled_OP_grads.resize(_n_OP_variables);

  for(unsigned int i=0; i< _n_OP_variables; i++)
  {
    _coupled_OP_vars[i] = &coupledValue("OP_variable_names", i);
    _coupled_OP_grads[i] = &coupledGradient("OP_variable_names", i);
  }
}

RealGradient
CHCoupledCalphad::computeGradDFDCons(PFFunctionType type, Real /*c*/, RealGradient /*grad_c*/)
{
  RealGradient grad_conserved_term, grad_nonconserved_term;
  Real heaviside;

  switch (type)
  {
  case Residual:
    //calculate h(n1,n2,..np) (numerical heaviside function)
    heaviside = computeHeaviside();

    //calculate the d/dx (dfchem/dx) grad x term
    grad_conserved_term = computeGradConservedTerm(heaviside);

    //loop to calculate the d/dn (dfchem/dx) grad n terms
    grad_nonconserved_term = computeGradNonconservedTerm();

    return grad_conserved_term + grad_nonconserved_term;

  case Jacobian:
    RealGradient op_sum;
    op_sum.zero();

    Real ideal_alpha, ideal_delta, excess, L0, L1;

    for(unsigned int i=0; i<_n_OP_variables; i++)
      op_sum += 6*(*_coupled_OP_vars[i])[_qp]*(1 - (*_coupled_OP_vars[i])[_qp])*(*_coupled_OP_grads[i])[_qp];

    if(_u[_qp] > 0.5 || _u[_qp] < 0)
      ideal_alpha = 0;
    else
      ideal_alpha = _R*_T/(_u[_qp]*(1 - 3*_u[_qp] + 2*_u[_qp]*_u[_qp]));

    L0 = _L0_a[_qp] - _L0_b[_qp]*_T;
    L1 = _L1_a[_qp] + _L1_b[_qp]*_T;

    //do some checking so that equations are only calculated over the valid region
    //here, 0<atomic fraction H <= 2/3
    if(_u[_qp] > 2/3 || _u[_qp] < 0)
      ideal_delta = 0;
    else
      ideal_delta = 2*_R*_T/(_u[_qp]*(2 - 5*_u[_qp] + 3*_u[_qp]*_u[_qp]));

    excess = (L0*(_u[_qp] - 1) + L1*(3 - 6*_u[_qp])) / (2*pow(_u[_qp] - 1, 4));

    Real temp1;

    temp1 = (1-computeHeaviside())*_R*_T*(-1-6*_u[_qp]*-6*_u[_qp]*_u[_qp])/(std::pow(_u[_qp]-3*_u[_qp]*_u[_qp]+2*_u[_qp]*_u[_qp]*_u[_qp],2));

    Real temp2;

    temp2 = computeHeaviside()*(_R*_T*(-4-20*_u[_qp]-6*_u[_qp]*_u[_qp])/(std::pow(2*_u[_qp]-5*_u[_qp]*_u[_qp] +3*_u[_qp]*_u[_qp]*_u[_qp],2)) +( (L0-6*L1)*(_u[_qp]-1) - 4*(L0*(_u[_qp]-1)+3*L1*(3-6*_u[_qp])))/(2*std::pow(_u[_qp]-1,5)));

    return 0;
    // return op_sum*(-1*ideal_alpha + (ideal_delta + excess))*_phi[_j][_qp] + (temp1 + temp2)*_grad_phi[_j][_qp];
  }

  mooseError("Invalid type passed in");
}

Real
CHCoupledCalphad::computeHeaviside()
{
  Real heaviside_first(0.0);
  Real heaviside_second(0.0);

  for(unsigned int i=0; i<_n_OP_variables; i++)
  {
    heaviside_first += std::pow((*_coupled_OP_vars[i])[_qp], 2);
    heaviside_second += std::pow((*_coupled_OP_vars[i])[_qp], 3);
  }

  return 3*heaviside_first - 2*heaviside_second;
}

RealGradient
CHCoupledCalphad::computeGradConservedTerm(Real & h)
{
  Real ideal, excess;
  Real L0, L1;
  Real d2Galphadx2, d2Gdeltadx2;

  //do some checking so that equations are only calculated over the valid region
  //here, 0<atomic fraction H <= 0.5
  if(_u[_qp] > 0.5 || _u[_qp] < 0)
    ideal = 0;
  else
    ideal = _R*_T/(_u[_qp]*(1 - 3*_u[_qp] + 2*_u[_qp]*_u[_qp]));

  d2Galphadx2 = (1-h)*ideal;


  L0 = _L0_a[_qp] - _L0_b[_qp]*_T;
  L1 = _L1_a[_qp] + _L1_b[_qp]*_T;

  //do some checking so that equations are only calculated over the valid region
  //here, 0<atomic fraction H <= 2/3
  if(_u[_qp] > 2/3 || _u[_qp] < 0)
    ideal = 0;
  else
    ideal = 2*_R*_T/(_u[_qp]*(2 - 5*_u[_qp] + 3*_u[_qp]*_u[_qp]));

  excess = (L0*(_u[_qp] - 1) + L1*(3 - 6*_u[_qp])) / (2*pow(_u[_qp] - 1, 4));

  d2Gdeltadx2 = h*(ideal + excess);

  return (d2Galphadx2 + d2Gdeltadx2)*_grad_u[_qp];
}

RealGradient
CHCoupledCalphad::computeGradNonconservedTerm()
{
  RealGradient op_sum;
  op_sum.zero();

  for(unsigned int i=0; i<_n_OP_variables; i++)
  {
    op_sum += 6*(*_coupled_OP_vars[i])[_qp]*(1 - (*_coupled_OP_vars[i])[_qp])*(*_coupled_OP_grads[i])[_qp];
  }

  return (-1*computeDGalphaDx() + computeDGdeltaDx())*op_sum;
}

Real
CHCoupledCalphad::computeDGalphaDx()
{
  Real reference, ideal;

  //do some checking so that equations are only calculated over the valid region
  //here, 0<atomic fraction H <= 0.5

  if(_u[_qp] > 0.5 || _u[_qp] < 0)
    return 0;

  reference = -2*computeGhcpZr() + computeGhcpZrH();

  //std::cout<<"_u[_qp] = "<<_u[_qp]<<std::endl;
  //std::cout<<"first = "<< ( _u[_qp]/(1-_u[_qp])) <<std::endl;
  //std::cout<<" second "<<((2*_u[_qp]-1)/(_u[_qp]-1)) << std::endl;

  ideal = + _R*_T*( std::log( _u[_qp]/(1-_u[_qp])) - 2*std::log((1-2*_u[_qp])/(1-_u[_qp])) );

  return reference + ideal;
}

Real
CHCoupledCalphad::computeDGdeltaDx()
{
  Real reference, ideal, excess;
  Real L0, L1;

  //do some checking so that equations are only calculated over the valid region
  //here, 0<atomic fraction H <= 0.5

  if(_u[_qp] > 2/3 || _u[_qp] <  0)
    return 0;

  L0 = _L0_a[_qp] - _L0_b[_qp]*_T;
  L1 = _L1_a[_qp] + _L1_b[_qp]*_T;

  reference = -(3/2)*computeGfccZr() + (1/2)*computeGfccZrH2();

  //std::cout<<"_u[_qp] = "<<_u[_qp]<<std::endl;

  ideal = + _R*_T*( std::log(-4*_u[_qp]/(_u[_qp]-1))- 3*std::log((3*_u[_qp]-2)/(_u[_qp]-1)) );

  excess = ( 1/(4*std::pow(_u[_qp]-1, 3)) )*( L0*(-2 + 8*_u[_qp] - 9*_u[_qp]*_u[_qp] + 3*std::pow(_u[_qp], 3))
                                              -2*L1*(-1 + 6*_u[_qp] - 9*_u[_qp]*_u[_qp]
                                                     + 3*std::pow(_u[_qp], 3)));

  return reference + ideal + excess;;
}

Real
CHCoupledCalphad::computeGhcpZr()
{
  //Calculates molar Gibbs free energy for Zr in HCP form - molar standard enthalpy for Zr

   return _G_hcp_Zr_a[_qp] + _G_hcp_Zr_b[_qp]*_T + _G_hcp_Zr_c[_qp]*_T*std::log(_T) + _G_hcp_Zr_d[_qp]*_T*_T
    + _G_hcp_Zr_e[_qp]/_T;
}

Real
CHCoupledCalphad::computeGhcpZrH()
{
  //Calculates molar Gibbs free energy for ZrH compound - molar standard enthalpies

  return _G_hcp_ZrH_a[_qp] + _G_hcp_ZrH_b[_qp]*_T + computeGhcpZr() + 0.5*computeGH2();
}

Real
CHCoupledCalphad::computeGfccZr()
{
  //Calculates molar Gibbs free energy for fcc Zr - molar standard enthalpy for Zr

 return  _G_fcc_Zr_a[_qp] + _G_fcc_Zr_b[_qp]*_T + _G_fcc_Zr_c[_qp]*_T*std::log(_T) + _G_fcc_Zr_d[_qp]*_T*_T
    + _G_fcc_Zr_e[_qp]/_T;
}

Real
CHCoupledCalphad::computeGfccZrH2()
{
  //Calculates molar Gibbs free energy for fcc ZrH2 - molar standard enthalpies

  return _G_fcc_ZrH2_a[_qp] + _G_fcc_ZrH2_b[_qp]*_T + _G_fcc_ZrH2_c[_qp]*_T*std::log(_T) + computeGhcpZr()
    + computeGH2();
}

Real
CHCoupledCalphad::computeGH2()
{
  //Calculates molar Gibbs free energy for H2 gas - molar standard enthalpy for H2

 return _G_H2_a[_qp] + _G_H2_b[_qp]*_T + _G_H2_c[_qp]*_T*std::log(_T) + _G_H2_d[_qp]*_T*_T + _G_H2_e[_qp]/_T;
}

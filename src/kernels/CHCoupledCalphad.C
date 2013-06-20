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
  Real L0, L1;

  //calculate h(n1,n2,..np) (numerical heaviside function)
  heaviside = computeHeaviside();

  L0 = _L0_a[_qp] + _L0_b[_qp]*_T;
  L1 = _L1_a[_qp] + _L1_b[_qp]*_T;

  switch (type)
  {
  case Residual:
    //calculate the d/dx (dfchem/dx) grad x term
    //std::cout<<"in residual"<<std::endl;
    grad_conserved_term = computeGradConservedTerm(heaviside, L0, L1);

    // std::cout<<"after computeGradConservedTerm"<<std::endl;
    //loop to calculate the d/dn (dfchem/dx) grad n terms
    grad_nonconserved_term = computeGradNonconservedTerm(L0, L1);

    // std::cout<<"After computeGradConservedTerm"<<std::endl;
    return grad_conserved_term + grad_nonconserved_term;

  case Jacobian:
    RealGradient nonconserved, conserved_first, conserved_second;
    nonconserved = computeGradOPDHeavisideDOP()*_phi[_j][_qp]*(-1*computeD2GalphaDx2() + computeD2GdeltaDx2(L0, L1));

    conserved_first = (1-heaviside)*(computeD3GalphaDx3()*_phi[_j][_qp]*_grad_u[_qp] + computeD2GalphaDx2()*_grad_phi[_j][_qp]);

    conserved_second = heaviside*(computeD3GdeltaDx3(L0, L1)*_phi[_j][_qp]*_grad_u[_qp] + computeD2GdeltaDx2(L0, L1)*_grad_phi[_j][_qp]);

    return nonconserved + conserved_first + conserved_second;
  }

  mooseError("Invalid type passed in");
}

Real
CHCoupledCalphad::computeHeaviside()
{
  Real heaviside_first(0.0);
  Real heaviside_second(0.0);

  //may need to put some checking in here so that OP fixed between 0 and 1
  for(unsigned int i=0; i<_n_OP_variables; i++)
  {
    heaviside_first += std::pow((*_coupled_OP_vars[i])[_qp], 2);
    heaviside_second += std::pow((*_coupled_OP_vars[i])[_qp], 3);
  }

  return 3*heaviside_first - 2*heaviside_second;
}

RealGradient
CHCoupledCalphad::computeGradOPDHeavisideDOP()
{
  RealGradient op_sum;
  op_sum.zero();

  //may need to put some checking in here so that OP fixed between 0 and 1
  for(unsigned int i=0; i<_n_OP_variables; i++)
    op_sum += 6*(*_coupled_OP_vars[i])[_qp]*(1 - (*_coupled_OP_vars[i])[_qp])*(*_coupled_OP_grads[i])[_qp];

  return op_sum;
}

RealGradient
CHCoupledCalphad::computeGradConservedTerm(Real & h, Real & L0, Real & L1)
{
  Real d2Galphadx2, d2Gdeltadx2;

  d2Galphadx2 = (1-h)*computeD2GalphaDx2();

  d2Gdeltadx2 = h*computeD2GdeltaDx2(L0, L1);
  // std::cout<<"in computeGradConservedTerm"<<std::endl;

  return (d2Galphadx2 + d2Gdeltadx2)*_grad_u[_qp];
}

RealGradient
CHCoupledCalphad::computeGradNonconservedTerm(Real & L0, Real & L1)
{
  //  std::cout<<"in computeGradNonconservedTerm"<<std::endl;
  return (-1*computeDGalphaDx() + computeDGdeltaDx(L0, L1))*computeGradOPDHeavisideDOP();
}

Real
CHCoupledCalphad::computeDGalphaDx()
{
  Real reference, ideal;

  //do some checking so that equations are only calculated over the valid region
  //THIS MAY NEED TO BE CHANGED
  if(_u[_qp] > 0.5 || _u[_qp] < 0)
    return 0;
  else
  {
    reference = -2*computeGhcpZr() + computeGhcpZrH();

    // std::cout<<"x/(1-x) = "<< _u[_qp]/(1-_u[_qp])<<std::endl;
    // std::cout<<"(2x-1)/(x-1)"<<(2*_u[_qp] - 1)/(_u[_qp] -1)<<std::endl;

    ideal = _R*_T*( std::log( _u[_qp]/(1-_u[_qp])) - 2*std::log( (2*_u[_qp] - 1)/(_u[_qp] -1) ) );

    return reference + ideal;
  }
}

Real
CHCoupledCalphad::computeDGdeltaDx(Real & L0, Real & L1)
{
  Real reference, ideal, excess;

  //do some checking so that equations are only calculated over the valid region
  //THIS MAY NEED TO BE CHANGED
  if(_u[_qp] > 2/3 || _u[_qp] <  0)
    return 0;
  else
  {
    reference = -(3/2)*computeGfccZr() + (1/2)*computeGfccZrH2();

    // std::cout<<"-4x/(x-1) = "<< -4*_u[_qp]/(_u[_qp]-1)<<std::endl;
    // std::cout<<"(3x-2)/(x-1)"<<(3*_u[_qp] - 2)/(_u[_qp] -1)<<std::endl;

    ideal = _R*_T*( std::log(-4*_u[_qp]/(_u[_qp]-1))- 3*std::log((3*_u[_qp]-2)/(_u[_qp]-1)) );

    excess = ( 1/(4*std::pow(_u[_qp]-1, 3)) )*( L0*(-2 + 8*_u[_qp] - 9*_u[_qp]*_u[_qp]
                                                    + 3*std::pow(_u[_qp], 3))
                                                - 2*L1*(-1 + 6*_u[_qp] - 9*_u[_qp]*_u[_qp]
                                                        + 3*std::pow(_u[_qp], 3)));

    return reference + ideal + excess;
  }
}

Real
CHCoupledCalphad::computeD2GalphaDx2()
{
  //do some checking so that equations are only calculated over the valid region
  //THIS MAY CHANGE
  if(_u[_qp] > 0.5 || _u[_qp] < 0)
    return 0;
  else
    return _R*_T/(_u[_qp] - 3*_u[_qp]*_u[_qp] + 2*_u[_qp]*_u[_qp]*_u[_qp]);
}

Real
CHCoupledCalphad::computeD2GdeltaDx2(Real & L0, Real & L1)
{
  Real ideal, excess;

  //do some checking so that equations are only calculated over the valid region
  //THIS MAY CHANGE
  if(_u[_qp] > 2/3 || _u[_qp] < 0)
    return 0;
  else
  {
    ideal = 2*_R*_T/(2*_u[_qp] - 5*_u[_qp]*_u[_qp] + 3*_u[_qp]*_u[_qp]*_u[_qp]);
    excess = (L0*(_u[_qp] - 1) + L1*(3 - 6*_u[_qp])) / (2*pow(_u[_qp] - 1, 4));

    return ideal + excess;
  }
}

Real
CHCoupledCalphad::computeD3GalphaDx3()
{
 return -1*_R*_T*(1 - 6*_u[_qp] + 6*_u[_qp]*_u[_qp])/std::pow(_u[_qp] - 3*_u[_qp]*_u[_qp]
                                                      + 2*_u[_qp]*_u[_qp]*_u[_qp], 2);
}

Real
CHCoupledCalphad::computeD3GdeltaDx3(Real & L0, Real & L1)
{
  Real ideal, excess;

  ideal = -2*_R*_T*(2 - 10*_u[_qp] + 9*_u[_qp]*_u[_qp])/std::pow(2*_u[_qp] - 5*_u[_qp]*_u[_qp]
                                                      + 3*_u[_qp]*_u[_qp]*_u[_qp], 2);

  excess =  3*(L0*(1-_u[_qp]) + L1*(6*_u[_qp]-2))/(2*std::pow((_u[_qp] - 1), 5));

  return ideal + excess;
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

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
  params.addParam<Real>("temperature", 600, "Simulation temperature in K");
  params.addParam<Real>("gas_constant", 8.3144621, "Universal gas constant");
  //universal gas constant supplied here in J/mol-K

  return params;
}

CHCoupledCalphad::CHCoupledCalphad(const std::string & name, InputParameters parameters)
    : CHBulk(name, parameters),
      _G_hcp_Zr(getMaterialProperty<Real>("G_hcp_Zr")),
      _G_hcp_ZrH(getMaterialProperty<Real>("G_hcp_ZrH")),
      _G_fcc_Zr(getMaterialProperty<Real>("G_fcc_Zr")),
      _G_fcc_ZrH2(getMaterialProperty<Real>("G_fcc_ZrH2")),
      _G_H2(getMaterialProperty<Real>("G_H2")),
      _L0(getMaterialProperty<Real>("L0")),
      _L1(getMaterialProperty<Real>("L1")),
      _molarVol_alpha_Zr(getMaterialProperty<Real>("molar_volume_alpha_Zr")),
      _molarVol_delta_ZrH2(getMaterialProperty<Real>("molar_volume_delta_ZrH2")),
      _T(getParam<Real>("temperature")),
      _R(getParam<Real>("gas_constant")),
      _n_OP_variables(getParam<int>("n_OP_variables"))
{
  // Create a vector of the coupled OP variables and gradients
  if(_n_OP_variables != coupledComponents("OP_variable_names"))
    mooseError("Please match the # of orientation variants to coupled OPs (CHCoupledCalphad)");

  _OP.resize(_n_OP_variables);
  _grad_OP.resize(_n_OP_variables);

  for(unsigned int i=0; i< _n_OP_variables; i++)
  {
    _OP[i] = &coupledValue("OP_variable_names", i);
    _grad_OP[i] = &coupledGradient("OP_variable_names", i);
  }
}

RealGradient
CHCoupledCalphad::computeGradDFDCons(PFFunctionType type, Real /*c*/, RealGradient /*grad_c*/)
{
  RealGradient grad_conserved_term, grad_nonconserved_term;

  _Heaviside = computeHeaviside();
  computeDHeaviside();

 switch (type)
  {
  case Residual:
    //calculate the d/dx (dfchem/dx) grad x term
    grad_conserved_term = computeGradConservedTerm();

    //loop to calculate the d/dn (dfchem/dx) grad n terms
    grad_nonconserved_term = computeGradNonconservedTerm();

    return grad_conserved_term + grad_nonconserved_term;

  case Jacobian:
    RealGradient nonconserved, conserved_first, conserved_second;
    nonconserved.zero();
    // nonconserved = computeGradOPDHeavisideDOP()*_phi[_j][_qp]*(-1*computeD2GalphaDx2() + computeD2GdeltaDx2(L0, L1));

    conserved_first = (1-_Heaviside)*(computeD3GalphaDx3()*_phi[_j][_qp]*_grad_u[_qp] + computeD2GalphaDx2()*_grad_phi[_j][_qp]);

    conserved_second = _Heaviside*(computeD3GdeltaDx3()*_phi[_j][_qp]*_grad_u[_qp] + computeD2GdeltaDx2()*_grad_phi[_j][_qp]);

    for(unsigned int i(0); i < _n_OP_variables; ++i)
    {
      nonconserved += _dHeaviside[i]*(*_grad_OP[i])[_qp]*_phi[_j][_qp];
    }
    nonconserved *= computeD2GdeltaDx2() - computeD2GalphaDx2();

    /* nonconserved = _GradOPDHeavisideDOP[_qp]*_phi[_j][_qp]*(-1*_D2GalphaDx2[_qp] + _D2GdeltaDx2[_qp]);

    conserved_first = (1.0 - _heaviside[_qp])*(_D3GalphaDx3[_qp]*_phi[_j][_qp]*_grad_u[_qp] + _D2GalphaDx2[_qp]*_grad_phi[_j][_qp]);

    conserved_second = _heaviside[_qp]*(_D3GdeltaDx3[_qp]*_phi[_j][_qp]*_grad_u[_qp] + _D2GdeltaDx2[_qp]*_grad_phi[_j][_qp]); */

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
    if( (*_OP[i])[_qp] > 1.0)
    {
      heaviside_first += 1.;
      heaviside_second += 1.;
    }
    else if( (*_OP[i])[_qp] < 0.0)
    {
      heaviside_first += 0.;
      heaviside_second += 0.;
    }
    else
    {
      heaviside_first += std::pow((*_OP[i])[_qp], 2.);
      heaviside_second += std::pow((*_OP[i])[_qp], 3.);
    }

  }

  return 3.*heaviside_first - 2.*heaviside_second;
}

void
CHCoupledCalphad::computeDHeaviside()
{
  _dHeaviside.resize(_n_OP_variables);

  for(unsigned int i(0); i<_n_OP_variables; ++i)
  {
    //piecewise
    if( (*_OP[i])[_qp] < 0.0 || (*_OP[i])[_qp] > 1.0)
      _dHeaviside[i] = 0.0;
    else
      _dHeaviside[i] = 6.0*(*_OP[i])[_qp]*(1.0 - (*_OP[i])[_qp] );
  }
}

RealGradient
CHCoupledCalphad::computeGradConservedTerm()
{
  Real d2Galphadx2, d2Gdeltadx2;

  d2Galphadx2 = (1 - _Heaviside)*computeD2GalphaDx2();
  d2Gdeltadx2 = _Heaviside*computeD2GdeltaDx2();
  // std::cout<<"in computeGradConservedTerm"<<std::endl;

  return (d2Galphadx2 + d2Gdeltadx2)*_grad_u[_qp];
}

RealGradient
CHCoupledCalphad::computeGradNonconservedTerm()
{
  Real dGalphadx, dGdeltadx;

  dGalphadx = computeDGalphaDx();
  dGdeltadx = computeDGdeltaDx();

  RealGradient sum;
  sum.zero();

  for (unsigned int i(0); i < _n_OP_variables; ++i)
  {
    sum += ((1. - _dHeaviside[i])*dGalphadx  + _dHeaviside[i]*dGdeltadx)*(*_grad_OP[i])[_qp];
  }

  return sum;
}

Real
CHCoupledCalphad::computeDGalphaDx()
{
  Real reference, ideal;

  //do some checking so that equations are only calculated over the valid region
  //THIS MAY NEED TO BE CHANGED
  if(_u[_qp] > 0.5)
    return 0.0;
  else if( _u[_qp] < 0.0)
    return -1e20;
  else
  {
    reference = -2.*_G_hcp_Zr[_qp] + _G_hcp_ZrH[_qp];

    // std::cout<<"x/(1-x) = "<< _u[_qp]/(1-_u[_qp])<<std::endl;
    // std::cout<<"(2x-1)/(x-1)"<<(2*_u[_qp] - 1)/(_u[_qp] -1)<<std::endl;

    ideal = _R*_T*( std::log( _u[_qp]/(1-_u[_qp]))
                    - 2.*std::log( (2.*_u[_qp] - 1.)/(_u[_qp] - 1.) ) );

    return (reference + ideal)/_molarVol_alpha_Zr[_qp];
  }
}

Real
CHCoupledCalphad::computeDGdeltaDx()
{
  Real reference, ideal, excess;

  //do some checking so that equations are only calculated over the valid region
  //THIS MAY NEED TO BE CHANGED
  if(_u[_qp] > 2.0/3.0)
    return 1e20;
  else if ( _u[_qp] <  0.0)
    return 0.0;
  else
  {
    reference = -1.5*_G_fcc_Zr[_qp] + 0.5*_G_fcc_ZrH2[_qp];

    // std::cout<<"-4x/(x-1) = "<< -4*_u[_qp]/(_u[_qp]-1)<<std::endl;
    // std::cout<<"(3x-2)/(x-1)"<<(3*_u[_qp] - 2)/(_u[_qp] -1)<<std::endl;

    ideal = _R*_T*( std::log(-4.*_u[_qp]/(_u[_qp] - 1.))
                    - 3.*std::log((3.*_u[_qp] - 2.)/(_u[_qp] - 1)) );

    excess = ( 1./(4.*std::pow(_u[_qp] - 1., 3.)) )*( _L0[_qp]*(-2. + 8.*_u[_qp] - 9.*_u[_qp]*_u[_qp] + 3.*std::pow(_u[_qp], 3.))
                                                      - 2.*_L1[_qp]*(-1. + 6.*_u[_qp] - 9.*_u[_qp]*_u[_qp] + 3.*std::pow(_u[_qp], 3.)));

    return (reference + ideal + excess)/_molarVol_delta_ZrH2[_qp];
  }
}

Real
CHCoupledCalphad::computeD2GalphaDx2()
{
  //do some checking so that equations are only calculated over the valid region
  //THIS MAY CHANGE
  if(_u[_qp] > 0.5)
    return 0.0;
  else if(_u[_qp] < 0.0)
    return 1e20;
  else
    return (_R*_T/(_u[_qp] - 3.*_u[_qp]*_u[_qp] + 2.*_u[_qp]*_u[_qp]*_u[_qp]))/_molarVol_alpha_Zr[_qp];
}

Real
CHCoupledCalphad::computeD2GdeltaDx2()
{
  Real ideal, excess;

  //do some checking so that equations are only calculated over the valid region
  //THIS MAY CHANGE
  if(_u[_qp] > 2./3.)
    return 1e20;
  else if ( _u[_qp] < 0.0)
    return 0.0;
  else
  {
    ideal = 2.*_R*_T/(2.*_u[_qp] - 5.*_u[_qp]*_u[_qp] + 3.*_u[_qp]*_u[_qp]*_u[_qp]);
    excess = (_L0[_qp]*(_u[_qp] - 1.) + _L1[_qp]*(3. - 6.*_u[_qp])) / (2.*pow(_u[_qp] - 1., 4.));

    return (ideal + excess)/_molarVol_delta_ZrH2[_qp];
  }
}

Real
CHCoupledCalphad::computeD3GalphaDx3()
{
  return (-1.*_R*_T*(1. - 6.*_u[_qp] + 6.*_u[_qp]*_u[_qp])/std::pow(_u[_qp] - 3.*_u[_qp]*_u[_qp]                                                    + 2.*_u[_qp]*_u[_qp]*_u[_qp], 2.))/_molarVol_alpha_Zr[_qp];
}

Real
CHCoupledCalphad::computeD3GdeltaDx3()
{
  Real ideal, excess;

  ideal = -2.*_R*_T*(2. - 10.*_u[_qp] + 9.*_u[_qp]*_u[_qp])/std::pow(2.*_u[_qp] - 5.*_u[_qp]*_u[_qp]
                                                      + 3.*_u[_qp]*_u[_qp]*_u[_qp], 2.);

  excess =  3.*(_L0[_qp]*(1.-_u[_qp]) + _L1[_qp]*(6.*_u[_qp] - 2.))/(2.*std::pow((_u[_qp] - 1.), 5.));

  return (ideal + excess)/_molarVol_delta_ZrH2[_qp];
}

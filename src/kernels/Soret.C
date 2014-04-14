/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  10 April 2014
*
*************************************************************************/

#include "Soret.h"
#include <iostream>

template<>
InputParameters validParams<Soret>()
{
  InputParameters params = validParams<Diffusion>();

  params.addRequiredCoupledVar("T", "temperature field variable");

  return params;
}


Soret::Soret(const std::string & name, InputParameters parameters) :
    Diffusion(name, parameters),
    _L1Q(getMaterialProperty<Real>("L1Q")),
    _T(coupledValue("T")),
    _grad_T(coupledGradient("T"))
{
}

Real
Soret::computeQpResidual()
{
  //  std::cout<<"soret_residual = "<< ( _grad_test[_i][_qp]*_L1Q[_qp]*_grad_T[_qp] ) / _T[_qp]<<std::endl;
  return ( _grad_test[_i][_qp]*_L1Q[_qp]*_grad_T[_qp] ) / _T[_qp];
}

Real
Soret::computeQpJacobian()
{
 return  _grad_phi[_j][_qp]*_L1Q[_qp]*_grad_T[_qp]/_T[_qp];
}


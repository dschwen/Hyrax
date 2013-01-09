/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  24 October 2012
*
*************************************************************************/

#include "CHBulkSimmons.h"

template<>
InputParameters validParams<CHBulkSimmons>()
{
  InputParameters params = validParams<CHBulkCoupled>();

  return params;
}

CHBulkSimmons::CHBulkSimmons(const std::string & name, InputParameters parameters)
  :CHBulkCoupled(name, parameters)
{
}

RealGradient
CHBulkSimmons::computeGradDFDCons(PFFunctionType type, Real /*c*/, RealGradient grad_c)
{
  switch (type)
  {
  case Residual:
    return _a1[_qp]*(grad_c) - _a2[_qp]*(_coupled_OP_var[_qp]*_coupled_OP_grad[_qp]) ;

  case Jacobian:
    return _a1[_qp]*_grad_phi[_j][_qp] ;
  }

  mooseError("Invalid type passed in");
}

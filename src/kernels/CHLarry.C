/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  2 October 2013
*
*************************************************************************/

#include "CHLarry.h"

template<>
InputParameters validParams<CHLarry>()
{
  InputParameters params = validParams<CHBulk>();
  params.addRequiredParam<Real>("W", "well height of the double-well functional");

  return params;
}

CHLarry::CHLarry(const std::string & name, InputParameters parameters)
    :CHBulk(name, parameters),
     _W(getParam<Real>("W"))
{
}

RealGradient
CHLarry::computeGradDFDCons(PFFunctionType type, Real c, RealGradient grad_c)
{
  switch (type)
  {
  case Residual:
    return _W*(2.0 - 12.0*c + 12.0*c*c)*grad_c;

  case Jacobian:
    //this isn't necessarily perfect
    return _W*_grad_phi[_j][_qp]*(-12.0*_phi[_j][_qp] + 24.0*c*_phi[_j][_qp]);
  }
  mooseError("invalid type passed in");
}

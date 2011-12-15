/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  14 December 2011
*
*  This code inherits from AuxKernel in MOOSE
*  
*  This code handles the supersaturation calculation in the concurrent
*  nucleation and growth algorithm first proposed by J.P. Simmons.
*  
*************************************************************************/

#include "AuxSupersaturation.h"

template<>
InputParameters validParams<AuxSupersaturation>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("coupled_var","coupled variable: concentration"); 
  // I'd have line in input file, coupled_var = concentration
  return params;
}

AuxSupersaturation::AuxSupersaturation(const std::string & name, InputParameters parameters) 
  : AuxKernel(name, parameters),
  _coupled_conc(coupledValue("coupled_var")),
  _c1(getMaterialProperty<Real>("C1"))
{
}

Real
AuxSupersaturation::computeValue()
{ 
  _supersaturation = _coupled_conc[_qp] - _c1[_qp];  
  if (_supersaturation <= 0.0)
  {
    _supersaturation = 1.0e-10;  // fixed to some arbitrary small value but preventing division by zero
  }
 return _supersaturation;
}


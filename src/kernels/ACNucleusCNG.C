/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  19 May 2012
*
*************************************************************************/

#include "ACNucleusCNG.h"

template<>
InputParameters validParams<ACNucleusCNG>()
{
  InputParameters params = validParams<ACNucleus>();
  params.addRequiredCoupledVar("nuc_variable", "Auxiliary variable for nucleation");

  return params;
}

ACNucleusCNG::ACNucleusCNG(const std::string & name, InputParameters parameters)
    : ACNucleus(name, parameters),
      _nucleation(coupledValue("nuc_variable"))
{
}

Real
ACNucleusCNG::computeQpResidual()
{
  if(_nucleation[_qp] > 1.0)
  {
    //do the nucleus introduction (modified from ACNucleus)
    return 0.0;
  }
  else
    return 0.0;
}



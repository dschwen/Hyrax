/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  22 December 2011
*
*************************************************************************/

#include "DiracForcedAMR.h"

#include "elem.h"

/**
 * DiracForcedAMR is designed to have the introduction of a point source after some amount of time
 * has passed.  This kernel is designed to work with the TransientMultiAMR executioner for multiple
 * refinement levels in one timestep.
 */

template<>
InputParameters validParams<DiracForcedAMR>()
{
 InputParameters params = validParams<ConstantPointSource>();
  params.addRequiredParam<int>("active_after","timestep number at which to turn on the Dirac kernel");
  // to simulate the Dirac function being activated at some random time for the CNG algorithm
 return params;
}

DiracForcedAMR::DiracForcedAMR(const std::string & name, InputParameters parameters) :
    ConstantPointSource(name, parameters),
    _active_after(getParam<int>("active_after"))
{
}

void
DiracForcedAMR::addPoints()
{  
  if (_t_step >= _active_after)
  {
    //get the element info for that point - for use with the refinement forcing if necessary.
    const Elem * dirac_elem = addPoint(_p);
  }
}

Real
DiracForcedAMR::computeQpResidual()
{
  return -_test[_i][_qp]*_value;
}

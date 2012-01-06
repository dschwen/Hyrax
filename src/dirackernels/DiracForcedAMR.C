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
  params.addRequiredParam<int>("active_for", "number of timesteps the kernel is active");
  // to simulate the Dirac function being activated at some random time for the CNG algorithm
 return params;
}

DiracForcedAMR::DiracForcedAMR(const std::string & name, InputParameters parameters) :
    ConstantPointSource(name, parameters),
    _active_after(getParam<int>("active_after")),
    _active_for(getParam<int>("active_for"))
{
}

void
DiracForcedAMR::addPoints()
{  
  if (_t_step >= _active_after) 
// if I want it to operate for however long, we could do a few statements here, where you get what the initial current time is, then if _t >= init_current_time+1.0 or whatever, call addPoints
  {
    if(_t_step >= _active_after+_active_for)
    {
      _value = _value/2.0;
    }
    //get the element info for that point - for use with the refinement forcing if necessary.
    const Elem * dirac_elem = addPoint(_p);
  }

//  // Do this several more times for more arbitrary points
//  if (_t_step >= 5)
//  {
//	addPoint(Point(1.2, 2.6, 0.0));
//  }

//  if (_t_step >= 7)
//  {
//	addPoint(Point(7.4, 3.3, 0.0));
//  }

}

Real
DiracForcedAMR::computeQpResidual()
{
  return -_test[_i][_qp]*_value;
}

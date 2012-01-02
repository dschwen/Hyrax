/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  20 December 2011
*  
*************************************************************************/

#include "TransientMultiAMR.h"

// C++ Includes
#include <iomanip>
#include <iostream>
#include <fstream>

/**
 * Transient executioner without adaptive timestepping that allows for multiple mesh adaptivity steps
 * within one timestep.  (Designed for use with Dirac kernels for modeling nucleation.)
 */

template<>
InputParameters validParams<TransientMultiAMR>()
{
  InputParameters params = validParams<Transient>();
  params.addRequiredParam<int>("num_refines","Number of mesh refine steps to perform at each timestep");
  // there will always be a variable in the input file called "num_refines"
  return params;
}

TransientMultiAMR::TransientMultiAMR(const std::string & name, InputParameters parameters) :
    Transient(name, parameters),
    _num_refines(getParam<int>("num_refines"))
{
}

void
TransientMultiAMR::endStep()
{
  // Most of this is lifted from Transient

  // if _reset_dt is true, force the output no matter what
  _problem.output(_reset_dt);
  _problem.outputPostprocessors(_reset_dt);

#ifdef LIBMESH_ENABLE_AMR
  if (_problem.adaptivity().isOn())
  {
   for (int i=1; i<= _num_refines; i++) // This is the loop for multiple adaptivity steps per timestep
   {
    _problem.adaptMesh();
    _problem.out().meshChanged();
   }
  }
#endif

  _t_step++;
  _time_old = _time;

}

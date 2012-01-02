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
 * This kernel is designed to work with the TransientMultiAMR executioner. The mesh refinement for the 
 * element in which the point at which the Dirac kernel was introduced needs to be forced several levels.
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
//  _my_refine(_mesh),  //part of the previous attempt at mesh refinement
    _active_after(getParam<int>("active_after"))
{
}

void
DiracForcedAMR::addPoints()
{  
  if (_t_step >= _active_after)
  {
    //get the element info for that point - for use with the refinement forcing
    const Elem * dirac_elem = addPoint(_p);
  }

/* This is what had been come up with, but fails.  
//    Elem * elem;
//    elem = const_cast<Elem *> (dirac_elem);
//    elem->set_refinement_flag(libMesh::Elem::REFINE);
//    _my_refine.refine_elements();  // This is fun - compiles; but breaks on runtime.
*/
}

Real
DiracForcedAMR::computeQpResidual()
{
  return -_test[_i][_qp]*_value;
}

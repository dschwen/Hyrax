/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  20 December 2011
*  
*************************************************************************/

#ifndef TRANSIENTMULTIAMR_H
#define TRANSIENTMULTIAMR_H

#include "Transient.h"

//Forward declarations
class TransientMultiAMR;

template<>
InputParameters validParams<TransientMultiAMR>();

/**
 * Transient executioner without adaptive timestepping that allows for multiple mesh adaptivity steps
 * within one timestep.  (Designed for use with Dirac kernels for modeling nucleation.)
 */
class TransientMultiAMR: public Transient
{
public: 
  TransientMultiAMR(const std::string & name, InputParameters parameters);

protected:
/**
 * endStep() is the only overridden function.  It's the same as Transient::endStep() with the addition 
 * of a loop to call the mesh adaptivity function multiple times.
 */
  virtual void endStep();
 
private:
  int _num_refines;

};

#endif //TRANSIENTMULTIAMR_H

/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  15 December 2011
*
*  This code inherits from DiracKernel in MOOSE
*  
*  This code handles the introduction of nuclei as delta functions for 
*  the concurrent nucleation and growth algorithm first proposed by 
*  J.P. Simmons.
*  
*************************************************************************/

#ifndef DIRACNUCLEATION_H
#define DIRACNUCLEATION_H

// Moose Includes
#include "DiracKernel.h"

/* Remember, in the input file, tell this Dirac Kernel that it's operating on the appropriate
   variable - in this case, on the order parameter variable.  When I've got multiple OPs, this is
   going to get interesting. */

//Forward Declarations
class DiracNucleation;

template<>
InputParameters validParams<DiracNucleation>();

class DiracNucleation : public DiracKernel
{
public:
  DiracNucleation(const std::string & name, InputParameters parameters);
  virtual void addPoints();
  virtual Real computeQpResidual();

protected:
  Point _p;
  Real _value;

private: 
  VariableValue & _coupled_nucleation;

  Elem _true_element;
  Point _true_point;

// probably some other stuff here.
};


#endif //DIRACNUCLEATION_H

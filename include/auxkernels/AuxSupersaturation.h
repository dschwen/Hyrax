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

#ifndef AUXSUPERSATURATION_H
#define AUXSUPERSATURATION_H

#include "AuxKernel.h"

//forward declaration
class AuxSupersaturation;

template<>
InputParameters validParams<AuxSupersaturation>();

class AuxSupersaturation : public AuxKernel
{
public:
  AuxSupersaturation(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  // probably some other things here
  // Concentration coupled in

  
private:

  // may or may not have anything here
  VariableValue & _coupled_conc;
  MaterialProperty<Real> & _c1;
  Real _supersaturation;
};

// Supersaturation = C - C1


#endif //AUXSUPERSATURATION_H

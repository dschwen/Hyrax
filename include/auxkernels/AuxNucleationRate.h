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
*  This code handles the nucleation rate calculation in the concurrent
*  nucleation and growth algorithm first proposed by J.P. Simmons.
*  
*************************************************************************/

#ifndef AUXNUCLEATIONRATE_H
#define AUXNUCLEATIONRATE_H

#include "AuxKernel.h"

//forward declaration
class AuxNucleationRate; //j_star

template<>
InputParameters validParams<AuxNucleationRate>();

class AuxNucleationRate : public AuxKernel
{
public:
  AuxNucleationRate(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  // probably some other things here

  
private:

  // may or may not have anything here

  VariableValue & _coupled_supersaturation;
  Real _Kn1;
  Real _Kn2;

};

// j_star = Kn1 * exp(-1*Kn2 / supersaturation) 

#endif //AUXNUCLEATIONRATE_H

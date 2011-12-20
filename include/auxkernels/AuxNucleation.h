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
*  This code handles the nucleation/no nucleation portion of the concurrent
*  nucleation and growth algorithm first proposed by J.P. Simmons.
*  
*************************************************************************/

#ifndef AUXNUCLEATION_H
#define AUXNUCLEATION_H

#include "AuxKernel.h"

//forward declaration
class AuxNucleation;

template<>
InputParameters validParams<AuxNucleation>();

/**
 * returns a boolean: true if nucleation occured; false if it didn't.
 */
class AuxNucleation : public AuxKernel
{
public:
  AuxNucleation(const std::string & name, InputParameters params);

protected: 
  // this might actually be a bool.  but i don't think I can make it a bool return type.
  virtual Real computeValue();

  // probably some other things here, or maybe nothing
  
private:

  // may or may not have anything here

  VariableValue & _coupled_probability;
  double _random_number;
  
};

#endif //AUXNUCLEATION_H

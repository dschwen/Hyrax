/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  19 May 2012
*
*************************************************************************/

#ifndef ACNUCLEUSCNG_H
#define ACNUCLEUSCNG_H

#include "ACNucleus.h"

class ACNucleusCNG;

template<>
InputParameters validParams<ACNucleusCNG>();

class ACNucleusCNG : public ACNucleus
{
public:
  ACNucleusCNG(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeQpResidual();

private:
  //aux variable for the nucleation true/false field.
  VariableValue & _nucleation;


};

#endif //ACNUCLEUSCNG_H

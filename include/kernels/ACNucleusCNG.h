/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  29 May 2012
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
  MaterialProperty<std::vector<RealVectorValue> > & _nucleation_locations;
  MaterialProperty<std::vector<Real> > & _start_times;
  MaterialProperty<std::vector<Real> > & _end_times;
};

#endif //ACNUCLEUSCNG_H

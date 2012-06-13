/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  12 June 2012
*
*************************************************************************/

#ifndef ACINTERFACENUCLEATION_H
#define ACINTERFACENUCLEATION_H

#include "ACInterface.h"

//forward declarations
class ACInterfaceNucleation;

template<>
InputParameters validParams<ACInterfaceNucleation>();

class ACInterfaceNucleation : public ACInterface
{
public:

  ACInterfaceNucleation(const std::string & name, InputParameters parameters);

protected:

  virtual RealGradient precomputeQpResidual();
  virtual RealGradient precomputeQpJacobian();

private:

  Real _start_time;
  Real _end_time;
  Real _radius;
//  Real _int_width;
  Point _nucleation_center;
};

#endif //ACINTERFACENUCLEATION_H

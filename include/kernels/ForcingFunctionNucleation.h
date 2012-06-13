/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  12 June 2012
*
*************************************************************************/

#ifndef FORCINGFUNCTIONNUCLEATION_H
#define FORCINGFUNCTIONNUCLEATION_H

#include "UserForcingFunction.h"

// forward declarations
class ForcingFunctionNucleation;

template<>
InputParameters validParams<ForcingFunctionNucleation>();

class ForcingFunctionNucleation : public UserForcingFunction
{
public:

  ForcingFunctionNucleation(const std::string & name, InputParameters parameters);

protected:

  virtual Real computeQpResidual();

private:

  Real _start_time;
  Real _end_time;
  Real _radius;
//  Real _int_width;
  Point _nucleation_center;
};

#endif //FORCINGFUNCTIONNUCLEATION_H

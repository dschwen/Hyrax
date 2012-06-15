/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  12 June 2012
*
*************************************************************************/

#ifndef VALUENUCLEATION_H
#define VALUENUCLEATION_H

#include "Value.h"

//forward declarations
class ValueNucleation;

/**
 * ValueNucleation is for creating a round nucleus of some value with a
 * diffuse interface. To be used with ForcingFunctionNucleation to perform
 * the L2 projection.
 */

template<>
InputParameters validParams<ValueNucleation>();

class ValueNucleation : public Value
{
public:

  ValueNucleation(const std::string & name, InputParameters parameters);

protected:

  virtual Real computeQpResidual();

private:

  Real _start_time;
  Real _end_time;
  Real _radius;
  //Real _int_width;
  Point _nucleation_center;
};

#endif //VALUENUCLEATION_H

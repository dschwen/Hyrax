/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  14 June 2012
*
*************************************************************************/

#ifndef TIMEDERIVATIVENUCLEATION_H
#define TIMEDERIVATIVENUCLEATION_H

#include "TimeDerivative.h"

//forward declarations
class TimeDerivativeNucleation;

template<>
InputParameters validParams<TimeDerivativeNucleation>();

class TimeDerivativeNucleation : public TimeDerivative
{
public:

  TimeDerivativeNucleation(const std::string & name, InputParameters parameters);

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:

  Real _start_time;
  Real _end_time;
  Real _radius;
  //Real _int_width;
  Point _nucleation_center;
};

#endif //TIMEDERIVATIVENUCLEATION_H

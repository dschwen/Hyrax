/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  17 May 2012
*
*************************************************************************/

#ifndef ACNUCLEUS_H
#define ACNUCLEUS_H

#include "Kernel.h"

class ACNucleus;

template<>
InputParameters validParams<ACNucleus>();

class ACNucleus : public Kernel
{
public:

  ACNucleus(const std::string & name, InputParameters parameters);

protected:

  virtual Real computeQpResidual();

private:

  Real _n_value;
  Real _radius;
  Real _int_width;
  Real _x_center;
  Real _y_center;
  Real _z_center;
  Real _start_time;
  Real _end_time;

};

#endif //ACNUCLEUS_H

/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  1 November 2013
*
*************************************************************************/

#ifndef HEAT_H
#define HEAT_H

#include "Diffusion.h"

//forward declaration
class Heat;

template<>
InputParameters validParams<Heat>();

class Heat : public Diffusion
{
public:

  Heat(const std::string & name, InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  MaterialProperty<Real> & _diffusivity;
  MaterialProperty<Real> & _dDiffusivity_dT;

  //anything else here that you may want to add, like coupled variables

private:

};

#endif //HEAT_H

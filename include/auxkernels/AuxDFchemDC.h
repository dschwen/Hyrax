/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
* 25 April 2015
*
*************************************************************************/

#ifndef AUXDFCHEMDC_H
#define AUXDFCHEMDC_H

#include "AuxKernel.h"

class AuxDFchemDC;

template<>
InputParameters validParams<AuxDFchemDC>();

class AuxDFchemDC : public AuxKernel
{
public:
    AuxDFchemDC(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

private:

  VariableValue & _X;
  VariableValue & _OP;

  MaterialProperty<Real> & _Omega;
  MaterialProperty<Real> & _dGalpha_dc;
  MaterialProperty<Real> & _dGdelta_dc;

//Real _scaling_factor;

};

#endif //AUXELASTICENERGY_H

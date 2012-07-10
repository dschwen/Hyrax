/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  7 July 2012
*
*************************************************************************/

#ifndef AUXTESTFLIP_H
#define AUXTESTFLIP_H

#include "AuxKernel.h"

class AuxTestFlip;

template<>
InputParameters validParams<AuxTestFlip>();

class AuxTestFlip : public AuxKernel
{

public:
  AuxTestFlip(const std::string & name, InputParameters parameters);

protected:

  virtual Real computeValue();

private:

  VariableValue & _coupled_c;
  Real _radius;
  Point _center_point;

};

#endif //AUXTESTFLIP_H

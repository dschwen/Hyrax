/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  2 October 2013
*
*************************************************************************/

#ifndef CHLARRY_H
#define CHLARRY_H

#include "CHBulk.h"

class CHLarry;

template<>
InputParameters validParams<CHLarry>();

class CHLarry : public CHBulk
{
public:
  CHLarry(const std::string & name, InputParameters parameters);

protected:
  virtual RealGradient computeGradDFDCons(PFFunctionType type, Real c, RealGradient grad_c);

private:

  Real _W;
};

#endif //CHLARRY_H
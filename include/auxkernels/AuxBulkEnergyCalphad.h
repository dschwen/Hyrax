/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  9 January 2014
*
*************************************************************************/

#ifndef AUXBULKENERGYCALPHAD_H
#define AUXBULKENERGYCALPHAD_H

#include "AuxKernel.h"
#include "CalphadAB1CD1.h"
#include "CalphadAB1CD2.h"

class AuxBulkEnergyCalphad;

template<>
InputParameters validParams<AuxBulkEnergyCalphad>();

class AuxBulkEnergyCalphad : public AuxKernel
{
public:
  AuxBulkEnergyCalphad(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  Real computeHeaviside();
  Real computeBarrier();

private:

  VariableValue & _C;
  std::vector<VariableValue *> _OP;
  unsigned int _n_OP_variables;

  Real _scaling_factor;

  MaterialProperty<Real> & _W;                            //Well height
  MaterialProperty<Real> & _Omega;                        //Molar volume
  MaterialProperty<Real> & _G_alpha;                      //Gmix_alpha
  MaterialProperty<Real> & _G_delta;                      //Gmix_delta

};

#endif //AUXBULKENERGYCALPHAD_H

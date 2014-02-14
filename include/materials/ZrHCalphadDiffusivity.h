/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  13 February 2014
*
*************************************************************************/

#ifndef ZRHCALPHADDIFFUSIVITY_H
#define ZRHCALPHADDIFFUSIVITY_H

#include "ZrHCalphad.h"

//forward declaration
class ZrHCalphadDiffusivity;

template<>
InputParameters validParams<ZrHCalphad>();

class ZrHCalphadDiffusivity : public ZrHCalphad
{
public:
  ZrHCalphadDiffusivity(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();
  Real computeHeaviside();


private:
  //Diffusion coefficient information
  Real _H_Zr_D0;
  Real _H_ZrH2_D0;
  Real _H_Zr_Q0;
  Real _H_ZrH2_Q0;
  Real _R;

  MaterialProperty<Real> & _d2Galpha_dc2;
  MaterialProperty<Real> & _d2Gdelta_dc2;

  MaterialProperty<Real> & _D_alpha;
  MaterialProperty<Real> & _D_delta;

  unsigned int _n_OP_variables;
  std::vector<VariableValue *> _OP;

  VariableValue & _c;
  VariableValue & _T;

};

#endif //ZRHCALPHADDIFFUSIVITY_H

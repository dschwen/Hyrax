/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*
*  19 November 2013
*
*************************************************************************/
#ifndef CHCOUPLEDCALPHADSPLIT_H
#define CHCOUPLEDCALPHADSPLIT_H

#include "SplitCHCRes.h"

class CHCoupledCalphadSplit;

template<>
InputParameters validParams<CHCoupledCalphadSplit>();

class CHCoupledCalphadSplit : public SplitCHCRes
{
public:
  CHCoupledCalphadSplit(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeDFDC(PFFunctionType type);

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  Real computeHeaviside();
  Real computeDHeaviside();

  // RealGradient computeGradConservedTerm();
  //RealGradient computeGradNonconservedTerm();

private:
  MaterialProperty<Real> & _W;                            //Well height
  MaterialProperty<Real> & _Omega;                        //Molar volume
  MaterialProperty<Real> & _dGalpha_dc;
  MaterialProperty<Real> & _d2Galpha_dc2;
  //MaterialProperty<Real> & _d3Galpha_dc3;
  MaterialProperty<Real> & _dGdelta_dc;
  MaterialProperty<Real> & _d2Gdelta_dc2;
  //MaterialProperty<Real> & _d3Gdelta_dc3;

  VariableValue & _OP;

  unsigned int _OP_var;

  //unsigned int _n_OP_variables;
  //std::vector<VariableValue *> _OP;
  //std::vector<VariableGradient *> _grad_OP;

//  Real _Heaviside;
//  Real _dHeaviside;
};

#endif //CHCOUPLEDCALPHADSPLIT_H

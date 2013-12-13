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
  Real computeDHeaviside(unsigned int i);

private:
  MaterialProperty<Real> & _W;                            //Well height
  MaterialProperty<Real> & _Omega;                        //Molar volume

  MaterialProperty<Real> & _dGalpha_dc;
  MaterialProperty<Real> & _d2Galpha_dc2;

  MaterialProperty<Real> & _dGdelta_dc;
  MaterialProperty<Real> & _d2Gdelta_dc2;



  //VariableValue & _OP;

  //unsigned int _OP_var;

  Real _scaling_factor;

  //unsigned int _n_OP_variables;
  //std::vector<VariableValue *> _OP;
  //std::vector<VariableGradient *> _grad_OP;

  // unsigned int _n_var;
  unsigned int _w_var;
  unsigned int _T_var;

  //VariableValue & _n;
  VariableValue & _w;
  VariableValue & _T;

  unsigned int _n_OP_vars;
  std::vector<unsigned int> _n_var;
  std::vector<VariableValue *> _OP;
};

#endif //CHCOUPLEDCALPHADSPLIT_H

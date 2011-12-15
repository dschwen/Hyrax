/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  16 November 2011
*
*  This code inherits from CHBulk in ELK
*  
*  This code handles the Cahn-Hilliard equation (conserved order parameter)
*  evolution coupled to a non-conserved order parameter.
*  
*************************************************************************/

#ifndef CHBulkCoupled_H
#define CHBulkCoupled_H

#include "CHBulk.h"

/* This class is for the concentration of hydrogen in the cladding, 
evolved using the Cahn-Hilliard equation.  It couples to an order 
parameter from the Alan-Cahn equation*/

//Forward Declarations
class CHBulkCoupled;

template<>
InputParameters validParams<CHBulkCoupled>();

class CHBulkCoupled : public CHBulk
{
public:

  CHBulkCoupled(const std::string & name, InputParameters parameters);
  
protected:

  std::string _a1_name;
  std::string _a2_name;
  std::string _c1_name;
 
  virtual RealGradient computeGradDFDCons(PFFunctionType type, Real c, RealGradient grad_c);
  
private:
  
  MaterialProperty<Real> & _a1;
  MaterialProperty<Real> & _a2;
  MaterialProperty<Real> & _c1;

  VariableValue & _coupled_var;
  VariableGradient & _coupled_grad;
};
#endif //CHBulkCoupled_H

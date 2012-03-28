/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  16 November 2011
*  
*************************************************************************/

#ifndef CHBULKCOUPLED_H
#define CHBULKCOUPLED_H

#include "CHBulk.h"

/** CHBulkCoupled handles the conserved order parameter(probably concentration), 
 * evolved using the Cahn-Hilliard equation.  It couples to an order 
 * parameter from the Alan-Cahn equation.
 */

//Forward Declarations
class CHBulkCoupled;

template<>
InputParameters validParams<CHBulkCoupled>();

class CHBulkCoupled : public CHBulk
{
public:

  CHBulkCoupled(const std::string & name, InputParameters parameters);
  
protected:
 
  /**
   * computeGradDFDCons()
   * @return returns the GRADIENT of the partial(bulk free energy)/partial(c).  Don't screw that up.
   */

  virtual RealGradient computeGradDFDCons(PFFunctionType type, Real c, RealGradient grad_c);
  
private:
  
  MaterialProperty<Real> & _a1;  ///< Landau polynomial parameters (see Guo, 2008)
  MaterialProperty<Real> & _a2;
  MaterialProperty<Real> & _c1;  ///< position-ish of 1st energy well in c-space (terminal solid solubility)

  VariableValue & _coupled_OP_var;  ///< Allen-Cahn equation variable (order parameter, probably)
  VariableGradient & _coupled_OP_grad;  ///< gradient of AC variable
};

#endif //CHBULKCOUPLED_H

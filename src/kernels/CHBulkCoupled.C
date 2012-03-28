/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  16 November 2011
*
*************************************************************************/

#include "CHBulkCoupled.h"

/** CHBulkCoupled handles the conserved order parameter(probably concentration), 
 * evolved using the Cahn-Hilliard equation.  It couples to an order 
 * parameter from the Alan-Cahn equation.
 */

template<>
InputParameters validParams<CHBulkCoupled>()
{
  InputParameters params = validParams<CHBulk>();
  params.addRequiredCoupledVar("coupled_OP_var", "The order parameter coupled to the CH eqn");
//params.addRequiredCoupledVar("coupled_var","The order parameter to be coupled to the CH equation");
  // this is a little extraneous and can get cleaned up
//  params.addParam<std::string>("first_landau","A1","First Landau coefficient"); 
//  params.addParam<std::string>("second_landau","A2","Second Landau coefficient"); 
//  params.addParam<std::string>("first_well","C1","First Landau well position"); 

  return params;
}

CHBulkCoupled::CHBulkCoupled(const std::string & name, InputParameters parameters)
  :CHBulk(name, parameters),

   // Pull in the names from the input file for the coefficients and then declare the materials properties
   // _a1_name(getParam<std::string>("first_landau")),
   //_a2_name(getParam<std::string>("second_landau")),
   //_c1_name(getParam<std::string>("first_well")),
   //
   // Get the material properties for the Landau coefficients
   _a1(getMaterialProperty<Real>("A1")),
   _a2(getMaterialProperty<Real>("A2")),
   _c1(getMaterialProperty<Real>("C1")),
   // Make the coupled value whatever is directed in the input file
   _coupled_OP_var(coupledValue("coupled_OP_var")),
   _coupled_OP_grad(coupledGradient("coupled_OP_var"))
{
}

RealGradient
CHBulkCoupled::computeGradDFDCons(PFFunctionType type, Real c, RealGradient grad_c)
{
  switch (type)
  {
  case Residual:
    return _a1[_qp]*(grad_c) + _a2[_qp]*(_coupled_OP_var[_qp]*_coupled_OP_grad[_qp]) ;
 
  case Jacobian: 
    return _a1[_qp]*_grad_phi[_j][_qp] ;
  }
  
  mooseError("Invalid type passed in");
}

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
*  This code handles the materials parameters for a coupled 
*  conserved order parameter, non-conserved order parameter
*  system. 
*************************************************************************/

#include "ACBulkCoupled.h"

template<>
InputParameters validParams<ACBulkCoupled>()
{
  InputParameters params = validParams<ACBulk>();
  params.addRequiredCoupledVar("coupled_var","The concentration to be coupled to the AC equation");
  params.addParam<std::string>("second_landau","A2","Second Landau coefficient"); 
  params.addParam<std::string>("third_landau","A3","Third Landau coefficient"); 
  params.addParam<std::string>("fourth_landau","A4","Fourth Landau coefficient"); 
  params.addParam<std::string>("second_well","C2","Second Landau well position"); 
  return params;
}

ACBulkCoupled::ACBulkCoupled(const std::string & name, InputParameters parameters)
  :ACBulk(name, parameters),

   // Pull in the names from the input file for the coefficients and then declare the materials properties
   _a2_name(getParam<std::string>("second_landau")),
   _a3_name(getParam<std::string>("third_landau")),
   _a4_name(getParam<std::string>("fourth_landau")),
   _c2_name(getParam<std::string>("second_well")),

   // Get the material property for the name that's been given
   _a2(getMaterialProperty<Real>(_a2_name)),
   _a3(getMaterialProperty<Real>(_a3_name)),
   _a4(getMaterialProperty<Real>(_a4_name)),
   _c2(getMaterialProperty<Real>(_c2_name)),

   // Make the coupled value whatever is directed in the input file
   _coupled_var(coupledValue("coupled_var"))

{}

Real
ACBulkCoupled::computeDFDOP(PFFunctionType type)
{
  switch (type)
  {
  case Residual:
    return _a2[_qp]*(_coupled_var[_qp]- _c2[_qp])*_u[_qp] - _a3[_qp]*_u[_qp]*_u[_qp]*_u[_qp] 
           + _a4[_qp]*_u[_qp]*_u[_qp]*_u[_qp]*_u[_qp]*_u[_qp] ;

  case Jacobian:
    return _phi[_j][_qp]*( _a2[_qp]*(_coupled_var[_qp]- _c2[_qp]) - 3.0*_a3[_qp]*_u[_qp]*_u[_qp] + 
           5.0*_a4[_qp]*_u[_qp]*_u[_qp]*_u[_qp]*_u[_qp] ) ;
  }

  mooseError("Invalid type passed in");
}

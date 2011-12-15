/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  14 December 2011
*
*  This code inherits from AuxKernel in MOOSE
*  
*  This code handles the nucleation probability calculation in the concurrent
*  nucleation and growth algorithm first proposed by J.P. Simmons.
*  
*************************************************************************/

#include "AuxNucleationProbability.h"

template<>
InputParameters validParams<AuxNucleationProbability>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("coupled_aux_var","coupled auxiliary variable: nucleation rate"); 
  // I'd have line in input file, coupled_aux_var = j_star
  return params;
}

//  p_nm = 1 - exp(-1*j_star*dt) 

AuxNucleationProbability::AuxNucleationProbability(const std::string & name, InputParameters parameters) 
  : AuxKernel(name, parameters),
  _coupled_nuc_rate(coupledValue("coupled_aux_var"))
{
}

Real
AuxNucleationProbability::computeValue()
{ 
 return 1.0 - exp(-1.0*_coupled_nuc_rate[_qp]*_dt);
}


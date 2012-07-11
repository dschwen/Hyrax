/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*************************************************************************/

#include "AuxNucleationProbability.h"

#include <ostream>

/**
 *  AuxNucleationProbability handles the nucleation probability (P_nm) calculation in the concurrent
 *  nucleation and growth algorithm first proposed by J.P. Simmons (2000).
 *  Returns the nucleation probability over the domain.
 */

template<>
InputParameters validParams<AuxNucleationProbability>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("coupled_aux_var","coupled auxiliary variable: nucleation rate");
  // I'd have line in input file, coupled_aux_var = j_star
  return params;
}

AuxNucleationProbability::AuxNucleationProbability(const std::string & name, InputParameters parameters)
  : AuxKernel(name, parameters),
  _coupled_nuc_rate(coupledValue("coupled_aux_var"))
{
}

Real
AuxNucleationProbability::computeValue()
{
 //  p_nm = 1 - exp(-1*j_star*dt)
 //  this maybe should just be for each element
  std::cout<<"probability = "<<1.0 - exp(-1.0*_coupled_nuc_rate[_qp]*_dt)<<std::endl;
 return 1.0 - exp(-1.0*_coupled_nuc_rate[_qp]*_dt);
}


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
  params.addRequiredCoupledVar("coupled_variables", "coupled order parameter variables");
  params.addRequiredParam<int>("n_OP_vars", "# of coupled OP variables");

  return params;
}

AuxNucleationProbability::AuxNucleationProbability(const std::string & name, InputParameters parameters)
  : AuxKernel(name, parameters),
    _coupled_nuc_rate(coupledValue("coupled_aux_var")),
    _n_OP_vars(getParam<int>("n_OP_vars"))
    //  _coupled_OP(coupledValue("coupled_variable"))
{
  if(_n_OP_vars != coupledComponents("coupled_variables"))
    mooseError("Please match the number of orientation variants to coupled OPs (AuxNucleationProbability).");

  _coupled_OP.resize(_n_OP_vars);

  for(unsigned int i=0; i<_n_OP_vars; i++)
    _coupled_OP[i] = &coupledValue("coupled_variables", i);
}

Real
AuxNucleationProbability::computeValue()
{
  for(unsigned int i=0; i<_n_OP_vars; i++)
  {
    if((*_coupled_OP[i])[_qp] > 0.1)
      return 0.0;
  }

  return 1.0 - exp(-1.0*_coupled_nuc_rate[_qp]*_dt);
}


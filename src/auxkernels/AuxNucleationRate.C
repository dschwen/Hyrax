/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*************************************************************************/

#include "AuxNucleationRate.h"

/**
 *  AuxNucleationRate handles the nucleation rate (j_star) calculation in the concurrent
 *  nucleation and growth algorithm first proposed by J.P. Simmons (2000).
 */

template<>
InputParameters validParams<AuxNucleationRate>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<Real>("Kn1", "First nucleation rate coefficient");
  params.addRequiredParam<Real>("Kn2", "Second nucleation rate coefficient");
  params.addRequiredCoupledVar("coupled_aux_var","coupled auxiliary variable: supersaturation");
  // I'd have line in input file, coupled_aux_var = supersaturation
  return params;
}

AuxNucleationRate::AuxNucleationRate(const std::string & name, InputParameters parameters)
  : AuxKernel(name, parameters),
  _coupled_supersaturation(coupledValue("coupled_aux_var")),
  _Kn1(getParam<Real>("Kn1")),
  _Kn2(getParam<Real>("Kn2"))
{
}

Real
AuxNucleationRate::computeValue()
{
 // j_star = Kn1 * exp(-1*Kn2 / supersaturation)

 return _Kn1*exp(-1.0*_Kn2/_coupled_supersaturation[_qp]);
}


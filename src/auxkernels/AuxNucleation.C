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
*  This code handles the nucleation/no nucleation portion of the concurrent
*  nucleation and growth algorithm first proposed by J.P. Simmons.
*
*************************************************************************/

#include "AuxNucleation.h"

template<>
InputParameters validParams<AuxNucleation>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("coupled_aux_var","Auxiliary coupled variable: nucleation probability");
  // I'd have line in input file, coupled_aux_var = p_nm, the nucleation probability
  return params;
}

AuxNucleation::AuxNucleation(const std::string & name, InputParameters parameters)
  : AuxKernel(name, parameters),
  _coupled_probability(coupledValue("coupled_aux_var"))
{
}

Real
AuxNucleation::computeValue()
{
//  //supersaturation
//  supersaturation = coupled_val[_qp] - C1;  // C1 a material value.  Should couple that in.
//  if supersaturation <= 0,
//   supersaturation = 1x10-10;  // fixed to some arbitrary small value but preventing division by zero

//  // nucleation rate equation - this  can be arbitrary
//  j_star = Kn1 * exp(-1*Kn2 / supersaturation)  // Kn1 and Kn2 are going to be supplied from the input file

//  // probability of nucleation in that location for this timestep
//  p_nm = 1 - exp(-1*j_star*dt)  // the timestep is going to have to be supplied from MOOSE

  // currently will use rand(), the intrinsic random number function.  need a better one, though.

  // CJP: You might try Moose::seed(unsigned int) and Moose::rand()

  _random_number = double (rand()%100000);
  _random_number = _random_number/100000;


  if (_random_number < _coupled_probability[_qp])
  {
   // return true
   return 2.0;
  }
  else
  {
   // return false
   return 0.0;
  }
}

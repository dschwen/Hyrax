/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*************************************************************************/

#include "AuxNucleation.h"
#include "MooseRandom.h"

/**
 *  AuxNucleation handles the nucleation/no nucleation b of the concurrent
 *  nucleation and growth algorithm first proposed by J.P. Simmons (2000).
 *  Returns a sort-of boolean: true if nucleation occured; false if it didn't.
 */

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
  /**
   * _random_number_seed will start at one and be incremented by one each timestep.
   */
  _random_number_seed = 1;
}

Real
AuxNucleation::computeValue()
{
/* AuxNucleation is the final step in the probabilistic nucleation method.  The steps to get here are
 * shown below:
   1) supersaturation
      supersaturation = C - C1;  // C1 is from Guo (2008) Landau polynomial.
      if supersaturation <= 0, supersaturation = 1x10-10; this prevents division by zero or negative in step 2
   2) nucleation rate equation - this  can be arbitrary
      currently: j_star = Kn1 * exp(-1*Kn2 / supersaturation)
      Kn1 and Kn2 are going to be supplied from the input file
   3) probability of nucleation in that location for this timestep
      p_nm = 1 - exp(-1*j_star*dt)
   4) stochastic testing of nucleation: nucleation probability vs random number between 0 and 1
*/

   // we are controlling the random number seeding this way for reproducibility
   _random_number_seed += 1;
   MooseRandom::seed(_random_number_seed);
   _random_number = MooseRandom::rand();

   //std::cout << "random number generated = " << _random_number << std::endl;
    _random_number = _random_number/5000;

   std::cout << "actual random number used = " << _random_number << std::endl;
   std::cout << "coupled probability[" << _qp << "] =" << _coupled_probability[_qp] <<std::endl;

   if (_random_number < _coupled_probability[_qp])
   {
     // return true
     return 2.0;
   }
   else
   {
     // return false
     return -2.0;
   }
}

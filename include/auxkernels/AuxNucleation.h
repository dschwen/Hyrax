/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*   
*  15 December 2011
*
*************************************************************************/

#ifndef AUXNUCLEATION_H
#define AUXNUCLEATION_H

#include "AuxKernel.h"

class AuxNucleation;

template<>
InputParameters validParams<AuxNucleation>();

/**
 *  AuxNucleation handles the nucleation/no nucleation portion of the concurrent
 *  nucleation and growth algorithm first proposed by J.P. Simmons (2000).
 *  Returns a sort-of boolean: true if nucleation occured; false if it didn't.
 */

class AuxNucleation : public AuxKernel
{
public:
  AuxNucleation(const std::string & name, InputParameters params);

protected: 

  /**
   * computeValue()
   * @return this is a sort-of boolean; sadly can't actually be a boolean.  0.0 if false, 
   * 2.0 if true for a nucleation event occuring (element average value).
   */
  
  virtual Real computeValue();
  
private:

  VariableValue & _coupled_probability; ///< nucleation probability
  double _random_number;		///< for stochastic nucleation
  unsigned int _random_number_seed;  	///< to seed the random number generator

};

#endif //AUXNUCLEATION_H

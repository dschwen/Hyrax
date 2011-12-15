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

#ifndef AUXNUCLEATIONPROBABILITY_H
#define AUXNUCLEATIONPROBABILITY_H

#include "AuxKernel.h"

//forward declaration
class AuxNucleationProbability;  //P_nm

template<>
InputParameters validParams<AuxNucleationProbability>();

class AuxNucleationProbability : public AuxKernel
{
public:
  AuxNucleationProbability(const std::string & name, InputParameters params);

protected:
  virtual Real computeValue();

private:
  VariableValue & _coupled_nuc_rate;
};

//  p_nm = 1 - exp(-1*j_star*dt) 

#endif //AUXNUCLEATIONPROBABILITY_H

/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  6 December 2011
*
*  This code inherits from SymmElasticityTensor in ELK
*  
*  This code handles the materials parameters for an HCP stiffness tensor.
*
*************************************************************************/


#ifndef SYMMHCPELASTICITYTENSOR_H
#define SYMMHCPELASTICITYTENSOR_H

#include "SymmElasticityTensor.h"

/* 
Defines and HCP stiffness tensor.  The input requires C1111, C1122, C1133, C3333, C2323, and C1212. 
* 
*Note that by default this tensor is constant, meaning it doesn't change once it's been initialized.
*To modify this, pass in "false" to the constructor.
*/

//Forward declaration
class SymmHCPElasticityTensor

template<>
InputParameters validParams<SymmHCPElasticityTensor>();

class SymmHCPElasticityTensor : public SymmElasticityTensor
{
public:
  SymmHCPElasticityTensor(const bool constant = true, InputParameters parameters);

protected:

  virtual void calculateEntries(unsigned int qp);

private: 

}

#endif //SYMMHCPELASTICITYTENSOR_H

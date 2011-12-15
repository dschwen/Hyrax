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

#include "SymmHCPElasticityTensor.h"

SymmHCPElasticityTensor::SymmHCPElasticityTensor(const bool constant)
   : SymmElasticityTensor(constant)

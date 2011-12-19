/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  16 November 2011
*
*  
*  This code handles the materials parameters for a coupled 
*  conserved order parameter, non-conserved order parameter
*  system. 
*************************************************************************/


#include "Moose.h"
#include "Factory.h"
#include "ActionFactory.h"

//Elk Includes
#include "Elk.h"

//Kernels
#include "ACBulkCoupled.h"
#include "CHBulkCoupled.h"

//Auxiliary Kernels
#include "AuxNucleation.h"
#include "AuxNucleationProbability.h"
#include "AuxNucleationRate.h"
#include "AuxSupersaturation.h"

//Dirac Kernels

//Boundary Conditions

//Materials
#include "PFMobilityLandau.h"

//Initial Conditions

//Dampers

//Executioners

//Post Processors

//Actions

namespace Hyrax
{
  void registerObjects()
  {
    Moose::registerObjects();
    Elk::registerObjects();
    //Kernels
    registerKernel(CHBulkCoupled);
    registerKernel(ACBulkCoupled);

    //Auxiliary Kernels

    //Dirac Kernels

    //Boundary Conditions

    //Materials
    registerMaterial(PFMobilityLandau);    

    //Initial Conditions

    //Dampers

    //Executioners
    
    //Post Processors

    // Actions

  }
}

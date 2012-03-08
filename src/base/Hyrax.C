/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  16 November 2011
*
*************************************************************************/


#include "Moose.h"
#include "Factory.h"
#include "ActionFactory.h"

//Elk Includes
#include "Elk.h"

//Kernels
#include "ACBulkCoupled.h"
#include "CHBulkCoupled.h"
#include "ACTransformElasticDF.h"

//Auxiliary Kernels
#include "AuxNucleation.h"
#include "AuxNucleationProbability.h"
#include "AuxNucleationRate.h"
#include "AuxSupersaturation.h"

//Dirac Kernels
#include "DiracNucleation.h"
#include "DiracForcedAMR.h"

//Boundary Conditions

//Materials
#include "PFMobilityLandau.h"
#include "LinearSingleCrystalPrecipitateMaterial.h"

//Initial Conditions
#include "PolySpecifiedSmoothCircleIC.h"

//Dampers

//Executioners
#include "TransientMultiAMR.h"

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
    registerKernel(ACTransformElasticDF);

    //Auxiliary Kernels
    registerAux(AuxNucleation);
    registerAux(AuxSupersaturation);
    registerAux(AuxNucleationRate);
    registerAux(AuxNucleationProbability);

    //Dirac Kernels
    registerDiracKernel(DiracNucleation);
    registerDiracKernel(DiracForcedAMR);

    //Boundary Conditions

    //Materials
    registerMaterial(PFMobilityLandau);
    registerMaterial(LinearSingleCrystalPrecipitateMaterial);

    //Initial Conditions
    registerInitialCondition(PolySpecifiedSmoothCircleIC);

    //Dampers

    //Executioners
    registerExecutioner(TransientMultiAMR);  
  
    //Post Processors

    // Actions

  }
}

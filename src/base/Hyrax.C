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
#include "SolidMechanicsModule.h"
#include "PhaseFieldModule.h"
#include "HeatConductionModule.h"

//Kernels
#include "ACBulkCoupled.h"
#include "CHBulkCoupled.h"
#include "ACTransformElasticDF.h"
#include "CHBulkPolyCoupled.h"
#include "ACBulkPolyCoupled.h"
#include "ACNucleus.h"
#include "ACNucleusCNG.h"

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
#include "MaterialCNG.h"

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

    // Register Elk Modules
    Elk::PhaseField::registerObjects();
    Elk::SolidMechanics::registerObjects();
    Elk::HeatConduction::registerObjects();

    // Associate Syntax from SolidMechanics Module
    Elk::SolidMechanics::associateSyntax();

    //Kernels
    registerKernel(CHBulkCoupled);
    registerKernel(ACBulkCoupled);
    registerKernel(ACTransformElasticDF);
    registerKernel(ACBulkPolyCoupled);
    registerKernel(CHBulkPolyCoupled);
    registerKernel(ACNucleus);
    registerKernel(ACNucleusCNG);

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
    registerMaterial(MaterialCNG);

    //Initial Conditions
    registerInitialCondition(PolySpecifiedSmoothCircleIC);

    //Dampers

    //Executioners
    registerExecutioner(TransientMultiAMR);

    //Post Processors

    // Actions

  }
}

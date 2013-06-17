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
#include "Hyrax.h"
#include "HyraxApp.h"
#include "Factory.h"
#include "AppFactory.h"
#include "ActionFactory.h"

//Kernels
#include "ACBulkCoupled.h"
#include "CHBulkCoupled.h"
#include "ACTransformElasticDF.h"
#include "CHBulkPolyCoupled.h"
#include "ACBulkPolyCoupled.h"
#include "Value.h"
#include "CHBulkSimmons.h"
#include "ACBulkSimmons.h"
#include "CHCoupledCalphad.h"
#include "ACCoupledCalphad.h"

//Auxiliary Kernels
#include "AuxNucleationProbability.h"
#include "AuxNucleationRate.h"
#include "AuxSupersaturation.h"
#include "AuxChemElastic.h"
#include "AuxDeltaGStar.h"
#include "AuxRateSimple.h"

//#include "AuxTestFlip.h"
#include "ReporterAux.h"

//Dirac Kernels

//Boundary Conditions
#include "StressBC.h"

//Materials
#include "PFMobilityLandau.h"
#include "LinearSingleCrystalPrecipitateMaterial.h"
#include "ZrHCalphad.h"

//Initial Conditions
#include "PolySpecifiedSmoothCircleIC.h"

//Dampers

//Executioners
#include "MeshSolutionModify.h"

//Post Processors
#include "NucleationPostprocessor.h"
#include "OneSeed.h"

//Actions

//UserObjects
#include "NucleationLocationUserObject.h"
#include "NucleusIntroductionSolutionModifier.h"
#include "OneNucleusUserObject.h"

//Markers
#include "NucleationMarker.h"
#include "ErrorFractionMaxHMarker.h"


namespace Hyrax
{
  void registerApps()
  {
    registerApp(HyraxApp);
  }

  void registerObjects(Factory & factory)
  {
    //Kernels
    registerKernel(CHBulkCoupled);
    registerKernel(ACBulkCoupled);
    registerKernel(ACTransformElasticDF);
    registerKernel(ACBulkPolyCoupled);
    registerKernel(CHBulkPolyCoupled);
    registerKernel(Value);
    registerKernel(ACBulkSimmons);
    registerKernel(CHBulkSimmons);
    registerKernel(CHCoupledCalphad);
    registerKernel(ACCoupledCalphad);

    //Auxiliary Kernels
    registerAux(AuxSupersaturation);
    registerAux(AuxNucleationRate);
    registerAux(AuxNucleationProbability);
    registerAux(AuxChemElastic);
    registerAux(AuxDeltaGStar);
    registerAux(ReporterAux);
    registerAux(AuxRateSimple);

    //Dirac Kernels

    //Boundary Conditions
    registerBoundaryCondition(StressBC);

    //Materials
    registerMaterial(PFMobilityLandau);
    registerMaterial(LinearSingleCrystalPrecipitateMaterial);
    registerMaterial(ZrHCalphad);

    //Initial Conditions
    registerInitialCondition(PolySpecifiedSmoothCircleIC);

    //Dampers

    //Executioners
    registerExecutioner(MeshSolutionModify);

    //Postprocessors
    registerPostprocessor(NucleationPostprocessor);
    registerPostprocessor(OneSeed);

    // Actions

    // UserObjects
    registerUserObject(NucleationLocationUserObject);
    registerUserObject(NucleusIntroductionSolutionModifier);
    registerUserObject(OneNucleusUserObject);

    // Markers
    registerMarker(NucleationMarker);
    registerMarker(ErrorFractionMaxHMarker);
  }
}

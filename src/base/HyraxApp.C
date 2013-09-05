/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  30 May 2012
*
*************************************************************************/


#include "Moose.h"
#include "HyraxApp.h"
#include "Factory.h"
#include "AppFactory.h"
#include "ActionFactory.h"

//Elk Includes
#include "SolidMechanicsModule.h"
#include "TensorMechanicsModule.h"
#include "PhaseFieldModule.h"
#include "HeatConductionModule.h"
#include "MiscModule.h"


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
//#include "ACFerroelectric.h"

//Auxiliary Kernels
#include "AuxNucleationProbability.h"
#include "AuxNucleationRate.h"
#include "AuxSupersaturation.h"
#include "AuxChemElastic.h"
#include "AuxDeltaGStar.h"
#include "AuxRateSimple.h"
#include "AuxChem.h"

//#include "AuxTestFlip.h"
#include "ReporterAux.h"

//Dirac Kernels

//Boundary Conditions
#include "StressBC.h"

//Materials
#include "PFMobilityLandau.h"
#include "LinearSingleCrystalPrecipitateMaterial.h"
#include "ZrHCalphad.h"
//#include "FerroelectricBulk.h"

//Initial Conditions
#include "PolySpecifiedSmoothCircleIC.h"

//Dampers

//Executioners
#include "MeshSolutionModify.h"

//Post Processors
#include "NucleationPostprocessor.h"
#include "OneSeed.h"
#include "NucleiInformation.h"

//Actions
#include "OPVariantKernelAction.h"
#include "NucleationAuxAction.h"
#include "NucleationPostprocessorAction.h"

//UserObjects
#include "NucleationLocationUserObject.h"
#include "NucleusIntroductionSolutionModifier.h"
#include "OneNucleusUserObject.h"

//Markers
#include "NucleationMarker.h"
#include "ErrorFractionMaxHMarker.h"


template<>
InputParameters validParams<HyraxApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

HyraxApp::HyraxApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  Moose::registerObjects(_factory);
  HyraxApp::registerObjects(_factory);

  // Register Elk Modules
  Elk::PhaseField::registerObjects(_factory);
  Elk::SolidMechanics::registerObjects(_factory);
  Elk::TensorMechanics::registerObjects(_factory);
  Elk::HeatConduction::registerObjects(_factory);

  // Associate Syntax from SolidMechanics Module
  Moose::associateSyntax(_syntax, _action_factory);
  Elk::SolidMechanics::associateSyntax(_syntax, _action_factory);
  Elk::TensorMechanics::associateSyntax(_syntax, _action_factory);
  Elk::HeatConduction::associateSyntax(_syntax, _action_factory);
  Elk::Misc::associateSyntax(_syntax, _action_factory);

  //Associate syntax for Hyrax Actions
  HyraxApp::associateSyntax(_syntax, _action_factory);
}


void
HyraxApp::registerApps()
{
  registerApp(HyraxApp);
}

void
HyraxApp::registerObjects(Factory & factory)
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
  //registerKernel(ACFerroelectric);

  //Auxiliary Kernels
  registerAux(AuxSupersaturation);
  registerAux(AuxNucleationRate);
  registerAux(AuxNucleationProbability);
  registerAux(AuxChemElastic);
  registerAux(AuxDeltaGStar);
  registerAux(ReporterAux);
  registerAux(AuxRateSimple);
  registerAux(AuxChem);

  //Dirac Kernels

  //Boundary Conditions
  registerBoundaryCondition(StressBC);

  //Materials
  registerMaterial(PFMobilityLandau);
  registerMaterial(LinearSingleCrystalPrecipitateMaterial);
  registerMaterial(ZrHCalphad);
  //registerMaterial(FerroelectricBulk);

  //Initial Conditions
  registerInitialCondition(PolySpecifiedSmoothCircleIC);

  //Dampers

  //Executioners
  registerExecutioner(MeshSolutionModify);

  //Postprocessors
  registerPostprocessor(NucleationPostprocessor);
  registerPostprocessor(OneSeed);
  registerPostprocessor(NucleiInformation);

  // UserObjects
  registerUserObject(NucleationLocationUserObject);
  registerUserObject(NucleusIntroductionSolutionModifier);
  registerUserObject(OneNucleusUserObject);

  // Markers
  registerMarker(NucleationMarker);
  registerMarker(ErrorFractionMaxHMarker);
}

void
HyraxApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  // Actions
  registerAction(OPVariantKernelAction, "add_kernel");
  registerAction(NucleationAuxAction, "add_aux_kernel");
  registerAction(NucleationPostprocessorAction, "add_postprocessor");

  syntax.registerActionSyntax("OPVariantKernelAction", "OPVariantKernel");
  syntax.registerActionSyntax("NucleationAuxAction", "NucleationAux");
  syntax.registerActionSyntax("NucleationPostprocessorAction", "NucleationPostprocessor");
}

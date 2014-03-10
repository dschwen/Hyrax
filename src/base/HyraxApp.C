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

//Module Includes
#include "SolidMechanicsApp.h"
#include "TensorMechanicsApp.h"
#include "PhaseFieldApp.h"
#include "HeatConductionApp.h"
#include "MiscApp.h"

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
#include "CHLarry.h"
#include "CHLarrySplit.h"
#include "Heat.h"
#include "CHCoupledSplit.h"
#include "CHCoupledCalphadSplit.h"
#include "SplitCoupledCHWRes.h"
#include "CHPrecipMatrixElasticity.h"
#include "ACPrecipMatrixElasticity.h"

//Auxiliary Kernels
#include "AuxNucleationProbability.h"
#include "AuxNucleationRate.h"
#include "AuxSupersaturation.h"
#include "AuxChemElastic.h"
#include "AuxDeltaGStar.h"
#include "AuxRateSimple.h"
#include "AuxChem.h"
#include "AuxTemperature.h"
#include "AuxGuoEnergy.h"
#include "AuxCalphadEnergy.h"
//#include "AuxTestFlip.h"
#include "ReporterAux.h"
#include "AuxBulkEnergyCalphad.h"
#include "AuxGradientEnergy.h"
#include "AuxElasticEnergy.h"
#include "AuxFullNucleationRate.h"

//Dirac Kernels

//Boundary Conditions
#include "StressBC.h"

//Materials
#include "PFMobilityLandau.h"
#include "LinearSingleCrystalPrecipitateMaterial.h"
//#include "FerroelectricBulk.h"
#include "ZrHCalphad.h"
#include "CalphadEnergyMaterial.h"
#include "CalphadAB1CD1Material.h"
#include "CalphadAB1CD2Material.h"
#include "PrecipitateMatrixMisfitMaterial.h"
#include "ZrHCalphadDiffusivity.h"

//Initial Conditions
#include "PolySpecifiedSmoothCircleIC.h"
#include "EllipsoidIC.h"

//Dampers

//Executioners
#include "MeshSolutionModify.h"

//Post Processors
#include "NucleationPostprocessor.h"
#include "OneSeed.h"
#include "NucleiInformation.h"

//TimeSteppers
#include "InitialSolutionAdaptiveDT.h"
#include "PhaseFractionDT.h"

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

  // Register Modules
  PhaseFieldApp::registerObjects(_factory);
  SolidMechanicsApp::registerObjects(_factory);
  TensorMechanicsApp::registerObjects(_factory);
  HeatConductionApp::registerObjects(_factory);
  MiscApp::registerObjects(_factory);

  // Associate Syntax from SolidMechanics Module
  Moose::associateSyntax(_syntax, _action_factory);
  PhaseFieldApp::associateSyntax(_syntax, _action_factory);
  SolidMechanicsApp::associateSyntax(_syntax, _action_factory);
  TensorMechanicsApp::associateSyntax(_syntax, _action_factory);
  HeatConductionApp::associateSyntax(_syntax, _action_factory);
  MiscApp::associateSyntax(_syntax, _action_factory);

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
  registerKernel(CHLarry);
  registerKernel(CHLarrySplit);
  registerKernel(Heat);
  registerKernel(CHCoupledSplit);
  registerKernel(CHCoupledCalphadSplit);
  registerKernel(SplitCoupledCHWRes);
  registerKernel(CHPrecipMatrixElasticity);
  registerKernel(ACPrecipMatrixElasticity);

  //Auxiliary Kernels
  registerAux(AuxSupersaturation);
  registerAux(AuxNucleationRate);
  registerAux(AuxNucleationProbability);
  registerAux(AuxChemElastic);
  registerAux(AuxDeltaGStar);
  registerAux(ReporterAux);
  registerAux(AuxRateSimple);
  registerAux(AuxChem);
  registerAux(AuxTemperature);
  registerAux(AuxGuoEnergy);
  registerAux(AuxCalphadEnergy);
  registerAux(AuxBulkEnergyCalphad);
  registerAux(AuxGradientEnergy);
  registerAux(AuxElasticEnergy);
  registerAux(AuxFullNucleationRate);

  //Dirac Kernels

  //Boundary Conditions
  registerBoundaryCondition(StressBC);

  //Materials
  registerMaterial(PFMobilityLandau);
  registerMaterial(LinearSingleCrystalPrecipitateMaterial);
  registerMaterial(ZrHCalphad);
  //registerMaterial(FerroelectricBulk);
  registerMaterial(CalphadEnergyMaterial);
  registerMaterial(CalphadAB1CD1Material);
  registerMaterial(CalphadAB1CD2Material);
  registerMaterial(PrecipitateMatrixMisfitMaterial);
  registerMaterial(ZrHCalphadDiffusivity);

  //Initial Conditions
  registerInitialCondition(PolySpecifiedSmoothCircleIC);
  registerInitialCondition(EllipsoidIC);

  //Dampers

  //Executioners
  registerExecutioner(MeshSolutionModify);

  //Postprocessors
  registerPostprocessor(NucleationPostprocessor);
  registerPostprocessor(OneSeed);
  registerPostprocessor(NucleiInformation);

  //TimeSteppers
  registerTimeStepper(InitialSolutionAdaptiveDT);
  registerTimeStepper(PhaseFractionDT);

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

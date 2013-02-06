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
#include "Factory.h"
#include "ActionFactory.h"

//Kernels
#include "ACBulkCoupled.h"
#include "CHBulkCoupled.h"
#include "ACTransformElasticDF.h"
#include "CHBulkPolyCoupled.h"
#include "ACBulkPolyCoupled.h"
//#include "ACNucleus.h"
//#include "ACNucleusCNG.h"
//#include "ACInterfaceNucleation.h"
//#include "ACBulkNucleation.h"
#include "Value.h"
//#include "ValueNucleation.h"
//#include "ForcingFunctionNucleation.h"
//#include "TimeDerivativeNucleation.h"
#include "CHBulkSimmons.h"
#include "ACBulkSimmons.h"

//Auxiliary Kernels
//#include "AuxNucleation.h"
#include "AuxNucleationProbability.h"
#include "AuxNucleationRate.h"
#include "AuxSupersaturation.h"
//#include "AuxTestFlip.h"
#include "ReporterAux.h"

//Dirac Kernels
//#include "DiracNucleation.h"
//#include "DiracForcedAMR.h"

//Boundary Conditions
#include "StressBC.h"

//Materials
#include "PFMobilityLandau.h"
#include "LinearSingleCrystalPrecipitateMaterial.h"
//#include "MaterialCNG.h"

//Initial Conditions
#include "PolySpecifiedSmoothCircleIC.h"

//Dampers

//Executioners
//#include "TransientMultiAMR.h"
#include "MeshSolutionModify.h"

//Post Processors
#include "NucleationPostprocessor.h"
//#include "ValuePlusOne.h"
//#include "MaxElementValue.h"
#include "OneSeed.h"

//Actions

//UserObjects
#include "NucleationLocationUserObject.h"
#include "NucleusIntroductionSolutionModifier.h"

//Markers
#include "NucleationMarker.h"
#include "ErrorFractionMaxHMarker.h"


namespace Hyrax
{
  void registerObjects(Factory & factory)
  {
    //Kernels
    registerKernel(CHBulkCoupled);
    registerKernel(ACBulkCoupled);
    registerKernel(ACTransformElasticDF);
    registerKernel(ACBulkPolyCoupled);
    registerKernel(CHBulkPolyCoupled);
    //registerKernel(ACNucleus);
    //registerKernel(ACNucleusCNG);
    //registerKernel(ACInterfaceNucleation);
    //registerKernel(ACBulkNucleation);
    registerKernel(Value);
    //registerKernel(ValueNucleation);
    //registerKernel(ForcingFunctionNucleation);
    //registerKernel(TimeDerivativeNucleation);
    registerKernel(ACBulkSimmons);
    registerKernel(CHBulkSimmons);

    //Auxiliary Kernels
    //registerAux(AuxNucleation);
    registerAux(AuxSupersaturation);
    registerAux(AuxNucleationRate);
    registerAux(AuxNucleationProbability);
    //registerAux(AuxTestFlip);
    registerAux(ReporterAux);

    //Dirac Kernels
    // registerDiracKernel(DiracNucleation);
    //registerDiracKernel(DiracForcedAMR);

    //Boundary Conditions
    registerBoundaryCondition(StressBC);

    //Materials
    registerMaterial(PFMobilityLandau);
    registerMaterial(LinearSingleCrystalPrecipitateMaterial);
    //registerMaterial(MaterialCNG);

    //Initial Conditions
    registerInitialCondition(PolySpecifiedSmoothCircleIC);

    //Dampers

    //Executioners
    //registerExecutioner(TransientMultiAMR);
    registerExecutioner(MeshSolutionModify);

    //Postprocessors
    registerPostprocessor(NucleationPostprocessor);
    //registerPostprocessor(ValuePlusOne);
    //registerPostprocessor(MaxElementValue);
    registerPostprocessor(OneSeed);


    // Actions

    // UserObjects
    registerUserObject(NucleationLocationUserObject);
    registerUserObject(NucleusIntroductionSolutionModifier);

    // Markers
    registerMarker(NucleationMarker);
    registerMarker(ErrorFractionMaxHMarker);
  }
}

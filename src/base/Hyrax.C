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
#include "ACNucleus.h"
#include "ACNucleusCNG.h"
#include "ACInterfaceNucleation.h"
#include "ACBulkNucleation.h"
#include "Value.h"
#include "ValueNucleation.h"
#include "ForcingFunctionNucleation.h"
#include "TimeDerivativeNucleation.h"

//Auxiliary Kernels
#include "AuxNucleation.h"
#include "AuxNucleationProbability.h"
#include "AuxNucleationRate.h"
#include "AuxSupersaturation.h"
#include "AuxTestFlip.h"
#include "ReporterAux.h"

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
#include "ChangeVariableData.h"
#include "NucleationPostprocessor.h"

//Actions


namespace Hyrax
{
  void registerObjects()
  {
    //Kernels
    registerKernel(CHBulkCoupled);
    registerKernel(ACBulkCoupled);
    registerKernel(ACTransformElasticDF);
    registerKernel(ACBulkPolyCoupled);
    registerKernel(CHBulkPolyCoupled);
    registerKernel(ACNucleus);
    registerKernel(ACNucleusCNG);
    registerKernel(ACInterfaceNucleation);
    registerKernel(ACBulkNucleation);
    registerKernel(Value);
    registerKernel(ValueNucleation);
    registerKernel(ForcingFunctionNucleation);
    registerKernel(TimeDerivativeNucleation);

    //Auxiliary Kernels
    registerAux(AuxNucleation);
    registerAux(AuxSupersaturation);
    registerAux(AuxNucleationRate);
    registerAux(AuxNucleationProbability);
    registerAux(AuxTestFlip);
    registerAux(ReporterAux);

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
    registerPostprocessor(ChangeVariableData);
    registerPostprocessor(NucleationPostprocessor);

    // Actions
  }
}

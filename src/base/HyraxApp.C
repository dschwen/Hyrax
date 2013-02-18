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
#include "Hyrax.h"

//Elk Includes
#include "SolidMechanicsModule.h"
#include "TensorMechanicsModule.h"
#include "PhaseFieldModule.h"
#include "HeatConductionModule.h"
#include "MiscModule.h"


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
  Hyrax::registerObjects(_factory);

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
}

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


HyraxApp::HyraxApp(int argc, char * argv []) :
    MooseApp(argc, argv)
{
  init();
  Hyrax::registerObjects();
  // Register Elk Modules
  Elk::PhaseField::registerObjects();
  Elk::SolidMechanics::registerObjects();
  Elk::TensorMechanics::registerObjects();
  Elk::HeatConduction::registerObjects();

  // Associate Syntax from SolidMechanics Module
  Elk::SolidMechanics::associateSyntax(_syntax);
  Elk::TensorMechanics::associateSyntax(_syntax);
  Elk::HeatConduction::associateSyntax(_syntax);
  Elk::Misc::associateSyntax(_syntax);
}

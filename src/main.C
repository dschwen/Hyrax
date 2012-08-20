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

#include "HyraxApp.h"
#include "MooseInit.h"
#include "Moose.h"

// libMesh includes
#include "perf_log.h"

// Create a performance log
PerfLog Moose::perf_log("HYRAX");

// Begin the main program.
int main (int argc, char** argv)
{
  // Create a MooseInit Object
  MooseInit init (argc, argv);
  HyraxApp app(argc, argv);

  app.run();

  return 0;
}

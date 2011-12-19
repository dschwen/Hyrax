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

#include "Hyrax.h"

//Moose includes
#include "MooseInit.h"
#include "Executioner.h"

// Parser includes
#include "Parser.h"
#include "MooseSyntax.h"

// libMesh includes
#include "perf_log.h"

// Create a performance log
PerfLog Moose::perf_log("HYRAX");

// Begin the main program.
int main (int argc, char** argv)
{
  // Create a MooseInit Object
  MooseInit init (argc, argv);

  Hyrax::registerObjects();

  // Associate Parser Syntax with specific MOOSE Actions
  Moose::associateSyntax();

  // Create a parser object
  Parser p(Moose::syntax);

  // Parse commandline and return inputfile filename if appropriate
  std::string input_filename = p.parseCommandLine();

  // Tell the parser to parse the given file to setup the simulation and execute
  p.parse(input_filename);
  p.execute();

  Moose::executioner->execute();

  return 0;
}

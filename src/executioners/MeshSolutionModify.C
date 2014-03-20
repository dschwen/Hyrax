/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  11 December 2012
*
*  This code is originally by Derek Gaston, in moose_test - just been
*  copied and minorly modified.
*
*************************************************************************/

#include "MeshSolutionModify.h"
#include "NucleationLocationUserObject.h"
#include "TimeStepper.h"

// C++
#include <iostream>
#include <sstream>
#include <iomanip>

template<>
InputParameters validParams<MeshSolutionModify>()
{
  InputParameters params = validParams<Transient>();
  params.addParam<unsigned int>("adapt_cycles", 1, "# of adaptivity cycles to do normally.");
  params.addParam<unsigned int>("adapt_nucleus", 1, "# of adaptivity cycles to do with nucleus introduction.");
  params.addParam<bool>("use_nucleation_userobject", false, "Whether to pull in the nucleation user object or not");
  params.addParam<UserObjectName>("nucleation_userobject","", "The name of the UserObject to use for nucleation event locations");

  return params;
}

MeshSolutionModify::MeshSolutionModify(const std::string & name, InputParameters parameters) :
    Transient(name, parameters),
    _adapt_cycles(getParam<unsigned int>("adapt_cycles")),
    _adapt_nucleus(getParam<unsigned int>("adapt_nucleus")),
    _use_nucleation_userobject(getParam<bool>("use_nucleation_userobject")),
    _new_nucleus(false)
{
}

void
MeshSolutionModify::endStep()
{
  //std::cout<<"in MeshSolutionModify::endStep()"<<std::endl;
  _last_solve_converged = lastSolveConverged();

  // bool new_nucleus = false;
  unsigned int num_cycles;

 if(_use_nucleation_userobject)
  _new_nucleus = _nucleation_userobject->hasNewNucleus();

  if(_new_nucleus)
    num_cycles = _adapt_nucleus;
  else
    num_cycles = _adapt_cycles;

  if (_last_solve_converged)
  {
    for(unsigned int i=0; i<num_cycles; i++)
    {
      // Compute the Error Indicators and Markers
      _problem.computeIndicatorsAndMarkers();

#ifdef LIBMESH_ENABLE_AMR
      if (_problem.adaptivity().isOn())
      {
        std::cout<<"_problem.adaptMesh()"<<std::endl;
        _problem.adaptMesh();
        _problem.out().meshChanged();
      }
#endif
    }
    _problem.computeUserObjects(EXEC_CUSTOM);

    // perform the output for the current timestep
    _output_warehouse.outputStep();

    // if _at_sync_point is true, force the output no matter what
    _problem.output(_at_sync_point);
    _problem.outputPostprocessors(_at_sync_point);
    _problem.outputRestart(_at_sync_point);
  }


  std::cout<<"end of MeshSolutionModify::endStep()\n\n"<<std::endl;
}

void
MeshSolutionModify::preExecute()
{
  Transient::preExecute();

  // extra add-on here for nasty hacking of the user object system
  if (_use_nucleation_userobject)
    _nucleation_userobject = &getUserObject<NucleationLocationUserObject>("nucleation_userobject");
}

Real
MeshSolutionModify::getDT()
{

  Real new_dt = _dt;

  //If there's a nucleation event, cut the timestep back down to the timestepper input timestep
  //so the quickly varying solution is captured properly
  if (!_new_nucleus) //no event
  {
    new_dt = Transient::getDT();
  }
  else
  {
    if( lastSolveConverged() )
      new_dt = _time_stepper->getParam<Real>("dt");
    else
      new_dt = Transient::getDT();
  }

  return new_dt;
}

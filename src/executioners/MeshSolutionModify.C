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
#include <ostream>
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

/*void
MeshSolutionModify::takeStep(Real input_dt)
{
  _dt_old = _dt;
  if (input_dt == -1.0)
    _dt = computeConstrainedDT();
  else
    _dt = input_dt;

  _problem.onTimestepBegin();
  if (_converged)
  {
    // Update backward material data structures
    _problem.updateMaterials();
  }

  // Increment time
  _time = _time_old + _dt;

  std::cout<<"DT: "<<_dt<<std::endl;

  std::cout << " Solving time step ";
  {
    std::ostringstream out;

    out << std::setw(2)
    << _t_step
    << ", time="
    << std::setw(9)
    << std::setprecision(6)
    << std::setfill('0')
    << std::showpoint
    << std::left
    << _time
    <<  "...";

    std::cout << out.str() << std::endl;
  }

  preSolve();

  _problem.timestepSetup();

  // Compute Pre-Aux User Objects (Timestep begin)
  _problem.computeUserObjects(EXEC_TIMESTEP_BEGIN, UserObjectWarehouse::PRE_AUX);

  // Compute TimestepBegin AuxKernels
  _problem.computeAuxiliaryKernels(EXEC_TIMESTEP_BEGIN);

  // Compute Post-Aux User Objects (Timestep begin)
  _problem.computeUserObjects(EXEC_TIMESTEP_BEGIN, UserObjectWarehouse::POST_AUX);

  _problem.solve();

  _converged = _problem.converged();

  // We know whether or not the nonlinear solver thinks it converged, but we need to see if the executioner concurs
  bool last_solve_converged = lastSolveConverged();

  std::cout << "Converged:" << last_solve_converged << "\n";

  if (last_solve_converged)
    _problem.computeUserObjects(EXEC_TIMESTEP, UserObjectWarehouse::PRE_AUX);

  // User definable callback
  postSolve();

  _problem.onTimestepEnd();

  if (last_solve_converged)
  {
    _problem.computeAuxiliaryKernels(EXEC_TIMESTEP);
    _problem.computeUserObjects(EXEC_TIMESTEP, UserObjectWarehouse::POST_AUX);
  }
  } */

void
MeshSolutionModify::endStep()
{
  std::cout<<"in MeshSolutionModify::endStep()"<<std::endl;

  // bool new_nucleus = false;
  unsigned int num_cycles;

 if(_use_nucleation_userobject)
  _new_nucleus = _nucleation_userobject->hasNewNucleus();

  if(_new_nucleus)
    num_cycles = _adapt_nucleus;
  else
    num_cycles = _adapt_cycles;

  if (lastSolveConverged())
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

    // if _reset_dt is true, force the output no matter what
    _problem.output(_reset_dt);
    _problem.outputPostprocessors(_reset_dt);

    _time_old = _time;
    _t_step++;

    _problem.copyOldSolutions();
  }
  else
  {
    _time_stepper->rejectStep();
    _problem.getNonlinearSystem()._time_scheme->rejectStep();
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
MeshSolutionModify::computeDT()
{

  Real new_dt = _dt;

  // if(_use_nucleation_userobject)
//  _new_nucleus = _nucleation_userobject->hasNewNucleus();

  //If there's a nucleation event, cut the timestep back down to the input timestep so the quickly varying
  //solution is captured properly
  if (!_new_nucleus) //no event
  {
    new_dt = Transient::computeDT();
  }
  else
    new_dt = getParam<Real>("dt");

  // I also had the idea of chopping the timestep way down when a nucleation event occurs and only doing
  // one adaptivity step per timestep, and doing a few tiny timesteps to get in the adaptivity, then
  // going back to a larger (and growing) timestep.  This is kind of messy though, so I'm going to go with
  // this method first.

  return new_dt;
}

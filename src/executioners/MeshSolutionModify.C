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

#include <ostream>

template<>
InputParameters validParams<MeshSolutionModify>()
{
  InputParameters params = validParams<Transient>();
  params.addParam<unsigned int>("adapt_cycles", 1, "# of adaptivity cycles to do normally.");
  params.addParam<unsigned int>("adapt_nucleus", 1, "# of adaptivity cycles to do with nucleus introduction.");
  params.addRequiredParam<UserObjectName>("nucleation_userobject", "The name of the UserObject to use for nucleation event locations");

  return params;
}

MeshSolutionModify::MeshSolutionModify(const std::string & name, InputParameters parameters) :
    Transient(name, parameters),
    _adapt_cycles(getParam<unsigned int>("adapt_cycles")),
    _adapt_nucleus(getParam<unsigned int>("adapt_nucleus"))//,
//   _nucleation_userobject(getUserObject<NucleationLocationUserObject>("nucleation_userobject"))
{
}

void
MeshSolutionModify::endStep()
{
  //bool new_nucleus = _nucleation_userobject->hasNewNucleus();
  unsigned int num_cycles;

  //if(new_nucleus)
  //  num_cycles = _adapt_nucleus;
  //else
    num_cycles = _adapt_cycles;

  if (lastSolveConverged())
  {
    for(unsigned int i=0; i<num_cycles; i++)
    {
      //std::cout<<"adapting mesh"<<std::endl;
      // Compute the Error Indicators and Markers
      _problem.computeIndicatorsAndMarkers();

#ifdef LIBMESH_ENABLE_AMR
      if (_problem.adaptivity().isOn())
      {
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
    _problem.restoreSolutions();
}

//void
//MeshSolutionModify::preExecute()
//{
//  Transient::preExecute();

  // extra add-on here for nasty hacking of the user object system
//  _nucleation_userobject = &getUserObject<NucleationLocationUserObject>("nucleation_userobject");
//}

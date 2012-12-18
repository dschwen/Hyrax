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
#include <ostream>

//Moose includes

template<>
InputParameters validParams<MeshSolutionModify>()
{
  InputParameters params = validParams<Transient>();
  params.addParam<unsigned int>("adapt_cycles", 1, "# of adaptivity cycles to do.");
//  params.addRequiredParam<unsigned int>("max_h_level", "maximum refinement level");
  return params;
}

MeshSolutionModify::MeshSolutionModify(const std::string & name, InputParameters parameters) :
    Transient(name, parameters),
    _adapt_cycles(getParam<unsigned int>("adapt_cycles"))//,
    // _max_h_level(getParam<unsigned int>("max_h_level"))
{
}

void
MeshSolutionModify::endStep()
{
  if (lastSolveConverged())
  {
    for(unsigned int i=0; i<_adapt_cycles; i++)
    {
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

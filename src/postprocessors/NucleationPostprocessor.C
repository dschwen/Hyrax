/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  10 July 2012
*
*************************************************************************/

#include "NucleationPostprocessor.h"
#include "MooseVariable.h"
#include "SubProblem.h"
#include "MooseMesh.h"
#include "NonlinearSystem.h"

#include <ostream>

/**
 * This postprocessor is designed to take the nodal aux variable for
 * nucleation probability, calculate whether nucleation occurs or not at
 * each node and put the positive hits into a list that also includes
 * start and "end" time, then use those positives as the center of
 * nucleus.
 */

template<>
InputParameters validParams<NucleationPostprocessor>()
{
  InputParameters params = validParams <ChangeVariableData>();

  //something here

  // This postprocessor will couple to the AuxVariable for nucleation
  // probability.  No other aux variabled need to be coupled.

  // This postprocessor will also couple to the Kernel Variable order
  // parameter, which is manipulated in the postprocessor.  No other
  // variables are needed.

  params.addRequiredParam<Real>("radius", "nucleus radius");
  params.addRequiredParam<Real>("dwell_time", "How long nucleation event is");
  params.addRequiredParam<Real>("seed_value", "The nucleus seed order parameter value,");

  return params;
}

NucleationPostprocessor::NucleationPostprocessor(const std::string & name, InputParameters parameters)
    : ChangeVariableData(name, parameters),
      _radius(getParam<Real>("radius")),
      _dwell_time(getParam<Real>("radius")),
      _seed_value(getParam<Real>("seed_value")),
      _counter(0)
{
}

void
NucleationPostprocessor::initialize()
{
  Moose::seed(8675309 + _counter);
  _counter++;
}

void
NucleationPostprocessor::execute()
{
  searchForNucleationEvents();

  changeValues();
}

Real
NucleationPostprocessor::getValue()
{
  return Real(_nucleation_locations.size());
}

void
NucleationPostprocessor::searchForNucleationEvents()
{
  MeshBase::const_node_iterator it_end = _mesh.local_nodes_end();
  MeshBase::const_node_iterator it = _mesh.local_nodes_begin();

  // node loop to pick up nucleation locations
  for ( ; it != it_end ; ++it)
  {
    Node *node = *it;

    _subproblem.reinitNode(node, 0);

    // Pull out the value of the nucleation probability at this node
    Real probability = _coupled.getNodalValue(*node);

    // Generate a random number based on the node id and timestep
    Real random_number = Moose::rand();

    // make sure we're not trying to nucleate 2nd phase in a pre-existing 2nd phase
    if(probability > 0 && random_number < probability)
    {
      /* The trick here is to resize the vector that holds the nucleation events
       * locations each time there is a new event and tack it on to the end. */

      int s(_nucleation_locations.size());

      // resize the nucleation locations vector
      _nucleation_locations.resize(s+1);

      // fill in with the point location of the current node
      _nucleation_locations[s] = node;

      // resize the time vectors and type vector
      _start_times.resize(s+1);
      _end_times.resize(s+1);
      _orientation_type.resize(s+1);

      // fill in the time vectors with the start and end times for the new point
      _start_times[s] = _t;
      _end_times[s] = _t + _dwell_time;

      // test to see which orientation type the new phase is; assuming equal
      // probability for each type
      bool set(false);
      Real r_num = Moose::rand();
      for(int j=0; j< _moose_variable.size(); j++)
      {
        Real bin = (Real(j)+1.0)/Real(_moose_variable.size());

        if(set != true && r_num <= bin)
        {
          std::cout<<"r_num="<<r_num<<std::endl;
          std::cout<<"bin="<<bin<<std::endl;

          _orientation_type[s] = j;
          set = true;

          std::cout<<"orientationtype="<<_orientation_type[s]<<std::endl;
        }
      }
    }
  }
}

void
NucleationPostprocessor::changeValues()
{
  MeshBase::const_node_iterator it_end = _mesh.local_nodes_end();
  MeshBase::const_node_iterator it = _mesh.local_nodes_begin();

  // node loop to introduce nuclei into order parameter field
  for ( ; it != it_end ; ++it)
  {
    Node *node = *it;

    _subproblem.reinitNode(node, 0);

    // check the node against each nucleation point and see if it lives
    // within distance of a nucleus
    Real distance;
    for(unsigned int j(0); j<_nucleation_locations.size(); j++)
    {
      distance = (*_nucleation_locations[j] - *node).size();
      if(distance <=_radius &&
         _t >= _start_times[j] &&
         _t < _end_times[j])
      {
        int orientation = _orientation_type[j];

        _moose_variable[_orientation_type[j]]->setNodalValue(_seed_value);

        // Not sure if we need this, but probably :)
        _moose_variable[orientation]->insert(_nl.solution());
      }
    }
  }
  _nl.solution().close();
  _nl.sys().update();
}

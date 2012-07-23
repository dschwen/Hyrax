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
#include "GeneratedMesh.h"
#include "FEProblem.h"

#include "periodic_boundaries.h"

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

  // This postprocessor will couple to the AuxVariable for nucleation
  // probability.  No other aux variabled need to be coupled.

  // This postprocessor will also couple to one ore more Kernel Variable order
  // parameters, which is manipulated in the postprocessor.  No other
  // variables are needed.

  params.addRequiredParam<Real>("radius", "nucleus radius");
  params.addRequiredParam<Real>("dwell_time", "How long nucleation event is");
  params.addRequiredParam<Real>("seed_value", "The nucleus seed order parameter value,");

  return params;
}

NucleationPostprocessor::NucleationPostprocessor(const std::string & name, InputParameters parameters) :
    ChangeVariableData(name, parameters),
      _radius(getParam<Real>("radius")),
      _dwell_time(getParam<Real>("radius")),
      _seed_value(getParam<Real>("seed_value")),
      _counter(0),
      _phase_gen_index(std::numeric_limits<unsigned int>::max()),
      _gen_mesh(dynamic_cast<GeneratedMesh *>(&_mesh))
{
}

void
NucleationPostprocessor::initialize()
{
  // Initialize the periodic variable for each of our dimensions to false
  _periodic_dim.resize(LIBMESH_DIM, false);

  if (_gen_mesh)
  {
    for (unsigned int i=0; i<_gen_mesh->dimension(); ++i)
      // Assumption: We are going to assume that all variables are periodic together
      _periodic_dim[i] = _gen_mesh->isPeriodic(_nl, _moose_variable[0]->number(), i);

    _half_range = Point(_gen_mesh->dimensionWidth(0)/2.0, _gen_mesh->dimensionWidth(1)/2.0, _gen_mesh->dimensionWidth(2)/2.0);
  }

  for (unsigned int i=0; i<LIBMESH_DIM; ++i)
    std::cout << "Periodic in " << i << ": " << _periodic_dim[i] << "\n";

  _counter++;
  _local_start_times.clear();
  _local_end_times.clear();
  _local_orientation_type.clear();
  _local_node_ids.clear();
}

void
NucleationPostprocessor::modifySolutionVector()
{
  searchForNucleationEvents();

  // Gather all the shared data onto each processor
  Parallel::allgather(_local_start_times, false);
  Parallel::allgather(_local_end_times, false);
  Parallel::allgather(_local_orientation_type, false);
  Parallel::allgather(_local_node_ids, false);

  std::copy(_local_start_times.begin(), _local_start_times.end(), std::back_inserter(_start_times));
  std::copy(_local_end_times.begin(), _local_end_times.end(), std::back_inserter(_end_times));
  std::copy(_local_orientation_type.begin(), _local_orientation_type.end(), std::back_inserter(_orientation_type));

  // retrieve the nodes out of the mesh for use in the changeValues function
  _nucleation_locations.reserve(_nucleation_locations.size() + _local_node_ids.size());
  for (unsigned int i=0; i<_local_node_ids.size(); ++i)
    _nucleation_locations.push_back(&_mesh.node(_local_node_ids[i]));

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
    unsigned int node_id = node->id();

    // Our seed is complicated, it needs to be different for each node, each timestep so we have to take strides of n_nodes
    // This could potentially overflow but we'll deal with that later
    _mrand.seed(node_id, node_id + (_counter * _mesh.n_nodes()));
    Real random_number = _mrand.rand(node_id);

    // make sure we're not trying to nucleate 2nd phase in a pre-existing 2nd phase
    if(probability > 0 && random_number < probability)
    {
      /* The trick here is to resize the vector that holds the nucleation events
       * locations each time there is a new event and tack it on to the end. */

      // resize the nucleation locations vector and fill with point location of current node
      _local_node_ids.push_back(node_id);

      // resize the time vectors and type vectors and fill with current information
      _local_start_times.push_back(_t);
      _local_end_times.push_back(_t + _dwell_time);

      // test to see which orientation type the new phase is; assuming equal
      // probability for each type.  To generate a consistent number across
      // processors we will simply seed another generator with the node id.
      _mrand.seed(_phase_gen_index, node_id);

      bool set(false);
      int r_num = _mrand.randl(_phase_gen_index);
      r_num = r_num%_moose_variable.size();

      for(int j=0; j< _moose_variable.size(); j++)
      {
        int bin =j;

        if(set != true && r_num <= bin)
        {
          std::cout<<"r_num="<<r_num<<std::endl;
          std::cout<<"bin="<<bin<<std::endl;

          _local_orientation_type.push_back(bin);
          set = true;

          int s = _local_orientation_type.size();
          std::cout<<"orientationtype="<<_local_orientation_type[s-1]<<std::endl;
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
      distance = minPeriodicDistance(*_nucleation_locations[j], *node); //  (*_nucleation_locations[j] - *node).size();
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
}

Real
NucleationPostprocessor::minPeriodicDistance(Point p, Point q)
{
  for (unsigned int i=0; i<LIBMESH_DIM; ++i)
  {
    // check to see if we're closer in real or periodic space in x, y, and z
    if (_periodic_dim[i])
    {
      if (std::pow(p(i) - q(i), 2) > std::pow(_half_range(i), 2))
      {
        p(i) = p(i) - _half_range(i);
        q(i) = q(i) + _half_range(i);
      }
    }
  }

  return (p-q).size();
}

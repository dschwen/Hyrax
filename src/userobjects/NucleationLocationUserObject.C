/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  7 December 2012
*
*************************************************************************/

#include "NucleationLocationUserObject.h"
#include <ostream>

template<>
InputParameters validParams<NucleationLocationUserObject>()
{
  InputParameters params = validParams<GeneralUserObject>();

  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "The variable(s) we want to change");
  params.addRequiredParam<std::string>("coupled_aux", "The aux variable that we want to couple in");
  params.addRequiredParam<Real>("dwell_time", "How long nucleation event is");
  params.set<MooseEnum>("execute_on") = "timestep_begin";

  return params;
}

NucleationLocationUserObject::NucleationLocationUserObject(const std::string & name, InputParameters parameters) :
    GeneralUserObject(name, parameters),
    _mesh(_subproblem.mesh()),
    _nl(static_cast<FEProblem &>(_subproblem).getNonlinearSystem()),
    _coupled(_subproblem.getVariable(0, getParam<std::string>("coupled_aux"))),
    _dwell_time(getParam<Real>("dwell_time")),
    _counter(0),
    _phase_gen_index(std::numeric_limits<unsigned int>::max())
{
  std::vector<NonlinearVariableName> vars = getParam<std::vector<NonlinearVariableName> >("variables");
  _moose_variable.resize(vars.size());

  // initialize our vector of variable pointers
  for (unsigned int i=0; i<vars.size(); ++i)
    _moose_variable[i] = &_subproblem.getVariable(0, vars[i]);
}

void
NucleationLocationUserObject::initialize()
{
  // Assumption: We are going to assume that all variables are periodic together
  _mesh.initPeriodicDistanceForVariable(_nl, _moose_variable[0]->number());

  _counter++;
  _local_start_times.clear();
  _local_end_times.clear();
  _local_orientation_type.clear();
  _local_node_ids.clear();
}

void
NucleationLocationUserObject::execute()
{
  // we have to go over the local elements
  // Meshbase::const_element_iterator elem_it = _mesh.active_local_elements_begin();
  // Meshbase::const_element iterator elem_it_end = _mesh.active_local_elements_end();

  // element loop to pick up nucleation locations
  // for (; elem_it != elem_it_end; ++elem_it)
  //{
  //  const Elem *elem = *elem_it;
  //  // do we need this if and prepare statement?(from ElementalVariableValue)
  ///// if(elem->processor_id() == libMesh::processor_id())
  // {
  // _subproblem.prepare(elem, _tid);
  // //  _subproblem.reinitElem(elem, _tid); //probably need this
// pull out the value of the nucleation probability for this element


  //  //} //end of the possible If from ElementalVariableValue
  //}


  MeshBase::const_node_iterator it_end = _mesh.local_nodes_end();
  MeshBase::const_node_iterator it = _mesh.local_nodes_begin();

  //node loop to pick up nucleation locations
  for ( ; it != it_end ; ++it)
   {
     Node *node = *it;

     _subproblem.reinitNode(node, 0);

    // Pull out the value of the nucleation probability at this node
    Real probability = _coupled.getNodalValue(*node);

    // Generate a random number based on the node id and timestep
    unsigned int node_id = node->id();

    // Our seed is complicated, it needs to be different for each node, each timestep
    // so we have to take strides of n_nodes
    // This could potentially overflow but we'll deal with that later
    _mrand.seed(node_id, node_id + (_counter * _mesh.n_nodes()));
    Real random_number = _mrand.rand(node_id);

    // make sure we're not trying to nucleate 2nd phase in a pre-existing 2nd phase
    if(probability > 0 && random_number < probability)
    {
      // The trick here is to resize the vector that holds the nucleation events locations
      // each time there is a new event and tack it on to the end.

      // resize the nucleation locations vector and fill with point location of current node
      _local_node_ids.push_back(node_id);

      // resize the time vectors and type vectors and fill with current information
      _local_start_times.push_back(_t);
      _local_end_times.push_back(_t + _dwell_time);

      // test to see which orientation type the new phase is; assuming equal
      // probability for each type.  To generate a consistent number across
      // processors we will simply seed another generator with the node id.
      _mrand.seed(_phase_gen_index, node_id);

      int r_num = _mrand.randl(_phase_gen_index);

      // randl supplies some integer random number, we want to be between 1 and n coupled
      // vars, so modulo size()
      r_num = r_num % _moose_variable.size();

      _local_orientation_type.push_back(r_num);
    }
  }
}

void
NucleationLocationUserObject::finalize()
{
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
}

bool
NucleationLocationUserObject::elementWasHit(const Elem * elem) const
{
  bool was_hit = false;
  for(unsigned int i=0; i<_nucleation_locations.size(); i++)
  {
    was_hit = elem->contains_point(*_nucleation_locations[i]);
  if(was_hit)
  break;
  }
  return was_hit;
}

//const std::vector<Point> &
//NucleationLocationUserObject::hits()
//{
//return _nucleation_locations;
//}

const std::vector<Node* > &
NucleationLocationUserObject::getNucleationLocations() const
{
  return _nucleation_locations;
}

const std::vector<Real> &
NucleationLocationUserObject::getStartTimes() const
{
  return _start_times;
}

const std::vector<Real> &
NucleationLocationUserObject::getEndTimes() const
{
  return _end_times;
}

const std::vector<int> &
NucleationLocationUserObject::getOrientationType() const
{
  return _orientation_type;
}

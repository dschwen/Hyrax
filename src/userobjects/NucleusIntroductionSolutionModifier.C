/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  6 December 2012
*
*************************************************************************/

#include "NucleusIntroductionSolutionModifier.h"
#include "NucleationLocationUserObject.h"
#include <ostream>

template<>
InputParameters validParams<NucleusIntroductionSolutionModifier>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<UserObjectName>("nucleation_userobject", "The name of the UserObject to use for nucleation event locations");
  params.addRequiredParam<std::vector<VariableName> >("variables", "The name of the variable(s) to be modified");
  params.addRequiredParam<int>("num_vars", "The # of variables to modify");
  params.addRequiredParam<Real>("seed_value", "The nucleus order parameter value");
  params.addRequiredParam<Real>("radius", "nucleus radius");
  params.addParam<Real>("int_width", 0.0, "interface width of nucleus");

  return params;
}

NucleusIntroductionSolutionModifier::NucleusIntroductionSolutionModifier(const std::string & name,
                                                                         InputParameters parameters) :
    GeneralUserObject(name, parameters),
    _nucleation_userobject(getUserObject<NucleationLocationUserObject>("nucleation_userobject")),
    _mesh(_subproblem.mesh()),
    _nl(static_cast<FEProblem &>(_subproblem).getNonlinearSystem()),
    //_point_locator(_mesh.getMesh().point_locator()),
    //_variable(_subproblem.getVariable(0, parameters.get<VariableName>("modify"))),
    _radius(getParam<Real>("radius")),
    _seed_value(getParam<Real>("seed_value")),
    _int_width(getParam<Real>("int_width")),
    _num_vars(getParam<int>("num_vars"))
{
  std::vector<VariableName> vars = getParam<std::vector<VariableName> >("variables");
  _moose_variable.resize(vars.size());

  if(vars.size() != _num_vars)
    mooseError("Double check the # of variables you want to modify (NucleusIntroductionSolutionModifier).");

  // initialize our vector of variable pointers
  for (unsigned int i=0; i<vars.size(); ++i)
  {
    _moose_variable[i] = &_subproblem.getVariable(0, vars[i]);
  }

}

void
NucleusIntroductionSolutionModifier::execute()
{
  const std::vector<Real> & start_times = _nucleation_userobject.getStartTimes();
  const std::vector<Real> & end_times = _nucleation_userobject.getEndTimes();
  const std::vector<int> & orientation_type = _nucleation_userobject.getOrientationType();
  const std::vector<Node *> & nucleation_locations = _nucleation_userobject.getNucleationLocations();

  //The nucleation postprocessor does THIS:
  MeshBase::const_node_iterator it_end = _mesh.local_nodes_end();
  MeshBase::const_node_iterator it = _mesh.local_nodes_begin();

  // I believe this needs to be changed to work with constant monomial (elemental) aux variables.
  // node loop to introduce nuclei into order parameter field
  for ( ; it != it_end ; ++it)
  {
    Node *node = *it;

    _subproblem.reinitNode(node, 0);

    // check the node against each nucleation point and see if it lives
    // within distance of a nucleus
    Real distance;
    for(unsigned int j(0); j<nucleation_locations.size(); j++)
    {
      distance = _mesh.minPeriodicDistance(*nucleation_locations[j], *node);

      if( _t >= start_times[j] &&
          _t < end_times[j] )
      {
        int orientation = orientation_type[j];
        if(distance <=_radius - _int_width/2.0)
        {
          _moose_variable[orientation_type[j]]->setNodalValue(_seed_value);
        }
        else if(distance < _radius + _int_width/2.0)
        {
          Real interface_position = (distance - _radius + _int_width/2.0)/_int_width;
          Real int_value = _seed_value*(1.0 - std::cos(-interface_position + libMesh::pi/2.0));
          _moose_variable[orientation_type[j]]->setNodalValue(int_value);
        }

        // Not sure if we need this, but probably :)
        _moose_variable[orientation]->insert(_nl.solution());
      }
    }
  }
  _nl.solution().close();
  _nl.sys().update();
}


/*   //copy in the _random_hits internal vector from the RandomHits user object you coupled in
   const std::vector<Point> & locations = _nucleation_locations.hits();

  // *having nodes that were hit be the same size as the random hits - this would be the center of the nucleus...
  _nodes_that_were_hit.resize(locations.size());


  for(unsigned int i=0; i<locations.size(); i++)
  {
    // pull out the point information for each hit (nucleus center)
    const Point & location = locations[i];

    AutoPtr<PointLocatorBase> pl = PointLocatorBase::build(TREE, _mesh);

    // First find the element the hit lands in
    const Elem * elem = (*pl)(location);

    if(elem->processor_id() == libMesh::processor_id())
    {
      Real closest_distance = std::numeric_limits<unsigned int>::max();
      Node * closest_node = NULL;

//need to change this so it covers a finite area.
      // Find the node on that element that is closest.
      for(unsigned int n=0; n<elem->n_nodes(); n++)
      {
        Node * cur_node = elem->get_node(n);
        Real cur_distance = (hit - *cur_node).size();

        if(cur_distance < closest_distance)
        {
          closest_distance = cur_distance;
          closest_node = cur_node;
        }
      }

      _subproblem.reinitNode(closest_node, 0);
      _variable.setNodalValue(_variable.getNodalValue(*closest_node) + _seed_value);
      _variable.insert(_fe_problem.getNonlinearSystem().solution());
    }
  }

  _fe_problem.getNonlinearSystem().solution().close();
  _fe_problem.getNonlinearSystem().sys().update();


  }*/

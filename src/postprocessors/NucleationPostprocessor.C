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
      _seed_value(getParam<Real>("seed_value"))
{
}

void
NucleationPostprocessor::initialize()
{
  Moose::seed(8675309);
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

 MeshBase::const_element_iterator it_end = _mesh.active_local_elements_end();
  MeshBase::const_element_iterator it = _mesh.active_local_elements_begin();
//element/node loop to pick up nucleation locations
  for ( ; it != it_end ; ++it)
  {
    Elem *elem = *it;

    for (unsigned int i=0; i<elem->n_nodes(); ++i)
    {
      Node *node = elem->get_node(i);

      _subproblem.reinitNode(node, 0);

      // Pull out the value of the nucleation probability at this node
      Real probability = _coupled.getNodalValue(*node);

     // Generate a random number
      Real random_number(0.0);

      random_number = Moose::rand();

      std::cout <<"rand# = "<<random_number <<std::endl;
      std::cout<<"prob# = "<<probability <<std::endl<<std::endl;

      if(random_number < probability)
      {
        /* The trick here is to resize the vector that holds the
         * nucleation events locations each time there is a new event and
         * tack it on to the end. */
        std::cout<<"in if"<<std::endl;
        int s(_nucleation_locations.size());
        std::cout<<"s="<<s<<std::endl;
        // resize the nucleation locations vector
        _nucleation_locations.resize(s+1);
        std::cout<<"resized nuc loc, size ="<<_nucleation_locations.size()<<std::endl;
        // fill in with the point location of the current node
        _nucleation_locations[s] = *node;
        std::cout<<"filled nuc loc "<<_nucleation_locations[s]<<std::endl;
        // resize the time vectors
        _start_times.resize(s+1);
          std::cout<<"resized start"<<std::endl;
        _end_times.resize(s+1);

          std::cout<<"resized end"<<std::endl;
        // fill in the time vectors with the start and end times for the new point
        _start_times[s] = _t;
          std::cout<<"filled start"<<std::endl;
        _end_times[s] = _t + _dwell_time;
          std::cout<<"filled end"<<std::endl;
      }
    }
  }
}

void
NucleationPostprocessor::changeValues()
{
 MeshBase::const_element_iterator it_end = _mesh.active_local_elements_end();
  MeshBase::const_element_iterator it = _mesh.active_local_elements_begin();

//element/node loop to introduce nuclei into order parameter field
  //iterate over elements
  for ( ; it != it_end ; ++it)
  {
    Elem *elem = *it;

    // iterate over the nodes in the element
    for (unsigned int i=0; i<elem->n_nodes(); ++i)
    {
      Node *node = elem->get_node(i);
      _subproblem.reinitNode(node, 0);

      // check the node against each nucleation point and see if it lives
      // within distance of a nucleus
      Real distance;
      for(unsigned int j(0); j<_nucleation_locations.size(); j++)
      {
        distance = (_nucleation_locations[j] - *node).size();
        if(distance <=_radius &&
           _t >= _start_times[j] &&
           _t < _end_times[j])
        {
          _moose_variable.setNodalValue(_seed_value);

          // Not sure if we need this, but probably :)
          _moose_variable.insert(_nl.solution());
        }
      }
    }
  }
}

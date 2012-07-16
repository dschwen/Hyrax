/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  7 July 2012
*
*************************************************************************/

#include "ChangeVariableData.h"
#include "MooseVariable.h"
#include "SubProblem.h"
#include "MooseMesh.h"
#include "NonlinearSystem.h"

template<>
InputParameters validParams<ChangeVariableData>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<std::string>("variable", "The variable we want to query (we can always add more or make a vector, etc)");
  params.addRequiredParam<std::string>("coupled_aux", "The aux variable that we want to couple in");
  // We really want to only run this at the end of each timestep, so we'll force that here
  params.set<std::string>("execute_on") = "timestep";
  return params;
}

ChangeVariableData::ChangeVariableData(const std::string & name, InputParameters parameters) :
    GeneralPostprocessor(name, parameters),
    _mesh(_subproblem.mesh()),
    // Get this variable out of the subproblem, 0th thread (Postprocessors
    // not threaded)
    _moose_variable(_subproblem.getVariable(0, getParam<std::string>("variable"))),
    //be careful with this cast, it'll break if you change the line above
    _nl(*static_cast<NonlinearSystem *>(&_moose_variable.sys())),
    _variable_number(_moose_variable.number()),
    // Get the coupled variable reference
     _coupled(_subproblem.getVariable(0, getParam<std::string>("coupled_aux")))
{}

void
ChangeVariableData::initialize()
{
  _foo = 0.0;
}

void
ChangeVariableData::execute()
{
  /**
   * This loop is currently overriding the solution for our primary variable with the values from our
   * secondary variable used in the control structure.  We should be able to couple
   * in additional variables as needed and adjust our primary variable value as needed. :)
   */
  MeshBase::const_element_iterator it_end = _mesh.active_local_elements_end();
  MeshBase::const_element_iterator it = _mesh.active_local_elements_begin();
  for ( ; it != it_end ; ++it)
  {
    Elem *elem = *it;

    for (unsigned int i=0; i<elem->n_nodes(); ++i)
    {
      Node *node = elem->get_node(i);
      // DOF from the zeoreth system, we will retrieve the zeoreth component
      // unsigned int dof_number = node->dof_number(0, _variable_number, 0);

      // This call populates the datastructures in MooseVariable with the correct information
      // so that we can get/set Nodal information
      _subproblem.reinitNode(node, 0);

      // This is how you retrieve a value from the current solution
      Real coupled_value = _coupled.getNodalValue(*node);

//      Real value = _moose_variable.getNodalValue(*node);

      if(coupled_value > 1.0)
      {
        _moose_variable.setNodalValue(0.0);
        _foo += 1.0;
      }

//      _moose_variable.setNodalValue(1-value);

      // This is how you set a value in the current solution
      _moose_variable.setNodalValue(coupled_value+1);

      // Not sure if we need this, but probably :)
      _moose_variable.insert(_nl.solution());
    }
  }
  _nl.solution().close();
  _nl.sys().update();
}

Real
ChangeVariableData::getValue()
{
  return _foo;
}

void
ChangeVariableData::threadJoin(const Postprocessor & y)
{
}


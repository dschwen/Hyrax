#include "ChangeVariableData.h"
#include "MooseVariable.h"
#include "SubProblem.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<ChangeVariableData>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<std::string>("variable", "The variable we want to query (we can always add more or make a vector, etc)");

  return params;
}

ChangeVariableData::ChangeVariableData(const std::string & name, InputParameters parameters) :
    GeneralPostprocessor(name, parameters),
    _mesh(_subproblem.mesh()),
    _moose_variable(_subproblem.getVariable(0, getParam<std::string>("variable"))),   // Get this variable out of the zeroeth system (Nonlinear)
    _variable_number(_moose_variable.number())
{}

void
ChangeVariableData::initialize()
{
  _foo = 0;
}

void
ChangeVariableData::execute()
{
  {
    MeshBase::const_element_iterator it_end = _mesh.active_local_elements_end();
    MeshBase::const_element_iterator it = _mesh.active_local_elements_begin();
    for ( ; it != it_end ; ++it)
    {
      Elem *elem = *it;

      for (unsigned int i=0; i<elem->n_nodes(); ++i)
      {
        Node *node = elem->get_node(i);
        // DOF from the zeoreth system, we will retrieve the zeoreth component
        unsigned int dof_number = node->dof_number(0, _variable_number, 0);

        _foo += dof_number;
      }
    }
  }

  // Get the DOF out of the zeroeth system (Nonlinear), and the zeroeth component for scalars

  // TODO: Loop over elements or nodes and fill up a node pointer
  //
  // Node * node;
  // unsigned int dof_number = node->dof_number(0, _variable_number, 0);
  //
  // This is how you can view a value at a node
  // Real foo = _moose_variable.getNodalValue(*node);
  //
  // TODO: We need to actually get the DofObject so we can modify the value
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


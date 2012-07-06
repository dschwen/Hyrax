#include "ChangeVariableData.h"
#include "MooseVariable.h"
#include "SubProblem.h"

template<>
InputParameters validParams<ChangeVariableData>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<std::string>("variable", "The variable we want to query (we can always add more or make a vector, etc)");

  return params;
}

ChangeVariableData::ChangeVariableData(const std::string & name, InputParameters parameters) :
    GeneralPostprocessor(name, parameters),
    _moose_variable(_subproblem.getVariable(0, getParam<std::string>("variable"))),   // Get this variable out of the zeroeth system (Nonlinear)
    _variable_number(_moose_variable.number())
{}

void
ChangeVariableData::initialize()
{
}

void
ChangeVariableData::execute()
{
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
  return 0;
}

void
ChangeVariableData::threadJoin(const Postprocessor & y)
{
}


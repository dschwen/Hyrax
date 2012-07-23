#include "ValuePlusOne.h"
#include "MooseVariable.h"
#include "SubProblem.h"
#include "MooseMesh.h"
#include "NonlinearSystem.h"
#include "FEProblem.h"

template<>
InputParameters validParams<ValuePlusOne>()
{
  InputParameters params = validParams<ChangeVariableData>();
  return params;
}

ValuePlusOne::ValuePlusOne(const std::string & name, InputParameters parameters) :
    ChangeVariableData(name, parameters)
{
}


void
ValuePlusOne::initialize()
{
  _foo = 0.0;
}


void
ValuePlusOne::modifySolutionVector()
{
  MeshBase::const_element_iterator it_end = _mesh.active_local_elements_end();
  MeshBase::const_element_iterator it = _mesh.active_local_elements_begin();

  for ( ; it != it_end ; ++it)
  {

    Elem *elem = *it;

    for(unsigned int i=0; i<elem->n_nodes(); ++i)
    {
      Node *node = elem->get_node(i);

      _subproblem.reinitNode(node, 0);
      Real coupled_value = _coupled.getNodalValue(*node);

      if(coupled_value > 1.0)
      {
        _moose_variable[0]->setNodalValue(0.0);
        _foo += 1.0;
      }
      _moose_variable[0]->setNodalValue(coupled_value+1);
      _moose_variable[0]->insert(_nl.solution());
    }
  }
}


Real
ValuePlusOne::getValue()
{
  return _foo;
}

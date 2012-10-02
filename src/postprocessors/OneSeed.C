/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  26 September 2012
*
*************************************************************************/

#include "OneSeed.h"

#include "MooseVariable.h"
#include "SubProblem.h"
#include "MooseMesh.h"
#include "NonlinearSystem.h"
#include "GeneratedMesh.h"

#include <cmath>

#include <ostream>


template<>
InputParameters validParams<OneSeed>()
{
  InputParameters params = validParams <ChangeVariableData>();
  params.addRequiredParam<Real>("radius", "nucleus radius");
  params.addRequiredParam<Real>("dwell_time", "how long the nucleation event is");
  params.addRequiredParam<Real>("seed_value", "Nucleus seed order parameter value");
  params.addParam<Real>("int_width", 0.0, "nucleus interface width");
  params.addRequiredParam<Real>("x_position", "x-position of nucleus center");
  params.addRequiredParam<Real>("y_position", "y-position of nucleus center");
  params.addParam<Real>("z_position", 0.0, "z-position of the nucleus center");

  return params;
}

OneSeed::OneSeed(const std::string & name, InputParameters parameters) :
    ChangeVariableData(name, parameters),
    _radius(getParam<Real>("radius")),
    _dwell_time(getParam<Real>("dwell_time")),
    _seed_value(getParam<Real>("seed_value")),
    _int_width(getParam<Real>("int_width")),
    _x_position(getParam<Real>("x_position")),
    _y_position(getParam<Real>("y_position")),
    _z_position(getParam<Real>("z_position")),
    _location(_x_position, _y_position, _z_position),
    _counter(0),
    _gen_mesh(dynamic_cast<GeneratedMesh *>(&_mesh))
{
  std::cout<<"in constructor"<<std::endl;
}

void
OneSeed::initialize()
{
  _counter = 0;
}


void
OneSeed::modifySolutionVector()
{
  //Grab the mesh information
  MeshBase::const_node_iterator it_end = _mesh.local_nodes_end();
  MeshBase::const_node_iterator it = _mesh.local_nodes_begin();

  // node loop to introduce nucleus into order parameter field
  for ( ; it != it_end ; ++it)
  {
    Node *node = *it;

    _subproblem.reinitNode(node, 0);

    // check the node against the nucleation point and see if it lives within distance the nucleus
    Real distance;

    distance = _gen_mesh->minPeriodicDistance(_location, *node);
    if( _t < _dwell_time)
    {
      if(distance <=_radius - _int_width/2.0)
      {
        _moose_variable[0]->setNodalValue(_seed_value);
        _counter++;
      }
      else if(distance < _radius + _int_width/2.0)
      {
        Real interface_position = (distance - _radius + _int_width/2.0)/_int_width;
        Real int_value = _seed_value*(1.0 - std::cos(-interface_position + libMesh::pi/2.0));
        _moose_variable[0]->setNodalValue(int_value);
        _counter++;
      }
      _moose_variable[0]->insert(_nl.solution());
    }
  }
}

Real
OneSeed::getValue()
{
  return Real(_counter);
}

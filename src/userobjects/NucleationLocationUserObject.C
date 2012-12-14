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

#include "MooseVariable.h"
#include "SubProblem.h"
#include "MooseMesh.h"
#include "NonlinearSystem.h"
#include "GeneratedMesh.h"
#include "parallel.h"

#include <ostream>

template<>
InputParameters validParams<NucleationLocationUserObject>()
{
  InputParameters params = validParams<ElementUserObject>();

  params.addRequiredCoupledVar("coupled_aux", "coupled elemental auxiliary variable");
  params.addRequiredParam<Real>("dwell_time", "How long nucleation event is");
  params.set<MooseEnum>("execute_on") = "timestep_begin";
  params.addRequiredParam<int>("num_orientations", "# of orientation variants");

  return params;
}

NucleationLocationUserObject::NucleationLocationUserObject(const std::string & name, InputParameters parameters) :
    ElementUserObject(name, parameters),
    _mesh(_subproblem.mesh()),
    _coupled_probability(coupledValue("coupled_aux")),
    _dwell_time(getParam<Real>("dwell_time")),
    _num_orientations(getParam<int>("num_orientations")),
    _counter(0),
    _phase_gen_index(std::numeric_limits<unsigned int>::max())
{
}

void
NucleationLocationUserObject::initialize()
{
   _counter++;

   //reinitialize the local vectors to size 0 with no data
   _local_nucleus.clear();
}

void
NucleationLocationUserObject::execute()
{
  // here, we just look at an individual element

  //generate random # based on element id and timestep
  // random number needs to  be different for each element, each timestep
  // take steps of n_elem
  unsigned int elem_id = _current_elem->id();
  _mrand.seed(elem_id, elem_id + (_counter * _mesh.n_elem()));
  Real random_number = _mrand.rand(elem_id);

  //test for nucleation
  if (_coupled_probability[_qp] > 0 && random_number < _coupled_probability[_qp])
  {
    _local_nucleus.resize(1);

    // get the centroid of the element as the center of the nucleus
    Point nucleus_center = _current_elem->centroid();

    _local_nucleus[0].setLocation(nucleus_center);
    _local_nucleus[0].setStartTime(_t);
    _local_nucleus[0].setEndTime(_t+_dwell_time);

    _mrand.seed(_phase_gen_index, elem_id);

    int r_num = _mrand.randl(_phase_gen_index);

    // randl supplies some integer random number, we want to be between 1 and n coupled
    // vars, so modulo size()
    r_num = r_num%_num_orientations;

    _local_nucleus[0].setOrientation(r_num);
  }
}

void
NucleationLocationUserObject::finalize()
{
  // Gather all the shared data onto each processor
  // not sure what this call needs to be here
  // This is failing on compile
  //Parallel::allgather(_local_nucleus, false);

  //copy all the data into the final structure - not sure if this is right.
  std::copy(_local_nucleus.begin(), _local_nucleus.end(), std::back_inserter(_nuclei));
}

bool
NucleationLocationUserObject::elementWasHit(const Elem * elem) const
{
  bool was_hit = false;

  for(unsigned int i(0); i<_nuclei.size(); i++)
  {
    was_hit = elem->contains_point(_nuclei[i].getLocation());
    if(was_hit)
      break;
  }
  return was_hit;
}

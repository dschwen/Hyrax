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
  // params.addRequiredCoupledVar("coupled_variables", "coupled order parameter variables");
  //params.addRequiredParam<unsigned int>("n_OP_vars", "# of coupled OP variables");
  //params.addParam<Real>("nucleus_radius", 0.0, "the radius of the nucleus");

  return params;
}

NucleationLocationUserObject::NucleationLocationUserObject(const std::string & name, InputParameters parameters) :
    ElementUserObject(name, parameters),
    _mesh(_subproblem.mesh()),
    _coupled_probability(coupledValue("coupled_aux")),
    _dwell_time(getParam<Real>("dwell_time")),
    _num_orientations(getParam<int>("num_orientations")),
    _counter(0),
    _phase_gen_index(std::numeric_limits<unsigned int>::max()),
    _nuclei(0)//,
//    _n_OP_vars(getParam<unsigned int>("n_OP_vars")),
//    _nucleus_radius(getParam<Real>("nucleus_radius"))
{
/*    if(_n_OP_vars != coupledComponents("coupled_variables"))
    mooseError("Please match the number of orientation variants to coupled OPs (NucleationLocationUserObject).");

  _coupled_OP.resize(_n_OP_vars);

  for(unsigned int i=0; i<_n_OP_vars; i++)
  _coupled_OP[i] = &coupledValue("coupled_variables", i); */
}

void
NucleationLocationUserObject::initialize()
{
   _counter++;

   //reinitialize the local vectors to size 0 with no data
   _local_nucleus.clear();
   _packed_data.clear();
}

void
NucleationLocationUserObject::execute()
{
  // here, we just look at an individual element

  /**
   * generate random # based on element id and timestep. random number needs to
   * be different for each element, each timestep.  take steps of n_elem
   */
  unsigned int elem_id = _current_elem->id();
  _mrand.seed(elem_id, elem_id + (_counter * _mesh.n_elem()));
  Real random_number = _mrand.rand(elem_id);



  //test for nucleation
  if ((_coupled_probability[_qp] > 0) && (random_number < _coupled_probability[_qp]))
  {
    // get the centroid of the element as the center of the nucleus
    Point nucleus_center = _current_elem->centroid();

    Nucleus current_nucleus;
    current_nucleus.setLocation(nucleus_center);
    current_nucleus.setStartTime(_t);
    current_nucleus.setEndTime(_t+_dwell_time);

    _mrand.seed(_phase_gen_index, elem_id);
    int r_num = _mrand.randl(_phase_gen_index);

    /**
     * randl supplies some integer random number, we want to be between 1 and n coupled
     * vars, so modulo size()
     */
    r_num = r_num%_num_orientations;

    current_nucleus.setOrientation(r_num);
    _local_nucleus.push_back(current_nucleus);
  }
}

void
NucleationLocationUserObject::finalize()
{
  //pack up the data so it can be communicated using allgather
  //pack(_packed_data);
  Nucleus::pack(_local_nucleus, _packed_data);

  // Gather all the shared data onto each processor
  Parallel::allgather(_packed_data, false);

  //unpack all the data into the "global" variable
  Nucleus::unpack(_packed_data, _nuclei);

  //clean out instances of nuclei being too close to each other (overlapping)
  //cleanUp();
}

void
NucleationLocationUserObject::threadJoin(const UserObject &a)
{
  const NucleationLocationUserObject & nluo = dynamic_cast<const NucleationLocationUserObject &>(a);

  //Pack up the data on both threads
  // Nucleus::pack(_local_nucleus, _packed_data);

  //std::vector<Real> nluo_packed_data;
  std::vector<Nucleus> nluo_local_nucleus = nluo._local_nucleus;
  //nluo.pack(nluo_packed_data);

  //stick the two pieces of data together like peanut butter and jelly
  std::copy(nluo_local_nucleus.begin(), nluo_local_nucleus.end(), std::back_inserter(_local_nucleus));
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


/*void
NucleationLocationUserObject::cleanUp()
{
  // find the nuclei that overlap
  for(unsigned int i=0; i<_nuclei.size(); i++)
  {
  // must do for only this timestep
  // save their index, then in another loop, delete them out
  // however
  }

  }*/


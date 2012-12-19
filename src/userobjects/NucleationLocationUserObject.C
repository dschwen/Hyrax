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
    _phase_gen_index(std::numeric_limits<unsigned int>::max()),
    _nuclei(0)//,
    // _stride(6)
{
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
  if (_coupled_probability[_qp] > 0 && random_number < _coupled_probability[_qp])
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
NucleationLocationUserObject::pack(std::vector<Real> & packed_data) const
{
  // Don't repack data if it's already packed - could cause data loss, and that would make me sad.
  if (!packed_data.empty())
    return;

  //get the size of the _local_nucleus and multiply by 6, the # of individual pieces of data for each nucleus
  const unsigned int stride = Nucleus::stride();
  unsigned int packed_size = _local_nucleus.size()*stride;
  packed_data.resize(packed_size);

  unsigned int j=0;
  for(unsigned int i=0; i<packed_data.size(); i+=stride)
  {
    Point this_location = _local_nucleus[j].getLocation();
    int this_orientation = _local_nucleus[j].getOrientation();

    packed_data[i] = this_location(0); //x-coord of point
    packed_data[i+1] = this_location(1); //y-coord
    packed_data[i+2] = this_location(2); //z-coord
    packed_data[i+3] = _local_nucleus[j].getStartTime();
    packed_data[i+4] = _local_nucleus[j].getEndTime();
    packed_data[i+5] = Real(this_orientation); //need to make sure we won't get any type conversion errors...

    j++;
  }
  } */

 /*void
NucleationLocationUserObject::unpack(const std::vector<Real> & packed_data)
{
  const unsigned int stride = Nucleus::stride();
  unsigned int number_nuclei = packed_data.size()/stride;
  std::vector<Nucleus> new_nuclei(number_nuclei);

  Real tol = 1e-5;
  unsigned int j(0);
  for(unsigned int i=0; i<packed_data.size(); i+=stride)
  {
    Point current_point(packed_data[i], packed_data[i+1], packed_data[i+2]);

    new_nuclei[j].setLocation(current_point);
    new_nuclei[j].setStartTime(packed_data[i+3]);
    new_nuclei[j].setEndTime(packed_data[i+4]);

    Real local_orientation = packed_data[i+5];

    //make sure there's no type conversion errors for the orientation
    if(int(local_orientation + tol) == int(local_orientation))
      new_nuclei[j].setOrientation(int(local_orientation));
    else
      new_nuclei[j].setOrientation(int(local_orientation)+1);

    j++;
  }

  std::copy(new_nuclei.begin(), new_nuclei.end(), std::back_inserter(_nuclei));
  //might want to include some error message in here in case unpacking screws up
}
 */


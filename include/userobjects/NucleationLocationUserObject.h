/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  7 December 2012
*
*************************************************************************/

#ifndef NUCLEATIONLOCATIONUSEROBJECT_H
#define NUCLEATIONLOCATIONUSEROBJECT_H

#include "GeneralUserObject.h"
#include "MooseRandom.h"

//Forward declarations
class NucleationLocationUserObject;

template<>
InputParameters validParams<NucleationLocationUserObject>();

class NucleationLocationUserObject : public GeneralUserObject
{
public:
  NucleationLocationUserObject (const std::string & name, InputParameters parameters);
  virtual ~NucleationLocationUserObject() {}

  virtual void initialize();
  virtual void execute();
  virtual void destroy() {}
  virtual void finalize();

  bool elementWasHit(const Elem * elem) const;

  //const std::vector<Point> & hits() const;
  const std::vector<Node *> & getNucleationLocations() const;
  const std::vector<Real> & getStartTimes() const;
  const std::vector<Real> & getEndTimes() const;
  const std::vector<int> & getOrientationType() const;

protected:
private:

  MooseMesh & _mesh;
  NonlinearSystem & _nl;
  MooseVariable & _coupled;
  std::vector<MooseVariable *> _moose_variable;

  Real _dwell_time;
  int _counter;

  /**
   * Global data used for calculating when the solution needs to be modified
   */
  std::vector<Real> _start_times;
  std::vector<Real> _end_times;
  std::vector<int> _orientation_type;
  std::vector<Node *> _nucleation_locations;

  /**
   * The following vectors are filled by each local process and concatenated
   * onto the "global" arrays at the end of each step
   */
  std::vector<Real> _local_start_times;
  std::vector<Real> _local_end_times;
  std::vector<int> _local_orientation_type;
  // Node ids used to retieve nucleation location information
  std::vector<unsigned int> _local_node_ids;

  /// The Moose stateful random number generator
  MooseRandom _mrand;

  /// The index of the phase orientation generator (we will use a high index that isn't used by the node generators)
  const unsigned int _phase_gen_index;
};

#endif //NUCLEATIONLOCATIONUSEROBJECT_H

/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  10 July 2012
*
*************************************************************************/

#ifndef NUCLEATIONPOSTPROCESSOR_H
#define NUCLEATIONPOSTPROCESSOR_H

#include "ChangeVariableData.h"
#include "MooseStatefulRandom.h"

//forward declarations
class NucleationPostprocessor;

template<>
InputParameters validParams<NucleationPostprocessor>();

class NucleationPostprocessor : public ChangeVariableData
{
public:
  NucleationPostprocessor (const std::string & name, InputParameters parameters);

  virtual void initialize();

  virtual void execute();

  virtual Real getValue();

protected:
  void searchForNucleationEvents();

  void changeValues();

private:

  Real _radius;
  Real _dwell_time;
  Real _seed_value;

  int _counter;

  std::vector<Real> _start_times;
  std::vector<Real> _end_times;
  std::vector<int> _orientation_type;

  std::vector<Node *> _nucleation_locations;

  MooseStatefulRandom _mrand;
};

#endif //NUCLEATIONPOSTPROCESSOR_H

/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  30 May 2012
*
*************************************************************************/

#ifndef MATERIALCNG_H
#define MATERIALCNG_H

#include "Material.h"

class MaterialCNG;

template<>
InputParameters validParams<MaterialCNG>();

/// template specialization for MaterialProperty Init
template <>
PropertyValue *
MaterialProperty<std::vector<Point> >::init(int size);
///

class MaterialCNG : public Material
{
public:
  MaterialCNG(const std::string & name, InputParameters parameters);

protected:
  virtual void computeProperties();

  virtual void calculateSupersaturation(unsigned int qp);

  virtual void calculateNucleationRate(unsigned int qp);

  virtual void calculateNucleationProbability(unsigned int qp);

  virtual void testForNucleation(unsigned int qp);

  virtual void appendNucleationList(unsigned int qp);

private:
  //std::string _nucleation_field_name;

  // coupled nucleation field
  //VariableValue & _nucleation_field;

  MaterialProperty<Real> & _c1;
  VariableValue & _conc;
  VariableValue & _OP;


  MaterialProperty<Real> & _supersaturation;

  Real _Kn1;
  Real _Kn2;

  MaterialProperty<Real> & _j_star;
  MaterialProperty<Real> & _p_nm;
  MaterialProperty<bool> & _nucleation;

  // vector of locations where nucleation returned true
  MaterialProperty<std::vector<Point> > & _nucleation_locations;

  // vectors (paired with _nucleation_locations) of when the nucleation
  // event occurred and when it should stop
  MaterialProperty<std::vector<Real> > & _start_times;
  MaterialProperty<std::vector<Real> > & _end_times;

  // input of how long the nucleation event needs to take
  Real _dwell_time;

  // list to hold nucleation locations, gets put into MP
  std::vector<Point> _locations;

  // list to hold nucleation start and stop times, put into MP
  std::vector<Real> _start;
  std::vector<Real> _end;

  unsigned int _random_number_seed;

  // holds timestep information
  Real _local_time;

};

#endif //MATERIALCNG_H

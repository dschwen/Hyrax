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


class MaterialCNG : public Material
{
public:
  MaterialCNG(const std::string & name, InputParameters parameters);

protected:
  virtual void computeProperties();

private:
  //std::string _nucleation_field_name;

  // coupled nucleation field
  VariableValue & _nucleation_field;

  // vector of locations where nucleation returned true
  MaterialProperty<std::vector<RealVectorValue> > & _nucleation_locations;

  // vectors (paired with _nucleation_locations) of when the nucleation
  // event happened and when the nucleation event ends
  MaterialProperty<std::vector<Real> > & _start_times;
  MaterialProperty<std::vector<Real> > & _end_times;

  // input of how long the nucleation event needs to take
  Real _dwell_time;

  std::vector<RealVectorValue> _locations;
  std::vector<Real> _start;
  std::vector<Real> _end;

};

#endif //MATERIALCNG_H

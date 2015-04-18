/*************************************************************************
*************************************************************************/

#ifndef FREEENERGY_H
#define FREEENERGY_H

#include "Material.h"

//forward declaration
class FreeEnergy;

template<>
InputParameters validParams<FreeEnergy>();

class FreeEnergy : public Material
{
public:
  FreeEnergy(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();

  Real _length_scale;
  Real _energy_scale;

  MaterialProperty<Real> & _kappa_c;      //CH gradient energy coefficient (isotropic)
  MaterialProperty<Real> & _kappa_n;      //AC gradient energy coefficient (isotropic)

  MaterialProperty<Real> & _W;            //well height
  MaterialProperty<Real> & _molar_vol;    //molar volume

  MaterialProperty<Real> & _alpha_energy;
  MaterialProperty<Real> & _delta_energy;

  //COUPLED VARIABLES
  VariableValue & _c;   //coupled concentration
  VariableValue & _n;   //coupled order parameter
  VariableGradient & _grad_c;
  VariableGradient & _grad_n;

  MaterialProperty<Real> & _free_energy_density;  //free energy functional energy density
  MaterialProperty<Real> & _constant_phases_energy_density;  //free energy density of ONLY the alpha+delta phases

private:

};

#endif //FREEENERGY_H

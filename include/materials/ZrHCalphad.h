/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  13 June 2013
*
*************************************************************************/

#ifndef ZRHCALPHAD_H
#define ZRHCALPHAD_H

#include "Material.h"

//forward declaration
class ZrHCalphad;

template<>
InputParameters validParams<ZrHCalphad>();

class ZrHCalphad : public Material
{
public:
  ZrHCalphad(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();

//  Real computeMobility();
//  RealGradient computeGradMobility();

private:

  //INPUT VARIABLES
  //Diffusion coefficient information
  //Real _H_Zr_D0;
  //Real _H_ZrH2_D0;
  //Real _H_Zr_Q0;
  //Real _H_ZrH2_Q0;

  Real _mobility_CH;
  Real _mobility_AC;

  Real _kappa_CH;
  Real _kappa_AC;

  Real _well_height;

  Real _molar_volume;

  Real _thermal_diffusivity;
  Real _dThermal_diffusivity_dT;

  //MATERIAL PROPERTIES
  MaterialProperty<Real> & _M;			//Cahn-Hilliard mobility (isotropic)
  MaterialProperty<RealGradient> & _grad_M;

  MaterialProperty<Real> & _L;			//Allen-Cahn kinetic coefficient (isotropic)

  MaterialProperty<Real> & _kappa_c;		//CH gradient energy coefficient (isotropic)
  MaterialProperty<Real> & _kappa_n;		//AC gradient energy coefficient (isotropic)

  MaterialProperty<Real> & _W;                  //well height

  MaterialProperty<Real> & _molar_vol;          //molar volume

  MaterialProperty<Real> & _thermal_diff;       //thermal diffusivity of the material
  MaterialProperty<Real> & _dThermDiff_dT;      //d(Thermal diffusivity)/dTemperature

  //COUPLED VARIABLES
  VariableValue & _temperature;                 //coupled Temperature field

};

#endif //ZRHCALPHAD_H

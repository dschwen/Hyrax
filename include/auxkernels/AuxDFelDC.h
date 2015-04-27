/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
* 25 April 2015
*
*************************************************************************/

#ifndef AUXDFELDC_H
#define AUXDFELDC_H

#include "AuxKernel.h"
#include "ElasticityTensorR4.h"
#include "RankTwoTensor.h"

class AuxDFelDC;

template<>
InputParameters validParams<AuxDFelDC>();

class AuxDFelDC : public AuxKernel
{
public:
    AuxDFelDC(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

private:

  MaterialProperty<ElasticityTensorR4> & _elasticity_tensor;
  MaterialProperty<RankTwoTensor> & _elastic_strain;
  MaterialProperty<RankTwoTensor> & _dc_misfit_strain;

//Real _scaling_factor;

};

#endif //AUXDFELDC_H

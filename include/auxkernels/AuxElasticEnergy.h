/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
* 9 January 2014
*
*************************************************************************/

#ifndef AUXELASTICENERGY_H
#define AUXELASTICENERGY_H

#include "AuxKernel.h"
#include "ElasticityTensorR4.h"
#include "RankTwoTensor.h"

class AuxElasticEnergy;

template<>
InputParameters validParams<AuxElasticEnergy>();

class AuxElasticEnergy : public AuxKernel
{
public:
    AuxElasticEnergy(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

private:

  MaterialProperty<ElasticityTensorR4> & _elasticity_tensor;
  MaterialProperty<RankTwoTensor> & _elastic_strain;

//Real _scaling_factor;

};

#endif //AUXELASTICENERGY_H

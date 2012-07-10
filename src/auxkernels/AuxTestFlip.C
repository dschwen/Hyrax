/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  7 July 2012
*
*************************************************************************/

#include "AuxTestFlip.h"
#include <ostream>

template<>
InputParameters validParams<AuxTestFlip>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("coupled_c", "coupled variable: concentration");
  params.addRequiredParam<Real>("radius", "Radius of the value flip circle");
  params.addRequiredParam<RealVectorValue>("center_point", "center of the value flip circle");

  return params;
}

AuxTestFlip::AuxTestFlip(const std::string & name, InputParameters parameters)
    : AuxKernel(name, parameters),
      _coupled_c(coupledValue("coupled_c")),
      _radius(getParam<Real>("radius")),
      _center_point(getParam<RealVectorValue>("center_point"))
{
}

Real
AuxTestFlip::computeValue()
{

  //Point current_point;
  //current_point = _q_point[_qp];

  //Real distance;
  //distance = (current_point - _center_point).size();

  // std::cout<<"current qp="<<_qp<<std::endl;
  //std::cout<<"current point="<<current_point<<std::endl;
  //std::cout<<"center point="<<_center_point<<std::endl;
  //std::cout<<"distance="<<distance<<std::endl<<std::endl;

  //return 0.0 outside the area to flip
  if(_coupled_c[_qp] > 0.4)
    return 2.0;
  else
    return 0.0;
}

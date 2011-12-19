/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  15 December 2011
*
*  This code inherits from DiracKernel in MOOSE
*  
*  This code handles the introduction of nuclei as delta functions for 
*  the concurrent nucleation and growth algorithm first proposed by 
*  J.P. Simmons.
*  
*************************************************************************/
#include "DiracNucleation.h"

/* Remember, in the input file, tell this Dirac Kernel that it's operating on the appropriate
   variable - in this case, on the order parameter variable.  When I've got multiple OPs, this is
   going to get interesting. */

template<>
InputParameters validParams<DiracNucleation>()
{
  InputParameters params = validParams<DiracKernel>();
  params.addRequiredParam<Real>("value", "The value of the point source");
  params.addRequiredCoupledVar("nucleation", "Auxiliary variable: nucleation or not");
//  in input file, this should always be called "nucleation"
  return params;
}

DiracNucleation::DiracNucleation(const std::string & name, InputParameters parameters) :
    DiracKernel(name, parameters),
    _value(getParam<Real>("value")),
    _coupled_nucleation(coupledValue("nucleation"))
{
}

void
DiracNucleation::addPoints()
{
//// search through field variable to find where it is true: is the following correct?
//   for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//   {
//     /* in this case, auxvariable nucleation consists of 0.0s and 2.0s, for true and false - 
//     *   would use boolean but auxvariable will only return real.*/
//     if (_coupled_nucleation[_qp] > 1.0) 

//     {
//     // pull out the element information -> how/what? into _true_element
//     // pull out the point information -> how/what? into _true_point
//     addPoint(_true_element, _true_point);
//     }
//   }
}  

// /* this came from ConstantPointSource, to set the position of the point for x,y,z dimensions;
// *  it was in the constructor.  Should do something like this probably for setting my points 
// *  locations here, but up in addPoints().  */
//  _p(0) = _point_param[0];

//  if(_point_param.size() > 1)
//  {
//    _p(1) = _point_param[1];

//    if(_point_param.size() > 2)
//    {
//      _p(2) = _point_param[2];
//    }
//  }

Real
DiracNucleation::computeQpResidual()
{
// need to figure out what the math is here. The below line came from ConstantPointSource
//  return -_test[_i][_qp]*_value;

  return 0.0;
}



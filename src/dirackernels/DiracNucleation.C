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
  params.addRequiredCoupledVar("nucleation", "Auxiliary variable: nucleation or not");
//  in input file, this should always be called "nucleation"
  params.addRequiredParam<Real>("value", "The value of the point source");
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
  MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end();

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;
    _problem.prepare(elem, 0);
    _problem.reinitElem(elem, 0);

    for (unsigned int qp = 0; qp < _qrule->n_points(); qp++)
      if (_coupled_nucleation[qp] > 0.0)
      {
//        std::cout << "###DEBUG AddPoints: " << _coupled_nucleation[_qp] << " at " << _qrule->qp(qp) << "\n";
        addPoint(elem, _qrule->qp(qp));
      }
  }
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

//void
//DerivedDirac::addPoints()
//{




//}

Real
DiracNucleation::computeQpResidual()
{
//  This is negative because it's a forcing function that has been brought over to the left side
//  this needs to get changed.  actually, crap.  keep _value.
  return -_test[_i][_qp]*_value;
}

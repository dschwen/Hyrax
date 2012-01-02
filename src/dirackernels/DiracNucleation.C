/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  15 December 2011
*
*************************************************************************/
#include "DiracNucleation.h"

#include "elem.h"

/* Remember, in the input file, tell this Dirac Kernel that it's operating on the appropriate
   variable - in this case, on the order parameter variable.  When I've got multiple OPs, this is
   going to get interesting. */

/**
 *  DiracNucleation works with the AuxNucleation etc. system to determine where and when Dirac delta "spikes"
 *  are introduced into the order parameter field variable.  This is for the simulation of the explicit
 *  introduction of nuclei as delta functions for the concurrent nucleation and growth algorithm first 
 *  proposed by J.P. Simmons (2000). 
 */


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
    _coupled_nucleation(coupledValue("nucleation")),
    _my_refine(_mesh)
{
}

void
DiracNucleation::addPoints()
{
  MeshBase::const_element_iterator       el     = _mesh.active_local_elements_begin();
  MeshBase::const_element_iterator end_el = _mesh.active_local_elements_end();

  // iterate over the elements
  for ( ; el != end_el ; ++el)
  {
    Elem* elem = *el;
    _problem.prepare(elem, 0);
    _problem.reinitElem(elem, 0);

    //for (unsigned int qp = 0; qp < _qrule->n_points(); qp++)
    // shouldn't need the above line for average element value
      if (_coupled_nucleation[0] > 0.0)  
      /* access one quadrature point for the element - should give us the same value since it's an 
      element average value*/
      //CJP changed the if to >0.0 for rapid nucleation.  This should be >1.0
      {
//        std::cout << "###DEBUG AddPoints: " << _coupled_nucleation[_qp] << " at " << _qrule->qp(qp) << "\n";

        addPoint(elem, elem->centroid());
	// refine the element around the delta function
	elem->set_refinement_flag(libMesh::Elem::REFINE);
	_my_refine.refine_elements();
//	_problem.out().meshChanged();
        // this should add the point to the center of the element, which is fine for now
      }
  }
}


Real
DiracNucleation::computeQpResidual()
{
//  This is negative because it's a forcing function that has been brought over to the left side
//  this needs to get changed.  actually, crap.  keep _value.
  return -_test[_i][_qp]*_value;
}

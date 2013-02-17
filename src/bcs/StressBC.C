/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  3 January 2013
*
*  StressBC is pulled from Ferret, thanks to O. Heinonen et al.
*
****************************************************************/

#include "StressBC.h"

template<>
InputParameters validParams<StressBC>()
{
    InputParameters params = validParams<IntegratedBC>();
    params.addRequiredParam<int>("component","Which component(0 for x, 1 for y, 2 for z) in traction is used");

    params.addRequiredParam<std::vector<Real> >("boundary_stress", "Boundary stress: s11, s22, s33, s23, s13, s12");

    //params.addRequiredParam<Real>("stress_xx", "stress_xx");
    //params.addRequiredParam<Real>("stress_xy", "stress_xy");
    //params.addRequiredParam<Real>("stress_yy", "stress_yy");
    //params.addRequiredParam<Real>("stress_yz", "stress_yz");
    //params.addRequiredParam<Real>("stress_zx", "stress_zx");
    //params.addRequiredParam<Real>("stress_zz", "stress_zz");

    return params;
}

StressBC::StressBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    _stress_vector(getParam<std::vector<Real> >("boundary_stress")),
    _Jacobian_mult(getMaterialProperty<ElasticityTensorR4>("Jacobian_mult"/* + getParam<std::string>("appended_property_name")*/)),

    _component(getParam<int>("component"))
// _stress_xx(getParam<Real>("stress_xx")),
    //_stress_xy(getParam<Real>("stress_xy")),
    //_stress_yy(getParam<Real>("stress_yy")),
    //_stress_yz(getParam<Real>("stress_yz")),
    //_stress_zx(getParam<Real>("stress_zx")),
    //_stress_zz(getParam<Real>("stress_zz"))
{
  _boundary_stress.fillFromInputVector(_stress_vector);
}

Real
StressBC::computeQpResidual()
{
//  Real values[3][3];
  // Real traction[3];

  /*
  values[0][0]=_stress_xx;
  values[0][1]=_stress_xy;
  values[0][2]=_stress_zx;

  values[1][0]=_stress_xy;
  values[1][1]=_stress_yy;
  values[1][2]=_stress_yz;

  values[2][0]=_stress_zx;
  values[2][1]=_stress_yz;
  values[2][2]=_stress_zz;
  */

/*
  for(int i=0; i<3; ++i)
  {
    traction[i] = 0;
    for(int j=0; j<3; ++j)
      traction[i] = traction[i] + _boundary_stress[i][j]*_normals[_qp](j);
  }
*/

  return -_test[_i][_qp]*(_boundary_stress.row(_component)*_normals[_qp]);
  //-_test[_i][_qp]*traction[_component]; //be careful with the sign
}

Real
StressBC::computeQpJacobian()
{
//  return -_test[_i][_qp]*(_Jacobian_mult[_qp].stressBCJacobian(_component, _component, _grad_phi[_j][_qp]));

  return 0;
}

/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  13 June 2013
*
*************************************************************************/

#include "ZrHCalphad.h"

template<>
InputParameters validParams<ZrHCalphad>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredParam<std::vector<Real> >("G_hcp_Zr_vector", "A, B, C, D, E Gibbs coefficients: hcp Zr");
  params.addRequiredParam<std::vector<Real> >("G_fcc_Zr_vector", "A, B, C, D, E Gibbs coefficients: fcc Zr");
  params.addRequiredParam<std::vector<Real> >("G_hcp_ZrH_vector", "A, B Gibbs coefficients: hcp ZrH");
  params.addRequiredParam<std::vector<Real> >("G_fcc_ZrH2_vector", "A, B, C Gibbs coefficients: fcc ZrH2");
  params.addRequiredParam<std::vector<Real> >("G_H2_vector", "A, B, C, D, E Gibbs coefficients: H2 gas");

  params.addRequiredParam<std::vector<Real> >("L0_vector", "A, B Gibbs coefficients: L0 term in fcc ZrH2");
  params.addRequiredParam<std::vector<Real> >("L1_vector", "A, B Gibbs coefficients: L1 term in fcc ZrH2");

  params.addRequiredParam<Real>("H_Zr_D0", "Diffusion prefactor for H in hcp Zr (isotropic)");
  params.addRequiredParam<Real>("H_ZrH2_D0", "Diffusion prefactor for H in fcc ZrH2 (isotropic)");
  params.addRequiredParam<Real>("H_Zr_Q0", "Activation energy for H in hcp Zr");
  params.addRequiredParam<Real>("H_ZrH2_Q0", "Activaton energy for H in fcc ZrH2");

  params.addRequiredParam<Real>("mobility_AC", "mobility coefficient for Allen-Cahn equation (isotropic)");
  params.addRequiredParam<Real>("kappa_CH", "CH gradient energy coefficient (isotropic)");
  params.addRequiredParam<Real>("kappa_AC", "AC gradient energy coefficient (isotropic)");

  params.addRequiredCoupledVar("conserved_var", "Atomic fraction of hydrogen");
  params.addRequiredCoupledVar("OP_variable_names", "Array of order parameters,");
  params.addRequiredParam<int>("n_OP_variables", "# of coupled OP variables, >=1");

  //Default in J/mol-K
  params.addParam<Real>("gas_constant", 8.3144621, "Universal gas constant");
  params.addParam<Real>("temperature", 600, "Simulation temperature in K");

  return params;
}

ZrHCalphad::ZrHCalphad(const std::string & name, InputParameters parameters)
    : Material(name, parameters),
      _V_hcp_Zr(getParam<std::vector<Real> >("G_hcp_Zr_vector")),
      _V_fcc_Zr(getParam<std::vector<Real> >("G_fcc_Zr_vector")),
      _V_hcp_ZrH(getParam<std::vector<Real> >("G_hcp_ZrH_vector")),
      _V_fcc_ZrH2(getParam<std::vector<Real> >("G_fcc_ZrH2_vector")),
      _V_H2(getParam<std::vector<Real> >("G_H2_vector")),
      _V_L0(getParam<std::vector<Real> >("L0_vector")),
      _V_L1(getParam<std::vector<Real> >("L1_vector")),

      _H_Zr_D0(getParam<Real>("H_Zr_D0")),
      _H_ZrH2_D0(getParam<Real>("H_ZrH2_D0")),
      _H_Zr_Q0(getParam<Real>("H_Zr_Q0")),
      _H_ZrH2_Q0(getParam<Real>("H_ZrH2_Q0")),

      _mobility_AC(getParam<Real>("mobility_AC")),
      _kappa_CH(getParam<Real>("kappa_CH")),
      _kappa_AC(getParam<Real>("kappa_AC")),

      _T(getParam<Real>("temperature")),
      _R(getParam<Real>("gas_constant")),

      _n_OP_variables(getParam<int>("n_OP_variables")),

      _G_hcp_Zr(declareProperty<Real>("G_hcp_Zr")),
      _G_fcc_Zr(declareProperty<Real>("G_fcc_Zr")),
      _G_hcp_ZrH(declareProperty<Real>("G_hcp_ZrH")),
      _G_fcc_ZrH2(declareProperty<Real>("G_fcc_ZrH2")),
      _G_H2(declareProperty<Real>("G_H2")),
      _L0(declareProperty<Real>("L0")),
      _L1(declareProperty<Real>("L1")),

      _Zr_molar_volume(declareProperty<Real>("molar_volume_alpha_Zr")),
      _ZrH2_molar_volume(declareProperty<Real>("molar_volume_delta_ZrH2")),

      _M(declareProperty<Real>("M")),
      _grad_M(declareProperty<RealGradient>("grad_M")),
      _kappa_c(declareProperty<Real>("kappa_c")),
      _L(declareProperty<Real>("L")),
      _kappa_n(declareProperty<Real>("kappa_n")),

      _X(coupledValue("conserved_var")),
      _grad_X(coupledGradient("conserved_var"))
{
  //Make sure everything is set up from input correctly
  if((_V_hcp_Zr.size() != 5) ||
     (_V_fcc_Zr.size() != 5) ||
     (_V_hcp_ZrH.size() != 2) ||
     (_V_fcc_ZrH2.size() != 3) ||
     (_V_H2.size() != 5) ||
     (_V_L0.size() != 2) ||
     (_V_L1.size() != 2) )
    mooseError("Please supply the correct # of values for the Gibbs coefficients (ZrHCalphad).");

  if(_n_OP_variables != coupledComponents("OP_variable_names"))
    mooseError("Please match the # of orientation variants to coupled OPs (ZrHCalphad)");

  _OP.resize(_n_OP_variables);
  _grad_OP.resize(_n_OP_variables);

  for(unsigned int i=0; i< _n_OP_variables; i++)
  {
    _OP[i] = &coupledValue("OP_variable_names", i);
    _grad_OP[i] = &coupledGradient("OP_variable_names", i);
  }
  std::cout<<"In ZrHCalphad constructor"<<std::endl;
}

void
ZrHCalphad::computeQpProperties()
{
  // ALL UNITS IN SI except for length: that's in microns!

  //Rudlich-Kitsler terms
  _L0[_qp] = _V_L0[0] + _V_L0[1]*_T;
  _L1[_qp] = _V_L1[0] + _V_L1[1]*_T;

  _G_hcp_Zr[_qp] = computeGhcpZr();
  _G_fcc_Zr[_qp] = computeGfccZr();
  _G_hcp_ZrH[_qp] = computeGhcpZrH();
  _G_fcc_ZrH2[_qp] = computeGfccZrH2();
  _G_H2[_qp] = computeGH2();

  //Molar volume of pure Zr, in microns^3/mol
  _Zr_molar_volume[_qp] = (1.4e-5)*(std::pow(1.0e6, 3));

  //Molar volume of ZrH(2-x), in microns^3/mol
  _ZrH2_molar_volume[_qp] = (1.65e-5)*(std::pow(1.0e6, 3));

  //THIS NEEDS TO BE CHANGED
  //Real D;
  _D_alpha = _H_Zr_D0*std::exp(-_H_Zr_Q0/(_R*_T));
  _D_delta = _H_ZrH2_D0*std::exp(-_H_ZrH2_Q0/(_R*_T));

  //  std::cout<<"before computeMobilities()"<<std::endl;
  //if (_t > 0.0006)

  //  std::cout<<"_X[_qp] = "<<_X[_qp]<<std::endl;
  // std::cout<<"computeD2fchemDx2 = "<<computeD2fchemDx2()<<std::endl;
  // std::cout<<"absolute computeD2FchemDx2 = "<<std::abs(computeD2fchemDx2())<<std::endl;

// THIS WILL NEED FIXING - get ANALYTICAL form of the limit!
  // if(std::abs( computeD2fchemDx2() ) < 0.0000000001)
  // _M[_qp] = _H_coeff/0.0000000001;
  //else
  // _M[_qp] = _H_coeff/computeD2fchemDx2();
  //_grad_M[_qp] = computeGradMobility();
  //_M[_qp] = computeMobility();
  //_grad_M[_qp] = computeGradMobility();
  _M[_qp] = _D_delta;
  _grad_M[_qp] = 0.0;

//  _L[_qp] = _mobility_AC;
  _L[_qp] = _D_delta;
  _kappa_c[_qp] = _kappa_CH;
  _kappa_n[_qp] = _kappa_AC;

  //  std::cout<<"end of computeQpProperties"<<std::endl;
}

Real
ZrHCalphad::computeGhcpZr()
{
  //Calculates molar Gibbs free energy for Zr in HCP form - molar standard enthalpy for Zr
  return _V_hcp_Zr[0] + _V_hcp_Zr[1]*_T + _V_hcp_Zr[2]*_T*std::log(_T) + _V_hcp_Zr[3]*_T*_T
    + _V_hcp_Zr[4]/_T;
}


Real
ZrHCalphad::computeGfccZr()
{
  //Calculates molar Gibbs free energy for fcc Zr - molar standard enthalpy for Zr
  return  _V_fcc_Zr[0] + _V_fcc_Zr[1]*_T + _V_fcc_Zr[2]*_T*std::log(_T) + _V_fcc_Zr[3]*_T*_T
    + _V_fcc_Zr[4]/_T;
}


Real
ZrHCalphad::computeGhcpZrH()
{
  //Calculates molar Gibbs free energy for ZrH compound - molar standard enthalpies
  return _V_hcp_ZrH[0] + _V_hcp_ZrH[1]*_T + computeGhcpZr() + 0.5*computeGH2();
}


Real
ZrHCalphad::computeGfccZrH2()
{
  //Calculates molar Gibbs free energy for fcc ZrH2 - molar standard enthalpies
  return _V_fcc_ZrH2[0] + _V_fcc_ZrH2[1]*_T + _V_fcc_ZrH2[2]*_T*std::log(_T) + computeGhcpZr()
    + computeGH2();
}


Real
ZrHCalphad::computeGH2()
{
  //Calculates molar Gibbs free energy for H2 gas - molar standard enthalpy for H2
  return _V_H2[0] + _V_H2[1]*_T + _V_H2[2]*_T*std::log(_T) + _V_H2[3]*_T*_T + _V_H2[4]/_T;
}

Real
ZrHCalphad::computeMobility()
{
/*  Real diffusion, curvature;

  diffusion = _D_alpha*(1.0 - computeHeaviside()) + _D_delta*computeHeaviside();

  curvature = computeD2GalphaDx2()*(1.0 - computeHeaviside()) + computeD2GdeltaDx2()*computeHeaviside();

  // Hack to deal with asymtotic behavior
  if (std::abs(curvature) < 1e-15)
    curvature = 1e-15;

  if (std::abs(curvature) > 1e15)
    curvature = 1e15;

    return diffusion/curvature; */
  return 0;

}

RealGradient
ZrHCalphad::computeGradMobility()
{
  /*Real h;
  h = computeHeaviside();

//  std::cout<<"before computeCurvatures"<<std::endl;

  Real D2Alpha, D2Delta;
  D2Alpha = computeD2GalphaDx2();
  D2Delta = computeD2GdeltaDx2();

//  std::cout<<"before if"<<std::endl;

  Real curvature;
  curvature = D2Alpha*(1.0 - h) + D2Delta*h;
  if( curvature > 1e15)
    curvature = 1e15;
  else if( std::abs(curvature) < 1e-15)
    curvature = 1e-15;

//  std::cout<<"before conserved_term"<<std::endl;

  RealGradient conserved_term;
  conserved_term = -1*_grad_X[_qp]*(_D_alpha*(1 - h) + _D_delta*h)*( computeD3GalphaDx3()*(1.0 - h)
                                                      + computeD3GdeltaDx3()*h)/(curvature*curvature);
  RealGradient nonconserved_term;
  nonconserved_term.zero();

//  std::cout<<"before nonconserved_term"<<std::endl;
  Real first_term, second_term;
//  std::cout<<"before for"<<std::endl;
  for(unsigned int i(0); i < _n_OP_variables; ++i)
  {
//    std::cout<<"i = "<<i<<std::endl;

    first_term = (D2Alpha*(1 - h) + D2Delta*h)*( (_D_delta - _D_alpha)*(_d_heaviside[_qp])[i] );
    second_term = (_D_alpha*(1 - h) + _D_delta*h)*( (D2Delta - D2Alpha)*(_d_heaviside[_qp])[i]);

    nonconserved_term += (*_grad_OP[i])[_qp]*(first_term - second_term)/(curvature*curvature);
  }
//  std::cout<<"after for"<<std::endl;

  return conserved_term + nonconserved_term;
  */
  return 0;

}

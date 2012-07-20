/*************************************************************************
*
*  Welcome to HYRAX!
*  Andrea M. Jokisaari
*  CASL/MOOSE
*
*  30 May 2012
*
*************************************************************************/

#include "MaterialCNG.h"
#include "MooseRandom.h"

#include <ostream>

template<>
InputParameters validParams<MaterialCNG>()
{
  InputParameters params = validParams<Material>();
  //params.addRequiredCoupledVar("nucleation_field_name", "Same as the variable used for the AuxNucleation auxkernel");


  params.addRequiredCoupledVar("concentration_var", "coupled variable: concentration for supersaturation calculation");
  params.addRequiredCoupledVar("OP_var", "coupled variable: order parameter");
  params.addParam<Real>("Kn1", 0.0, "first nucleation rate coefficient");
  params.addParam<Real>("Kn2", 0.0, "second nucleation rate coefficient");
  params.addRequiredParam<Real>("dwell_time", "How long nucleation kernel should be applied to an event");

  return params;
}

MaterialCNG::MaterialCNG(const std::string & name, InputParameters parameters)
    : Material(name, parameters),
      _c1(getMaterialProperty<Real>("C1")),
      _conc(coupledValue("concentration_var")),
      _OP(coupledValue("OP_var")),
      _supersaturation(declareProperty<Real>("supersaturation")),
      _Kn1(getParam<Real>("Kn1")),
      _Kn2(getParam<Real>("Kn2")),
      _j_star(declareProperty<Real>("nucleation_rate")),
      _p_nm(declareProperty<Real>("nucleation_probability")),
      _nucleation(declareProperty<bool>("nucleation_true_false")),
      _nucleation_locations(declareProperty<std::vector<Point> >("nucleation_locations")),
      _start_times(declareProperty<std::vector<Real> >("start_times")),
      _end_times(declareProperty<std::vector<Real> >("end_times")),
      _dwell_time(getParam<Real>("dwell_time")),
      _local_time(-1.0)
{
  _random_number_seed = 1;
}


void
MaterialCNG::computeProperties()
{
  // operate only at the beginning of the timestep
  if(_local_time != _t)
  {
    for(unsigned int qp=0; qp<_qrule->n_points(); ++qp)
    {
      // operate only if the order parameter is low (not in 2nd phase)
      if(_OP[qp] < 0.2)
      {
        calculateSupersaturation(qp);
        calculateNucleationRate(qp);
        calculateNucleationProbability(qp);
        testForNucleation(qp);
      }
    }

    // fill in the list of points where nucleation occurred and the start and stop times
    for(unsigned int qp=0; qp<_qrule->n_points(); ++qp)
    {
      int s(_locations.size());

      // resize all the materials properties vectors
      _nucleation_locations[qp].resize(s);
      _start_times[qp].resize(s);
      _end_times[qp].resize(s);

      // fill in the vectors
      _nucleation_locations[qp] = _locations;
      _start_times[qp] = _start;
      _end_times[qp] = _end;
    }
    _local_time = _t;
  }
}


void
MaterialCNG::calculateSupersaturation(unsigned int qp)
{
  //supersaturation = concentration - terminal solid solubility of 2nd material in matrix
  _supersaturation[qp] = _conc[qp] - _c1[qp];
}


void
MaterialCNG::calculateNucleationRate(unsigned int qp)
{
  // nucleation rate = Kn1 * exp(-1*Kn2 / supersaturation)
  if(_supersaturation[qp] <= 0.0)
  {
    _j_star[qp] = 0.0;
  }
  else
  {
    _j_star[qp] = _Kn1*exp(-1.0*_Kn2/_supersaturation[qp]);
  }
}


void
MaterialCNG::calculateNucleationProbability(unsigned int qp)
{
   //  nucleation probability = 1 - exp(-1*nucleation rate*dt)
  _p_nm[qp] = 1.0 - exp(-1.0*_j_star[qp]*_dt);
}


void
MaterialCNG::testForNucleation(unsigned int qp)
{
  Real random_number(0.0);


  //reset the random number seed so that you keep generating the same
  // random number sequence each timestep.Although this will change the nucleation rate -
  // so make sure this random sequence is only tacking on info for the first iteration
  // of the timestep calculation.

  // compare large timestep to small timesteps - the number of nuclei.

  // there will be some small variation in the output data of j_star and p_nm because
  // they are calculated for each residual vs the nucleation/points info which will
  // only be calculated once.  this gonna be ok?

  _random_number_seed = _t;
  MooseRandom::seed(_random_number_seed);

  random_number = MooseRandom::rand();

  if(random_number < _p_nm[qp])
  {
    _nucleation[qp] = true;
    appendNucleationList(qp);
  }
  else
  {
    _nucleation[qp] = false;
  }
}


void
MaterialCNG::appendNucleationList(unsigned int qp)
{
  /* The trick here is to resize the vector that holds the nucleation events locations each time
   * there is a new event and tack it on to the end.  This completed list will then get put
   * into every QP for another material property. */

  int s(_locations.size());

  // resize the locations vector
  _locations.resize(++s);

  // fill in with the point location of the current qp
  _locations[s-1] = _q_point[qp];

  // resize the time vectors
  _start.resize(s);
  _end.resize(s);

  // fill in the time vectors with the start and end times for the new point
  _start[s-1] = _t;
  _end[s-1] = _t + _dwell_time;
}

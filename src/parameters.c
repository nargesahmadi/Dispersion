#include <stdio.h>
#include <math.h>
#include "dispersion.h"

void
init_parameters(void)
{
  
  // Number of species 
  // (up to 10 species are possible)
  params.nr_kinds = 2;

  //Maximum number of iterations in the Newton method
  params.numiter = 100;

  // Threshold for the determinant of the dispersion tensor:
  // If det D <= det_D_threshold, the Newton iteration will be stopped
  params.det_D_threshold = 1.d-20;

  // Maximum of sum in Bessel function
  // can be very low (e.g., 3) for quasi-parallel propagation
  params.nmax = 1000;

  //If I_n is less than this value, higher n are neglected:
  params.Bessel_zero = 1.e-50;

  // Parallel beta of the species
  double beta[] = { 1. , 1. , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

  // Temperature anisotropy (Tperp/Tparallel)
  double aniso[] = { 2. , 1. , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

  // Charge of the species in units of the first ion charge
  double charge[] = { 1. , -1. , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

  // Mass of the species in units of ion mass
  double mass[] = { 1. , 1./1836. , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

  // Density of the species in units of proton density
  double density[] = { 1. , 1. , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

  // Drift speed of the species in units of proton Alfven speed
  double vdrift[] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

  for (int d = 0; d < 10; d++) {
    params.beta[d] = beta[d];
    params.aniso[d] = aniso[d];
    params.charge[d] = charge[d];
    params.mass[d] = mass[d];
    params.density[d] = density[d];
    params.vdrift[d] = vdrift[d];
  }

  // Angle of propagation (in radian)
  params.theta = 0.001*M_PI/180.;

  // Alfven speed divided by speed of light 
  params.vAc=1.e-4;

  // Range and number of steps for NHDS_counter and NHDS_starter:
  params.omegarange[0] = 0.001;
  params.omegarange[1] = 1.1;
  params.gammarange[0] = 1.e-5;
  params.gammarange[1] = 1.e-1;


  for ( int i = 0; i < params.nr_kinds; i++) {
    params.Omega[i]=params.charge[i]/params.mass[i];
    params.ell[i]=sqrt(params.mass[i]/(params.density[i]*params.charge[i]*params.charge[i]));
    params.vtherm[i]=sqrt(params.beta[i]/(params.density[i]*params.mass[i]));
  }
 
}

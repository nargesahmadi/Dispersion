#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <complex.h>
#include "dispersion.h"

//--------------------------------------------------------------------------------------

void
epsilon(double complex w, double k, double complex D[3][3])
{
  double complex ep[3][3], chi[3][3];

  double kpar = k * cos(params.theta);
  double kperp = k * sin(params.theta);

  for ( int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {

      ep[j][k] = 0.0;
      
    }
  }

  for ( int kind = 0; kind < params.nr_kinds; kind++) {
    calc_chi(w,k,kind,chi);
    for ( int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {

	ep[j][k] = ep[j][k] + chi[j][k];
      
      }
    }
  }


  // vAc is the Alfven speed in units of speed of light

  for ( int i = 0; i < 3; i++) {
 
    ep[i][i] = ep[i][i] + params.vAc*params.vAc*(w*w);
  }

  for ( int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {

      D[j][k] = ep[j][k];
      
    }
  }

  D[0][0] = D[0][0] - kpar*kpar;
  D[0][1] = D[0][1];
  D[0][2] = D[0][2] + kpar*kperp;
  D[1][0] = D[1][0];
  D[1][1] = D[1][1] - (kperp*kperp+kpar*kpar);
  D[1][2] = D[1][2];
  D[2][0] = D[2][0] + kpar*kperp;
  D[2][1] = D[2][1];
  D[2][2] = D[2][2] - kperp*kperp;

}

//--------------------------------------------------------------------------------------

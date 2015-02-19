#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include "dispersion.h"


//--------------------------------------------------------------------------------------

double complex
determinant(double complex w, double k)
{
  double complex D[3][3], det;

  epsilon(w,k,D);

  det = D[0][0]*D[1][1]*D[2][2];
  det = det + D[0][1]*D[1][2]*D[2][0] + D[0][2]*D[1][0]*D[2][1];
  det = det - D[2][0]*D[1][1]*D[0][2] - D[2][1]*D[1][2]*D[0][0] - D[2][2]*D[1][0]*D[0][1];

  return det;
 
}

//--------------------------------------------------------------------------------------

double
absval(double complex z)
{
  double realpart;
  double imagpart;
  double dbval;

  realpart = creal(z);
  imagpart = cimag(z);

  dbval=sqrt(realpart*realpart+imagpart*imagpart);

  return dbval;

}

//--------------------------------------------------------------------------------------

double complex
newton(double complex w, double k, double complex *Dwk)
{

  double complex w_old, detw_old, detw, dw;
  int iter;
  bool go_for_newton;
  
  w_old = w-1.e-5-1.e-5*I;

  //newton method
  detw = 0;
  detw_old = determinant(w_old,k);
  iter = 0;
  go_for_newton = true;

  while ( (iter < params.numiter) && (go_for_newton) ) {
    iter = iter + 1;
    detw = determinant(w,k);
    if ( (absval(detw-detw_old) < 1e-50) || (absval(detw) < params.det_D_threshold ) ) {
      dw = 0.0;
      go_for_newton = false;
    } else {
      dw = detw*(w-w_old)/(detw-detw_old);
    }
    w_old = w;
    w = w - dw;
    detw_old = detw;
  }

  *Dwk = determinant(w,k);

  return w;
}

//--------------------------------------------------------------------------------------


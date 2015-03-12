#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include "dispersion.h"
#include "config.h"

#define wofz_F77 F77_FUNC_(wofz, WOFZ)
#define bessel_FC FC_FUNC_(bessel, BESSEL)

extern void wofz_(double *, double *, double *, double *, bool *);
extern void bessel_(int *, double *, double *, double *);
//extern void wofz_(double *XI, double *YI, double *U, double *V, bool *flag);
//extern void bessel_(int *n, double *lambda, double *BInz, double *dBInzdz);

//--------------------------------------------------------------------------------------
double complex
dispfunc(double complex zeta, bool kpos)
{
  double U,V;
  bool flag;
  double XI = creal(zeta);
  double YI = cimag(zeta);
  double complex DF;
  
  if ( kpos ) {
    wofz_F77(&XI, &YI, &U, &V,&flag);
    DF = 1.* I * sqrt(M_PI)* ( U + I * V);
  } else {
    double XXI = -XI;
    double YYI = -YI;
    wofz_F77(&XXI, &YYI, &U, &V,&flag);
    DF = -1.* I * sqrt(M_PI)* ( U + I * V);
  }
  return DF;
}
//--------------------------------------------------------------------------------------

void
calc_Y(double complex w, double k, int kind, int n, double complex Y[3][3])
{
  double complex zeta, An, Bn, resfac;
  double BInz, dBInzdz;
  double lambda, kpar, kperp;
  bool kpos;
  double *Omega = params.Omega;
  double *vdrift = params.vdrift;
  double *aniso = params.aniso;
  double *vtherm = params.vtherm;
  
  kpar = k * cos(params.theta);
  kperp = k * sin(params.theta);
  
  // positive k
  kpos= true;
  if (kpar < 0. ) kpos=false;
 

  zeta=(w-kpar*vdrift[kind]-1.*n*Omega[kind])/(kpar*vtherm[kind]);
  resfac=w-kpar*vdrift[kind]-1.*n*Omega[kind];
  lambda=0.5*(kperp*vtherm[kind]/Omega[kind])*(kperp*vtherm[kind]/Omega[kind])*aniso[kind];

  
  // calculating An and Bn
  An=(aniso[kind]-1.);
  An=An+(1./(kpar*vtherm[kind])) *( aniso[kind]*resfac + 1.*n*Omega[kind])*dispfunc(zeta,kpos);

  Bn=(aniso[kind]*(w-n*Omega[kind])-(kpar*vdrift[kind]-1.*n*Omega[kind]))/(kpar);
  Bn=Bn+( (w-1.*n*Omega[kind])*(aniso[kind]*resfac+1.*n*Omega[kind])/(kpar*kpar*vtherm[kind]) )*dispfunc(zeta,kpos);


  // calculating the modified Bessel Function and its derivative
  bessel_FC(&n, &lambda, &BInz, &dBInzdz);
  //printf("BInz = %g\n", BInz);
  
  // calculating the Y[j][k]
  Y[0][0] = 1. * (n*n) * BInz * An / lambda ;
  Y[0][1] = -I * n * (BInz - dBInzdz) * An;
  Y[0][2] = kperp *n *BInz* Bn/ ( Omega[kind] * lambda);
  Y[1][0] = I * n * (BInz - dBInzdz) * An; // -Y[1][2]
  Y[1][1] = (1.* (n*n) * BInz / lambda + 2.* lambda * BInz - 2.* lambda * dBInzdz)*An;
  Y[1][2] = I * kperp * (BInz - dBInzdz) * Bn / Omega[kind];
  Y[2][0] = kperp  *BInz * n * Bn / (Omega[kind] * lambda);
  Y[2][1] = -I * kperp * (BInz - dBInzdz) * Bn / Omega[kind];
  Y[2][2] = 2. * (w - 1.* n * Omega[kind]) * BInz * Bn / (kpar * vtherm[kind]*vtherm[kind]*aniso[kind]);
  
}
//--------------------------------------------------------------------------------------

void
calc_chi(double complex w, double k, int kind, double complex chi[3][3])
{
  double complex Y[3][3], Ynew[3][3];
  double lambda, kperp, kpar;
  int nmaxrun;
  double BInz, dBInzdz;
  bool Bessel_run;

  double *Omega = params.Omega;
  double *vdrift = params.vdrift;
  double *aniso = params.aniso;
  double *vtherm = params.vtherm;
  double *ell = params.ell;
  
  kpar = k * cos(params.theta);
  kperp = k * sin(params.theta);
  lambda=0.5*(kperp*vtherm[kind]/Omega[kind])*(kperp*vtherm[kind]/Omega[kind])*aniso[kind];
  
  
  for ( int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {

      Y[j][k] = 0.0;
    }
  }

  //FIXME; find a smarter way to determine the nmax for Bessels

  // determine maximum n for Besselfunction:
  nmaxrun = params.nmax;
  int m = 0;
  Bessel_run = true;
  while (Bessel_run) {
    bessel_FC(&m, &lambda, &BInz, &dBInzdz);
    if ((m <= params.nmax) && ( BInz < params.Bessel_zero)) {
      nmaxrun = m;
      Bessel_run = false;
    }
    m = m+1;
  }

  for (int n = -nmaxrun; n < nmaxrun + 1; n++) {
    calc_Y(w,k,kind,n,Ynew);
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
	Y[j][k] = Y[j][k] + exp(-lambda) * Ynew[j][k];
      }
    }
  }
  
  chi[0][0] = Y[0][0]/(ell[kind]*ell[kind]);
  chi[0][1] = Y[0][1]/(ell[kind]*ell[kind]);
  chi[0][2] = Y[0][2]/(ell[kind]*ell[kind]);
  chi[1][0] = Y[1][0]/(ell[kind]*ell[kind]);
  chi[1][1] = Y[1][1]/(ell[kind]*ell[kind]);
  chi[1][2] = Y[1][2]/(ell[kind]*ell[kind]);
  chi[2][0] = Y[2][0]/(ell[kind]*ell[kind]);
  chi[2][1] = Y[2][1]/(ell[kind]*ell[kind]);
  chi[2][2] = w*2.*vdrift[kind]/(ell[kind]*ell[kind]*kpar*vtherm[kind]*vtherm[kind]*aniso[kind])+Y[2][2]/(ell[kind]*ell[kind]);
}
//--------------------------------------------------------------------------------------

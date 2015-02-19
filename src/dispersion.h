#include <complex.h>

void epsilon(double complex w,double k, double complex D[3][3]);
void calc_chi(double complex w, double k, int kind, double complex chi[3][3]);
double complex newton(double complex w, double k, double complex *Dwk);

struct parameters {
  double density[10];             // Density of the species in units of proton density
  double charge[10];              // Charge of the species in units of the first ion charge
  double mass[10];                // Mass of the species in units of ion mass
  double aniso[10];               // Temperature anisotropy (Tperp/Tparallel)
  double beta[10];                // Parallel beta of the species
  double Omega[10];
  double ell[10];
  double vtherm[10];
  double vdrift[10];              // Drift speed of the species in units of proton Alfven speed
  double vAc;                     // Alfven speed divided by speed of light 
  double theta;                   // Angle of propagation (in degrees)
  double det_D_threshold;         // Threshold for the determinant of the dispersion tensor:
                                  // If det D <= det_D_threshold, the Newton iteration will be stopped
  double Bessel_zero;             // If I_n is less than this value, higher n are neglected:
  double omegarange[2];
  double gammarange[2];
  int nr_kinds;                   // Number of species 
  int numiter;                    // Maximum number of iterations in the Newton method
  int nmax;                       // Maximum of sum in Bessel function
                                  // can be very low (e.g., 3) for quasi-parallel propagation
  //int numb_counter;              // Range and number of steps for NHDS_counter and NHDS_starter:
  //int numb_starter;
};
void init_parameters(void);

extern struct parameters params;

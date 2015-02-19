
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "dispersion.h"

struct parameters params;

int
main(int argc, char **argv)
{
  double k ;
  double dk = 0.01;
  double complex w;
  double complex Dwk = 0;  // the quality of D(w,k) = 0
 

  // initialize parameters
  init_parameters();

  
  //w = 1e-20 + 1e-5* I;   // Starting frequency
  w = 1.e-1 + 1e-3*I;

  FILE *file = fopen("output.asc", "w");
  fprintf(file, "#k  wr  wi  Dr  Di\n");    

  for ( int i = 5; i < 101; i++ ) {

    k = dk * i;
    w = newton(w,k,&Dwk);

    printf("%6.4g %16.6g %16.6g %12.6g %12.6g\n", k, creal(w),cimag(w), creal(Dwk), cimag(Dwk));    
    fprintf(file, "%6.4g %16.6g %16.6g %12.6g %12.6g\n", k, creal(w),cimag(w), creal(Dwk), cimag(Dwk));    
  }

  fclose(file);
  return 0;
} // end main

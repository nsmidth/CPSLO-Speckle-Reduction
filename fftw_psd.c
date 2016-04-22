# include <stdlib.h>
# include <stdio.h>
# include <fftw3.h>

#define IMG_SIZE 512

// Compile with : 
//  


double * psd(double in[]);

double * psd(double in[]) {
  // Counters
  int i;
  int j;

  // Img Dimensions
  int nx = IMG_SIZE;
  int ny = IMG_SIZE;
  // PSD has symmetry about the center, 
  //  only need to calculate half the image
  int nyh = ( ny / 2 ) + 1;;

  // Allocate memory for output/plans
  fftw_complex *out;
  static double out_psd[IMG_SIZE*(IMG_SIZE/2+1)];
  fftw_plan plan_forward;

/*
  Create the output array OUT, which is of type FFTW_COMPLEX,
  and of a size NX * NYH
*/
  out = fftw_malloc ( sizeof ( fftw_complex ) * nx * nyh );

/*
  Set up and execute FFTW Plan
*/
  plan_forward = fftw_plan_dft_r2c_2d ( nx, ny, in, out, FFTW_ESTIMATE );

  fftw_execute ( plan_forward );
/*
  Calculate PSD of complex output
*/
  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < nyh; j++ )
    {
      out_psd[i*nyh+j] = ((out[i*nyh+j][0]*out[i*nyh+j][0]) + (out[i*nyh+j][1]*out[i*nyh+j][1]));    
    }
  }


/*
  Free up the allocated memory and return PSDs
*/
  fftw_destroy_plan ( plan_forward );

  fftw_free ( out );

  return out_psd;
}

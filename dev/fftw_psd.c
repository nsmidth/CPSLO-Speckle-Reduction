# include <stdlib.h>
# include <stdio.h>
# include <fftw3.h>
# include <time.h>

//#define PRINT_DATA

#define IMG_SIZE 512

// Compile with : 
//  gcc -o fftw_psd.so -shared -fPIC fftw_psd.c


float * psd(float * in_img);
float * test(float * in_img);

// Main tests the timing of PSD function
void main(void) {
  int i;
  float *in;
  float *out_psd;
  int j;
  int nx = IMG_SIZE;
  int ny = IMG_SIZE;
  int nyh = ( ny / 2 ) + 1;
  unsigned int seed = 123456789;
  unsigned int max_input = 1024;

/*
  Create the input array, an NX by NY array of floats.
*/
  in = ( float * ) malloc ( sizeof ( float ) * nx * ny );
  out_psd = ( float * ) malloc ( sizeof (float ) * nx * nyh );
/*
  srand ( seed );

  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {
      in[i*ny+j] = rand ( )%max_input;
    }
  }
*/  
  
  in[4*ny+5] = 1;
  in[4*ny+3] = 1;

/*
  Save start time
*/
  clock_t tic = clock();

/* 
  Execute
*/
  out_psd = psd(in);

/* 
Save FFT time
*/
  clock_t toc = clock();

/*
  Print Data
*/

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Demonstrate FFTW3 on a %d by %d array of real data.\n",
    nx, ny );
  printf ( "\n" );
  printf ( "  Transform data to FFT coefficients.\n" );
  printf ( "  Backtransform FFT coefficients to recover data.\n" );
  printf ( "  Compare recovered data to original data.\n" );

  printf("\n\n");
  printf("  Total Computation Time : \n");
  printf("  Elapsed: %f seconds\n", (float)(toc - tic) / CLOCKS_PER_SEC);

  #ifdef PRINT_DATA
      printf ( "\n" );
      printf ( "  Input Data:\n" );
      printf ( "\n" );

      for ( i = 0; i < nx; i++ )
      {
        for ( j = 0; j < ny; j++ )
        {
          printf ( "  %4d  %4d  %12f\n", i, j, in[i*ny+j] );
        }
      }

      printf ( "\n" );
      printf ( "  Output PSD:\n" );
      printf ( "\n" );

      for ( i = 0; i < nx; i++ )
      {
        for ( j = 0; j < nyh; j++ )
        {
          printf ( "  %4d  %4d  %12f  \n", 
          i, j, out_psd[i*nyh+j]);
        }
      }
  #endif


/*
  Free up the allocated memory.
*/

  free( in );
//  free( out_psd );
}

float * psd(float * in_img) {
  // Counters
  int i;
  int j;

  // Time holders
  clock_t tic,toc;

  // Img Dimensions
  int nx = IMG_SIZE;
  int ny = IMG_SIZE;
  // PSD has symmetry about the center, 
  //  only need to calculate half the image
  int nyh = ( ny / 2 ) + 1;
  double psd;

  // Allocate memory for fft output/plans
  fftw_complex *out;
  double in_img_double[IMG_SIZE*IMG_SIZE];
  static float out_psd[IMG_SIZE*(IMG_SIZE/2+1)];
  fftw_plan plan_forward;

  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {
       in_img_double[i*ny+j] = (double)in_img[i*ny+j] ;
    }
  }  

/*
  Create the output array OUT, which is of type FFTW_COMPLEX,
  and of a size NX * NYH
*/
  out = fftw_malloc ( sizeof ( fftw_complex ) * nx * nyh );
/*
  Set up and execute FFTW Plan
*/
  plan_forward = fftw_plan_dft_r2c_2d ( nx, ny, in_img_double, out, FFTW_ESTIMATE );

  fftw_execute ( plan_forward );
  
/*
  Calculate PSD of complex output
*/
  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < nyh; j++ )
    {
      //out_psd[i*nyh+j] = 3.0;
      psd = ((out[i*nyh+j][0]*out[i*nyh+j][0]) + (out[i*nyh+j][1]*out[i*nyh+j][1]));    
//      printf("Current index = %d, ", i*nyh+j);
//      printf("Present PSD (double) = %f, ", psd);
//      printf("Present PSD (float) = %f, ", (float)psd);
//      printf("\n");
      out_psd[i*nyh+j] = (float)psd;
    }
  }

/*
  Free up the allocated memory and return PSDs
*/  
  
  fftw_destroy_plan ( plan_forward );
  fftw_free ( out );

  return out_psd;
}

float * test(float * in_img) {
  // Counters
  int i;
  int j;
  static float out[IMG_SIZE*(IMG_SIZE/2+1)];
  
  for (i=0;i<IMG_SIZE;i++){
    for (j=0;j<((IMG_SIZE/2)+1);j++) {
      out[i*(IMG_SIZE/2+1)+j] = in_img[i*(IMG_SIZE/2+1)+j];
    }
  }
  
  return out;
  
}
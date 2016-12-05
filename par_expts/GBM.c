/*
 * Example SDE code, implementing the Euler-Maruyama scheme for 
 * Geometric Brownian Motion using MKL/VSL random number generation
 * and OpenMP vectorisation and parallelisation
 *
 * author: Mike Giles, based on previous code written by 
 *         David J. Warne
 *         School of Mathematical Sciences 
 *         Queensland University of Technology
 *
 * Date: Nov 23 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl.h>
#include <mkl_vsl.h>
#include <memory.h>
#include <omp.h>

// vector length for vectorisation
#define VECTOR_LENGTH 4   

// macros for rounding x up or down to a multiple of y
#define ROUND_UP(x, y) ( ( ((x) + (y) - 1) / (y) ) * (y) )
#define ROUND_DOWN(x, y) ( ((x) / (y)) * (y) )

/* each OpenMP thread has its own VSL RNG and storage */
#define RV_BYTES 65536
double           *dW;
VSLStreamStatePtr stream;
#pragma omp threadprivate(stream, dW)

void pathcalc(double,double,double,double,double, int,int,
              double *,double *);

int main(int argc, char **argv)
{
  double T=1.0, X0=1.0, mu=0.05, sigma=0.2, sum1=0.0, sum2=0.0, dt;
  int    M, N, N2;

  M  = 20;      /* number of timesteps */
  N  = 5120000; /* total number of MC samples */

  dt  = T / (double)M;

#pragma omp parallel shared(T,X0,mu,sigma,dt,M,N)	\
                     reduction(+:sum1,sum2)
  {
    double sum1_t = 0.0, sum2_t = 0.0;
    int    num_t  = omp_get_num_threads();
    int    tid    = omp_get_thread_num();

    printf("tid=%d, creating RNG generator and allocating memory \n",tid); 

    /* create RNG, then give each thread a unique skipahead */
    vslNewStream(&stream, VSL_BRNG_MRG32K3A,1337);
    long long skip = ((long long) (tid+1)) << 48;
    vslSkipAheadStream(stream,skip);

    dW = (double *)malloc(RV_BYTES);

    int N3 = ROUND_UP(((tid+1)*N)/num_t,VECTOR_LENGTH)
           - ROUND_UP(( tid   *N)/num_t,VECTOR_LENGTH);

    pathcalc(T,X0,mu,sigma,dt, M, N3, &sum1_t, &sum2_t);
    sum1 += sum1_t;
    sum2 += sum2_t;
  }
    
  printf("Exact solution E[X_T] = %g\n",X0*exp(mu*T));
  printf("Monte Carlo estimate  = %g +/- %g \n",sum1/N,
         3.0*sqrt((sum2/N-(sum1/N)*(sum1/N))/N));
  printf("Reminder: Monte Carlo estimate has discretisation bias\n");

  /* delete generator and storage */
#pragma omp parallel 
  {
    vslDeleteStream(&stream);
    free(dW);
  }
}

void pathcalc(double T,double X0,double mu,double sigma,double dt,
              int M, int N, double *sum1_t, double *sum2_t) {

  double sum1=0.0, sum2=0.0;

  /* work out max number of paths in a group */
  int bytes_per_path = M*sizeof(double);
  int N2 = ROUND_DOWN(RV_BYTES/bytes_per_path,VECTOR_LENGTH);

  /* loop over all paths in groups of size N2 */
  for (int n0=0; n0<N; n0+=N2) {
    /* may have to reduce size of final group */
    if (N2>N-n0) N2 = N-n0;

    /* generate required random numbers for this group */
    vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2,
                  stream,M*N2,dW,0,sqrt(dt));

    /* loop over paths within group in increments of VECTOR_LENGTH */
    for (int n1=0; n1<N2; n1+=VECTOR_LENGTH) {
      int offset = n1*M; /* number of random numbers already used */

      /* vectorised path calculation */
#pragma omp simd reduction(+:sum1,sum2)
      for (int n2=0; n2<VECTOR_LENGTH; n2++) {
        double X = X0;
	
        for (int m=0; m<M; m++) {
          double delW  = dW[offset+m*VECTOR_LENGTH+n2]; 
          X = X*(1.0 + mu*dt + sigma*delW);
        }

        sum1 += X;
        sum2 += X*X;
      }
    }
  }

  *sum1_t = sum1;
  *sum2_t = sum2;
}

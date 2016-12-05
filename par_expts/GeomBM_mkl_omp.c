/**
 * Example SDE code. Implements the Euler-Maruyama Scheme for
 * the Geometric Brownian Motion Stochastic Differential Equation.
 * exploits vector operations via MKL functions and OpenMP pragmas
 *
 * author: David J. Warne (david.warne@qut.edu.au)
 *         School of Mathematical Sciences 
 *         Queensland University of Technology
 *
 * Date: Oct 27 2016
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <mkl.h>
#include <mkl_vsl.h>

#define NTHREADS 2

int main(int argc, char **argv)
{
    double T, X_0, *X_t, mu, sigma, dt, *dW, E,*temp,*temp2;
    unsigned int N, N_t,N_par,i,k;
    VSLStreamStatePtr stream[NTHREADS];

    N_t = 10000;   
    N = (unsigned int)atoi(argv[1]);    

    N_par = N/NTHREADS; /*for now assume N is a multiple of NTHREADS*/

    X_0 = 1.0; mu = 1.0; sigma = 0.2; T = 1.0;    
    dt = T/(double)N_t;
    X_t = (double *)malloc(N*sizeof(double));
    temp = (double *)malloc(N*sizeof(double));
    temp2 = (double *)malloc(N*sizeof(double));
    dW = (double *)malloc(N*sizeof(double));

    /*Basic generator just to be consistent with non-MKL code*/
    vslNewStream(&(stream[0]), VSL_BRNG_MCG31,1337);

    /*split into independent subsequences using leap frog*/
    for (k=1;k<NTHREADS;k++){
        vslCopyStream(&(stream[k]), stream[0]);
    }
    for (k=0;k<NTHREADS;k++){
        vslLeapfrogStream(stream[k], k,NTHREADS);
    }


    E = 0;
    for (i=0;i<N;i++){
        X_t[i] = X_0;
    }
    omp_set_num_threads(NTHREADS);
    #pragma omp parallel for schedule(static)
    for (k=0;k<NTHREADS;k++){ 
        unsigned int offset;
        offset = k*N_par;
        /* Compute N Euler-Maruyama simulations using MKL*/
        for (int t=0;t<N_t;t++){
            /*dW[i] ~ N(0,sqrt(dt))*/
            vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,stream[k],N_par,dW+offset,0,sqrt(dt));
            /* temp2[i] = dW[i]*sigma + mu*dt */
            vdLinearFrac(N_par,dW+offset,temp+offset,sigma,mu*dt,0,1,temp2+offset);
            /*temp[i] = X_t[i]*temp2[i]*/
            vdMul(N_par,temp2+offset,X_t+offset,temp+offset);
            /*temp2[i] = X_t[i] + temp[i]*/
            vdAdd(N_par,temp+offset,X_t+offset,temp2+offset);
            memcpy(X_t+offset,temp2+offset,N_par*sizeof(double));
        }
    }
        
    for (i=0;i<N;i++){
        E += X_t[i];
    }
    E /= (double)N;
    
    printf("Exact solution E[X_T] = %g\n",X_0*exp(mu*T));
    printf("Monte Carlo estimate = %g\n",E);
    for (k=0;k<NTHREADS;k++){
        vslDeleteStream(&(stream[k]));
    }
}


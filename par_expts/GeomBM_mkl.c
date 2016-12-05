/**
 * Example SDE code. Implements the Euler-Maruyama Scheme for
 * the Geometric Brownian Motion Stochastic Differential Equation.
 * exploits vector operations via MKL fuctions
 *
 * author: David J. Warne (david.warne@qut.edu.au)
 *         School of Mathematical Sciences 
 *         Queensland University of Technology
 *
 * Date: Oct 13 2016
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl.h>
#include <mkl_vsl.h>
#include <memory.h>

int main(int argc, char **argv)
{
    double T, X_0, *X_t, mu, sigma, dt, *dW, E,*temp,*temp2;
    unsigned int N, N_t,i,t;
    VSLStreamStatePtr stream;

    N_t = 10000;   
    N = (unsigned int)atoi(argv[1]);    

    X_0 = 1.0; mu = 1.0; sigma = 0.2; T = 1.0;    
    dt = T/(double)N_t;
    X_t = (double *)malloc(N*sizeof(double));
    temp = (double *)malloc(N*sizeof(double));
    temp2 = (double *)malloc(N*sizeof(double));
    dW = (double *)malloc(N*sizeof(double));

    /*Basic generator just to be consistent with non-MKL code*/
    vslNewStream(&stream, VSL_BRNG_MCG31,1337);

    E = 0;
    for (i=0;i<N;i++){
        X_t[i] = X_0;
    }
    
    /* Compute N Euler-Maruyama simulations using MKL*/
    for (t=0;t<N_t;t++){
        /*dW[i] ~ N(0,sqrt(dt))*/
        vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,stream,N,dW,0,sqrt(dt));
        /* temp2[i] = dW[i]*sigma + mu*dt */
        vdLinearFrac(N,dW,temp,sigma,mu*dt,0,1,temp2);
        /*temp[i] = X_t[i]*temp2[i]*/
        vdMul(N,temp2,X_t,temp);
        /*temp2[i] = X_t[i] + temp[i]*/
        vdAdd(N,temp,X_t,temp2);
        memcpy(X_t,temp2,N*sizeof(double));
    }
        
    for (i=0;i<N;i++){
        E += X_t[i];
    }
    E /= (double)N;
    
    printf("Exact solution E[X_T] = %g\n",X_0*exp(mu*T));
    printf("Monte Carlo estimate = %g\n",E);
}


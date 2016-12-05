/**
 * Example SDE code. Implements the Euler-Maruyama Scheme for
 * the Geometric Brownian Motion Stochastic Differential Equation.
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

/* generate normal random variates*/
double normrnd(double, double);

int main(int argc, char **argv)
{
    double T, X_0, X_t, mu, sigma, dt, dW, E;
    unsigned int N, N_t,i,t;
    X_0 = 1.0; mu = 1.0; sigma = 0.2; T = 1.0;    
    N_t = 10000;   
    N = (unsigned int)atoi(argv[1]);    
    dt = T/(double)N_t;

    /*compute expectation*/
    E = 0;
    for (i=0;i<N;i++){
        /*Euler-Maruyama simulation*/
        X_t = X_0;
        for (t=0;t<N_t;t++){
            dW = normrnd(0,sqrt(dt));
            X_t += X_t*(mu*dt + sigma*dW);
        }
        E += X_t;
    }
    E /= (double)N;
    
    printf("Exact solution E[X_T] = %g\n",X_0*exp(mu*T));
    printf("Monte Carlo estimate = %g\n",E);
}

/**
 *
 * Generates n i.i.d random variates Z_i ~ N(mu,sigma)
 *
 * parameters:
 *
 * n - number of variates 
 * Z - array to store variates
 * mu - mean of the distribution,
 * sigma - the standard deviation.
 */
double normrnd(double mu, double sigma)
{
    int i;
    double U_1, U_2;
    double r,theta;
    /*uses the Box-Muller transform*/
    U_1 = ((double)rand())/((double)RAND_MAX);
    U_2 = ((double)rand())/((double)RAND_MAX);
    r = sigma*sqrt(-2.0*log(U_1));
    theta = 2.0*M_PI*U_2; 
    return r*cos(theta) + mu; 
}


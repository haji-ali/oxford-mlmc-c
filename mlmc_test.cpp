/*

  mlmc_test(mlmc_l, M, N,L, N0,Eps,Lmin,Lmax, fp)

  multilevel Monte Carlo test routine

   mlmc_l(l,N,sums)     low-level routine

   inputs:  l = level
            N = number of paths

   output: sums[0] = sum(cost)
           sums[1] = sum(Pf-Pc)
           sums[2] = sum((Pf-Pc).^2)
           sums[3] = sum((Pf-Pc).^3)
           sums[4] = sum((Pf-Pc).^4)
           sums[5] = sum(Pf)
           sums[6] = sum(Pf.^2)

   M      = refinement cost factor (if zero, gamma is estimated)

   N      = number of samples for convergence tests
   L      = number of levels for convergence tests

   N0     = initial number of samples
   Eps    = desired accuracy array (terminated by value 0)
   Lmin   = minimum level of refinement
   Lmax   = maximum level of refinement

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <mlmc.cpp>

// https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
// variadic macro to print to both file and stdout
#define PRINTF2(fp, ...) {printf(__VA_ARGS__);fprintf(fp,__VA_ARGS__);}


void mlmc_test(void (*mlmc_l)(int, int, double *), int M,int N,int L,
               int N0, float *Eps, int Lmin, int Lmax, FILE *fp) {

//
// first, convergence tests
//

  PRINTF2(fp,"\n");
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"*** Convergence tests, kurtosis, telescoping sum check ***\n");
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"\n l   ave(Pf-Pc)    ave(Pf)   var(Pf-Pc)    var(Pf)");
  PRINTF2(fp,"    kurtosis     check        cost \n-------------------------");
  PRINTF2(fp,"------------------------------------------------------------\n");

  double sums[7];
  float *cost = (float *)malloc((L+1)*sizeof(float));
  float *del1 = (float *)malloc((L+1)*sizeof(float));
  float *del2 = (float *)malloc((L+1)*sizeof(float));
  float *var1 = (float *)malloc((L+1)*sizeof(float));
  float *var2 = (float *)malloc((L+1)*sizeof(float));
  float *chk1 = (float *)malloc((L+1)*sizeof(float));
  float *kur1 = (float *)malloc((L+1)*sizeof(float));

  for (int l=0; l<=L; l++) {
    mlmc_l(l,N,sums);

    for (int m=0; m<7; m++) sums[m] = sums[m]/N;

    if (M>0)
      cost[l] = powf((float)M,(float)l);
    else
      cost[l] = sums[0];
    del1[l] = sums[1];
    del2[l] = sums[5];
    var1[l] = fmax(sums[2]-sums[1]*sums[1], 1e-10);
    var2[l] = fmax(sums[6]-sums[5]*sums[5], 1e-10);

    kur1[l]  = (      sums[4]
                - 4.0*sums[3]*sums[1]
                + 6.0*sums[2]*sums[1]*sums[1]
                - 3.0*sums[1]*sums[1]*sums[1]*sums[1] )
             / (var1[l]*var1[l]);

    if (l==0)
      chk1[l] = 0.0f;
    else
      chk1[l] = sqrtf((float) N) * 
                fabsf(  del1[l]  +       del2[l-1]  -       del2[l] )
         / (3.0f*(sqrtf(var1[l]) + sqrtf(var2[l-1]) + sqrtf(var2[l])));

    PRINTF2(fp,"%2d   %8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e \n",
    l,del1[l],del2[l],var1[l],var2[l],kur1[l],chk1[l],cost[l]);
  }

//
// print out a warning if kurtosis or consistency check looks bad
//

  if (kur1[L] > 100.0f) {
    PRINTF2(fp,"\n WARNING: kurtosis on finest level = %f \n",kur1[L]);
    PRINTF2(fp," indicates MLMC correction dominated by a few rare paths; \n");
    PRINTF2(fp," for information on the connection to variance of sample variances,\n");
    PRINTF2(fp," see http://mathworld.wolfram.com/SampleVarianceDistribution.html \n");
  }

  float max_chk = 0.0f;
  for (int l=0; l<=L; l++) max_chk = fmaxf(max_chk,chk1[l]);
  if (max_chk > 1.0f) {
    PRINTF2(fp,"\n WARNING: maximum consistency error = %f \n",max_chk);
    PRINTF2(fp," indicates identity E[Pf-Pc] = E[Pf] - E[Pc] not satisfied \n");
  }

//
// use linear regression to estimate alpha, beta, gamma
//

  float alpha, beta, gamma, foo;
  float *x = (float *)malloc(L*sizeof(float));
  float *y = (float *)malloc(L*sizeof(float));

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = - log2f(fabsf(del1[l]));
  } 
  regression(L,x,y,alpha,foo);

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = - log2f(var1[l]);
  } 
  regression(L,x,y,beta,foo);

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = log2f(cost[l]);
  } 
  regression(L,x,y,gamma,foo);

  PRINTF2(fp,"\n******************************************************\n");
  PRINTF2(fp,"*** Linear regression estimates of MLMC parameters ***\n");
  PRINTF2(fp,"******************************************************\n");
  PRINTF2(fp,"\n alpha = %f  (exponent for MLMC weak convergence)\n",alpha);
  PRINTF2(fp," beta  = %f  (exponent for MLMC variance) \n",beta);
  PRINTF2(fp," gamma = %f  (exponent for MLMC cost) \n",gamma);

//
// second, mlmc complexity tests
//

  PRINTF2(fp,"\n");
  PRINTF2(fp,"***************************** \n");
  PRINTF2(fp,"*** MLMC complexity tests *** \n");
  PRINTF2(fp,"***************************** \n\n");
  PRINTF2(fp,"  eps       value   mlmc_cost   std_cost  savings     N_l \n");
  PRINTF2(fp,"--------------------------------------------------------- \n");
 
  int i=0;
  int   *Nl = (int *)malloc((Lmax+1)*sizeof(int));
  float *Cl = (float *)malloc((Lmax+1)*sizeof(float));

  while (Eps[i]>0) {
    float eps = Eps[i++];

    float P = mlmc(Lmin,Lmax,N0,eps,mlmc_l, alpha,beta,gamma, Nl,Cl);

    float std_cost = 0.0f, mlmc_cost = 0.0f, theta=0.25f;

    for (int l=0; l<=Lmax; l++) {
      if (Nl[l]>0) {
        // printf(" l, Cl, cost = %d  %f  %f \n",l,Cl[l],cost[l]);
        mlmc_cost += Nl[l]*Cl[l];
        if (l<=L)
          std_cost = var2[l]*Cl[l] / ((1.0f-theta)*eps*eps);
        else
          std_cost = var2[L]*Cl[l] / ((1.0f-theta)*eps*eps);
      }
    }

    PRINTF2(fp,"%.4f  %.4e  %.3e  %.3e  %7.2f ",
	    eps, P, mlmc_cost, std_cost, std_cost/mlmc_cost);
    for (int l=0; Nl[l]>0; l++) PRINTF2(fp,"%9d",Nl[l]);
    PRINTF2(fp,"\n");
  }
  PRINTF2(fp,"\n");
}

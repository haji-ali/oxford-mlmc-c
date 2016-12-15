/**
 * @file mlmc.c
 * @brief Source file for Multilevel Monte Carlo driver.
 * @details Defines a set of useful types and functions for typical
 * applications of Multilevel Monte Carlo.
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include "mlmc.h"

typedef int bool;

mlmc_options*
mlmc_create_option(double eps, fn_mlmc_sample_levels_t fn_mlmc_sample_levels)
{
    mlmc_options* opt = malloc(sizeof(mlmc_options));
    opt->eps = eps;
    opt->fn_mlmc_sample_levels = fn_mlmc_sample_levels;
    opt->Lmin = 2;
    opt->Lmax = 32;
    opt->N0 = 100;
    opt->alpha = opt->beta = opt->gamma = -1;
    opt->per_sample = 1;
    return opt;
}

mlmc_options*
mlmc_create_option_simple(double eps, fn_mlmc_sample_level_t fn_mlmc_sample_level)
{
    mlmc_options* opt = malloc(sizeof(mlmc_options));
    opt->eps = eps;
    opt->fn_mlmc_sample_level = fn_mlmc_sample_level;
    opt->fn_mlmc_sample_levels = NULL;
    opt->Lmin = 2;
    opt->Lmax = 32;
    opt->N0 = 100;
    opt->alpha = opt->beta = opt->gamma = -1;
    opt->per_sample = 1;
    opt->theta = 0.25;
    return opt;
}

void
mlmc_free_options(mlmc_options* opt)
{
    free(opt);
}

void
mlmc_free_output(mlmc_output* out)
{
    free(out->sums_moments);
    free(out->sums_work);
    free(out->Nl);
    free(out);
}

void regression(int N, double *x, double *y, double *a, double *b){

    float sum0=0.0f, sum1=0.0f, sum2=0.0f, sumy0=0.0f, sumy1=0.0f;

    for (int i=1; i<N; i++) {
        sum0  += 1.0f;
        sum1  += x[i];
        sum2  += x[i]*x[i];

        sumy0 += y[i];
        sumy1 += y[i]*x[i];
    }

    *a = (sum0*sumy1 - sum1*sumy0) / (sum0*sum2 - sum1*sum1);
    *b = (sum2*sumy0 - sum1*sumy1) / (sum0*sum2 - sum1*sum1);
}

void
mlmc_run_all_levels(unsigned int L,
    const mlmc_options* opt, const unsigned long long *dNl,
    mlmc_output *out)
{
    if (opt->fn_mlmc_sample_levels)
    {
        // Multi-level mode
        double *sums_moments = malloc(sizeof(double)*opt->Lmax*
                                      opt->per_sample*MLMC_MOMENTS_COUNT);
        double *sums_work = malloc(sizeof(double)*opt->Lmax);
        unsigned long long *dN_done = malloc(sizeof(unsigned int)*L*
                                       opt->per_sample*MLMC_MOMENTS_COUNT);
        unsigned long long *dN_todo = malloc(sizeof(unsigned int)*L*
                                       opt->per_sample*MLMC_MOMENTS_COUNT);
        for (unsigned int l=0;l<L;l++)
        {
            dN_done[l] = 0;
            dN_todo[l] = dNl[0];
        }

        while (1){
            opt->fn_mlmc_sample_levels(L, dN_todo, sums_moments,
                                       opt->Lmax*opt->per_sample
                                       *MLMC_MOMENTS_COUNT,
                                       sums_work, opt->user_data);
            unsigned int sum=0;
            for (unsigned int l=0;l<L;l++)
            {
                dN_done[l] += dN_todo[l];
                // TODO: Assimilate into output
                dN_todo[l] = dNl[l]>dN_done[l] ? dNl[l] - dN_done[l] : 0 ;
                sum += dN_todo[l];
            }
            if (sum == 0)
                break;
        }
        free(sums_moments);
        free(sums_work);
    }
    else
    {
        // Single level mode
        double *sums_moments = malloc(sizeof(double)*
                                      opt->per_sample*MLMC_MOMENTS_COUNT);
        double sums_work;
        for (unsigned int l=0;l<L;l++){
            unsigned long long dN_todo = dNl[l];
            unsigned long long dN_done = 0;
            while(dN_todo > 0){
                opt->fn_mlmc_sample_level(l, &dN_todo, sums_moments,
                                          opt->per_sample*MLMC_MOMENTS_COUNT,
                                          &sums_work, opt->user_data);
                dN_done += dN_todo;
                // TODO: Assimilate into output

                dN_todo = dNl[l]>dN_done ? dNl[l] - dN_done : 0 ;
            }
        }
        free(sums_moments);
    }
}

mlmc_output*
mlmc_run(const mlmc_options* opt)
{
    unsigned int Lmax = opt->Lmax;

    mlmc_output *out = malloc(sizeof(mlmc_output));
    out->result = 0;
    out->sums_moments = malloc(sizeof(double)*Lmax*opt->per_sample*
                               MLMC_MOMENTS_COUNT);
    out->Nl = malloc(sizeof(unsigned long long)*Lmax);
    out->sums_work = malloc(sizeof(double)*Lmax);

    bool fit_alpha = opt->alpha < 0;
    bool fit_beta = opt->beta < 0;
    bool fit_gamma = opt->gamma < 0;

    unsigned int L = opt->Lmin;
    bool converged = 0;

    unsigned long long *dNl = malloc(sizeof(long long) * (Lmax+1));
    double *El = malloc(sizeof(double) * (Lmax+1));
    double *Vl = malloc(sizeof(double) * (Lmax+1));
    double *Cl = malloc(sizeof(double) * (Lmax+1));
    double *x = malloc(sizeof(double) * (Lmax+1));
    double *y = malloc(sizeof(double) * (Lmax+1));

    double alpha = fabs(opt->alpha);
    double beta = fabs(opt->beta);
    double gamma = fabs(opt->gamma);

    bool diag = 1;

    for(unsigned int l=0; l<=Lmax; l++) {
        dNl[l]  = 0;
        Cl[l]   = powf(2.0f,(float)l*gamma);
        for(unsigned int n=0; n<opt->per_sample * MLMC_MOMENTS_COUNT; n++)
            out->sums_moments[n + l*opt->per_sample * MLMC_MOMENTS_COUNT] = 0.0;
    }

    for(int l=0; l<=opt->Lmin; l++)
        dNl[l] = opt->N0;

    while (!converged){
        //
        // update sample sums
        //
        mlmc_run_all_levels(L, opt, dNl, out);

        //
        // compute absolute average, variance and cost,
        // correct for possible under-sampling,
        // and set optimal number of new samples
        //

        double sum = 0.0f;
        for (unsigned int l=0; l <= Lmax; l++) {
            if (!out->Nl[l])     // No samples on this level
                continue;

            El[l] =
                fabs(out->sums_moments[MLMC_MOMENTS_COUNT* opt->per_sample*l]
                     / out->Nl[l]);
            Vl[l] =
                fmaxf(out->sums_moments[MLMC_MOMENTS_COUNT *opt->per_sample *l
                                + 1] / out->Nl[l] - El[l]*El[l], 0.0f);

            for (unsigned int i=1;i<opt->per_sample;i++){
                El[l] =
                    fmaxf(fabs(out->sums_moments[MLMC_MOMENTS_COUNT*
                                                 opt->per_sample*l + i]
                               / out->Nl[l]), El[l]);
                Vl[l] =
                    fmaxf(out->sums_moments[MLMC_MOMENTS_COUNT*opt->per_sample*l
                                            + i + 1] / out->Nl[l] - El[l]*El[l],
                          Vl[l]);
            }
            Cl[l] = out->sums_work[l] / out->Nl[l];

            // This is limit the convergence of the absolute error and variance
            // TODO: Check this
            if (l>1) {
                El[l] = fmaxf(El[l],  0.5f*El[l-1]/powf(2.0f,alpha));
                Vl[l] = fmaxf(Vl[l],  0.5f*Vl[l-1]/powf(2.0f,beta));
            }
            sum += sqrtf(Vl[l]*Cl[l]);
        }

        for (unsigned int l=0; l<=Lmax; l++) {
            dNl[l] = ceilf(fmaxf(0.0f,
                                 sqrtf(Vl[l]/Cl[l])*sum
                                 / ((1.0f-opt->theta)*opt->eps*opt->eps)
                                 - out->Nl[l] ));
        }

        //
        // use linear regression to estimate alpha, beta, gamma if not given
        //
        if (fit_alpha) {
            for (unsigned int l=1; l<=L; l++) {
                x[l-1] = l;
                y[l-1] = - log2f(El[l]);
            }
            regression(L,x,y,&alpha, &sum);
        }
        if (fit_beta) {
            for (unsigned int l=1; l<=L; l++) {
                x[l-1] = l;
                y[l-1] = - log2f(Vl[l]);
            }
            regression(L,x,y,&beta, &sum);
        }
        if (fit_gamma) {
            for (unsigned int l=1; l<=L; l++) {
                x[l-1] = l;
                y[l-1] = log2f(Cl[l]);
            }
            regression(L,x,y,&gamma, &sum);
        }

        //
        // if (almost) converged, estimate remaining error and decide
        // whether a new level is required
        //
        sum = 0.0;
        for (int l=0; l<=Lmax; l++)
            sum += fmaxf(0.0f, (double)dNl[l]-0.01*out->Nl[l]);

        if (sum == 0) {
            if (diag) printf(" achieved variance target \n");

            converged = 1;
            float rem = El[L] / (powf(2.0f,alpha)-1.0f);

            if (rem > sqrtf(opt->theta)*opt->eps) {
                if (L==Lmax)
                    out->result |= MLMC_RES_NOT_CONVERGED;
                else {
                    converged = 0;
                    L++;
                    Vl[L] = Vl[L-1] / powf(2.0f,beta);
                    Cl[L] = Cl[L-1] * powf(2.0f,gamma);

                    if (diag) printf(" L = %d \n",L);

                    sum = 0.0f;
                    for (int l=0; l<=L; l++) sum += sqrtf(Vl[l]*Cl[l]);
                    for (int l=0; l<=L; l++)
                        dNl[l] = ceilf( fmaxf( 0.0f,
                                               sqrtf(Vl[l]/Cl[l])*sum/((1.0f-opt->theta)*opt->eps*opt->eps)
                                               - out->Nl[l] ) );
                }
            }
        }
    }

    free(El);
    free(Vl);
    free(Cl);
    free(x);
    free(y);
    return out;
}

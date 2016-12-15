/**
 * @file mlmc.c
 * @brief Source file for Multilevel Monte Carlo driver.
 * @details Defines a set of useful types and functions for typical
 * applications of Multilevel Monte Carlo.
 *
 */
#include <memory.h>
#include <math.h>

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
    free(out);
}

void regression(int N, float *x, float *y, float *a, float *b){

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

mlmc_output*
mlmc_run(const mlmc_options* opt)
{
    mlmc_output *out = malloc(sizeof(mlmc_output));
    opt->sum = malloc(sizeof(double) * Lmax*opt->per_sample*2);

    bool fit_alpha = opt->alpha < 0;
    bool fit_beta = opt->beta < 0;
    bool fit_gamma = opt->gamma < 0;

    unsigned int L = Lmin;
    bool converged = 0;
    unsigned long long Nl[Lmax+1];
    unsigned long long dNl[Lmax+1];
    double El[Lmax+1], Vl[Lmax+1];
    double Cl[Lmax+1];
    double NlCl[Lmax+1];
    double x[Lmax+1], y[Lmax+1];

    double alpha = fabs(opt->alpha);
    double beta = fabs(opt->beta);
    double gamma = fabs(opt->gamma);

    bool diag = 1;

    for(unsigned int l=0; l<=Lmax; l++) {
        dNl[l]  = 0;
        Nl[l]   = 0;
        Cl[l]   = powf(2.0f,(float)l*gamma);
        NlCl[l] = 0.0f;
        for(int n=0; n<3; n++) suml[n][l] = 0.0;
    }

    for(int l=0; l<=Lmin; l++) dNl[l] = N0;

    while (!converged){
        //
        // update sample sums
        //
        mlmc_run_all_levels(opt, Lmax, dNl, suml, sizeof(suml));
        for (int l=0; l<=maxL; l++)
            out->Nl[l] += dNl[l];

        //
        // compute absolute average, variance and cost,
        // correct for possible under-sampling,
        // and set optimal number of new samples
        //

        double sum = 0.0f;
        for (int l=0; l<=maxL; l++) {
            if (!out->Nl[l])     // No samples on this level
                continue;

            El[l] = fabs(out->sums[2*l]/out->Nl[l]);
            Vl[l] = fmaxf(out->sums[2*l+1]/out->Nl[l] - El[l]*El[l], 0.0f);
            Cl[l] = out->sums[2*l+2] / out->Nl[l];

            // This is taken from the original code. It seems to limit the
            // convergence of the absolute error and variance
            if (l>1) {
                El[l] = fmaxf(El[l],  0.5f*El[l-1]/powf(2.0f,alpha));
                Vl[l] = fmaxf(Vl[l],  0.5f*Vl[l-1]/powf(2.0f,beta));
            }
            sum += sqrtf(Vl[l]*Cl[l]);
        }

        for (int l=0; l<=maxL; l++) {
            dNl[l] = ceilf(fmaxf(0.0f,
                                 sqrtf(Vl[l]/Cl[l])*sum/((1.0f-opt->theta)*opt->eps*opt->eps)
                                 - opt->Nl[l] ));
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
            sum += fmaxf(0.0f, (float)dNl[l]-0.01f*suml[0][l]);

        if (sum == 0) {
            if (diag) printf(" achieved variance target \n");

            converged = 1;
            float rem = El[L] / (powf(2.0f,alpha)-1.0f);

            if (rem > sqrtf(theta)*eps) {
                if (L==Lmax)
                    opt->results |= MLMC_RES_NOT_CONVERGED;
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
                                               sqrtf(Vl[l]/Cl[l])*sum/((1.0f-theta)*eps*eps)
                                               - opt->Nl[l] ) );
                }
            }
        }
    }

    return out;
}

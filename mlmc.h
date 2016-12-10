#ifndef __MLMC_H__
#define __MLMC_H__

typedef void (*fn_mlmc_sample_levels_t)(unsigned int L,
                                        const unsigned long long* M,
                                        double *sums,
                                        unsigned int sums_size,
                                        void *user_data);
typedef void (*fn_mlmc_sample_level_t)(unsigned int ell,
                                       unsigned long long M,
                                       double *sums,
                                       unsigned int sums_size,
                                       void *user_data);


typedef struct s_mlmc_options{
    unsigned int Lmin;  // Minimum number of levels
    unsigned int Lmax;  // Maximum number of levels

    unsigned long long N0;   // Initial number of samples to compute variance estimates

    double alpha;             // Weak convergence rate. Must be positive
    double beta;              // Variance convergence rate. Must be positive
    double gamma;             // Work rate. Must be positive

    // It will be difficult to support variable length vectors efficiently.
    // I think we should focus on fixed length vectors for now.
    unsigned int per_sample;  // Vector length for vector quantities

    // user supplied function to sample all levels
    fn_mlmc_sample_levels_t fn_mlmc_sample_levels;
    // user supplied function to sample a single level
    fn_mlmc_sample_level_t fn_mlmc_sample_level;

    void *user_data;
} mlmc_options;

typedef struct s_mlmc_output{
    unsigned int L;
    unsigned int moments;
    unsigned int per_moment;

    unsigned long long* M;   // Size L
    double* sums;            // Size L*moments*per_moments
}mlmc_output;

// commandline like options creation
mlmc_options*
mlmc_create_option_cli(int argc, char ** argv);

// explicit MLMC options creation
mlmc_options*
mlmc_create_option(unsigned int Lmin, unsigned int LMax,
                   unsigned int N0, double alpha, double beta, double gamma,
                   unsigned int per_sample,
                   fn_mlmc_sample_levels_t fn_mlmc_sample_levels);

mlmc_options*
mlmc_create_option_simple(unsigned int Lmin, unsigned int LMax,
                          unsigned int N0, double alpha, double beta, double gamma,
                          unsigned int per_sample,
                          fn_mlmc_sample_level_t fn_mlmc_sample_level);
void
mlmc_free_options(mlmc_options*);

mlmc_output*
mlmc_run(double eps, const mlmc_options*); // Options can be NULL for default

void
mlmc_free_output(mlmc_options*);

// This function calls a simple single-level sampler for every level.
void
mlmc_sample_levels(const mlmc_options*,
                   unsigned int L,
                   const unsigned long long* M,
                   double *sums,
                   unsigned int sums_size,
                   void *user_data);

// More functions to manipulate output of MLMCRun
void
mlmc_estimate_error(const mlmc_output* error, double* discretization_error,
                    double *stat_error);

void
mlmc_compute_estimate(const mlmc_output* error,
                      double *estimate,
                      unsigned int estimate_max_size);

#endif   // __MLMC_H__

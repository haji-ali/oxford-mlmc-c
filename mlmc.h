#ifndef __MLMC_H__
#define __MLMC_H__

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
    void (*fn_mlmc_sample_levels)(unsigned int L,
                                   const unsigned long long* M,
                                   double *sums,
                                   unsigned int sums_size,
                                   void *user_data);

    // user supplied function to sample a single level
    void (*fn_mlmc_sample_level)(unsigned int ell,
                                  unsigned long long M,
                                  double *sums,
                                  unsigned int sums_size,
                                  void *user_data);
   
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
MLMCCreateOptionsForCLI(int argc, char ** argv);

// explicit MLMC options creation
mlmc_options* 
MLMCCreateOptions(unsigned int Lmin, unsigned int LMax,
    unsigned int N0, double alpha, double beta, double gamma,
    unsigned int per_sample, 
    void (*fn_mlmc_sample_levels)(unsigned int L,
                                   const unsigned long long* M,
                                   double *sums,
                                   unsigned int sums_size,
                              void *user_data),
    void (*fn_mlmc_sample_level)(unsigned int ell,
                                  unsigned long long M,
                                  double *sums,
                                  unsigned int sums_size,
                                  void *user_data));

void 
MLMCFreeOptions(mlmc_options*);

mlmc_output* 
MLMCRun(double eps, const mlmc_options*); // Options can be NULL for default

void 
MLMCFreeOutput(mlmc_options*);

// This function call a simple single-level sampler
// for every level.
void 
MLMCSampleLevels(const mlmc_options*,
                      unsigned int L,
                      const unsigned long long* M,
                      double *sums,
                      unsigned int sums_size,
                      void *user_data);

// More functions to manipulate output of MLMCRun
void 
MLMCEstimateError(const mlmc_output* error, double* discretization_error,
                       double *stat_error);

void 
MLMCComputeEstimate(const mlmc_output* error,
                         double *estimate,
                         unsigned int estimate_max_size);
#endif   // __MLMC_H__

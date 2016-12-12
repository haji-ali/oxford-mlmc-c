/**
 * @file mlmc.h
 * @brief Header file for Multilevel Monte Carlo driver.
 * @details Defines a set of useful types and functions for typical
 * applications of Multilevellevel Monte Carlo.
 *
 */
#ifndef __MLMC_H__
#define __MLMC_H__

/**
 * @typedef fn_mlmc_sample_levels_t
 * @brief User supplied sampling function pointer for all levels.
 * @details Users must define a functionof this type to enable the 
 * MLMC driver to sample all levels at once. This is generally only 
 * an adventage to the user application if better utilisation of computer
 * hardware can be achieved in this manner.
 *
 * @param L the final fine-grain level index.
 * @param M an L length array of  samples sizes, 
 *          M[l] is the number of sample to take of level l. 
 * @param sums a array of length sum_size that will contain the 
 *             the function cost, and Monte Carlo estimator raw-moments.
 * @param sums_size total size of the sums array (must take into account
 *                   the number of levels, thenumber of moments and 
 *                   the number of elements per moment).
 * @param user_data pointer to users application specific data.
 */
typedef void (*fn_mlmc_sample_levels_t)(unsigned int L,
                                        const unsigned long long* M,
                                        double *sums,
                                        unsigned int sums_size,
                                        void *user_data);
/**
 * @typedef fn_mlmc_sample_level_t
 * @brief User supplied sampling function pointer for single level.
 * @details Users must define a functionof this type to enable the 
 * MLMC driver to sample a single level. 
 *
 * @param ell the level index to sample.
 * @param M an number of samples to generate. 
 * @param sums a array of length sum_size that will contain the 
 *             the function cost, and Monte Carlo estimator raw-moments.
 * @param sums_size total size of the sums array (must take into account
 *                   number of moments and number of elements per moment).
 * @param user_data pointer to users application specific data.
 */
typedef void (*fn_mlmc_sample_level_t)(unsigned int ell,
                                       unsigned long long M,
                                       double *sums,
                                       unsigned int sums_size,
                                       void *user_data);


/**
 * @typedef s_mlmc_options
 * @brief MLMC input options structure, 
 * @details defines the configuration the MLMC driver will execute under.
 */
typedef struct s_mlmc_options{
    /** @var Lmin
     * @brief Minimum number of levels 
     */
    unsigned int Lmin;  
    /** @var Lmax
     * @brief Maximum number of levels 
     */
    unsigned int Lmax; 

    /** @var N0
     * @brief Initial number of samples to compute variance estimates.

     */
    unsigned long long N0;   

    /** @var alpha
     * @brief Weak convergence rate. 
     * @note Must be positive
     */
    double alpha; 
    /** @var beta
     * @brief Variance convergence rate.
     * @note Must be positive
     */
    double beta;  
    /** @var gamma  
     * @brief Work rate. 
     * @note Must be positive
     */
    double gamma;             
    // It will be difficult to support variable length vectors efficiently.
    // I think we should focus on fixed length vectors for now.
    /** @var per_sample
     * @brief Vector length for vector quantities
     * @todo I feel this variable should be re-named, isn't this the dimensionality of
     * the random variable? 
     */
    unsigned int per_sample; 

    /** @var fn_mlmc_sample_levels
     * @brief user supplied function to sample all levels.
     */
    fn_mlmc_sample_levels_t fn_mlmc_sample_levels;
    /** @var fn_mlmc_sample_levels
     * @brief user supplied function to sample a  single level.
     * @note This function will be ignored if fn_mlmc_sample_levels 
     * is set.
     */
    fn_mlmc_sample_level_t fn_mlmc_sample_level;

    /**
     * @var user_data
     * @brief pointer to application specific data/parameters.
     */
    void *user_data;
} mlmc_options;

/**
 * @typedef s_mlmc_output
 * @brief MLMC output structure, 
 * @details containts estimator moments and compuational cost.
 */
typedef struct s_mlmc_output{
    /**
     * @var L
     * @brief finest-grain level index.
     */
    unsigned int L;
    /**
     * @var moments 
     * @brief number of moments generated
     */
    unsigned int moments;
    /**
     * @var per_moment
     * @brief number of elements per moment
     */
    unsigned int per_moment;

    /**
     * @var M
     * @brief array of sample sizes size L
     */
    unsigned long long* M;  
    /**
     * @var sums
     * @brief array containing moment and cost data. 
     * The array is of length L*moments*per_moment
     */
    double* sums;            
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

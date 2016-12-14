/**
 * @file mlmc.h
 * @brief Header file for Multilevel Monte Carlo driver.
 * @details Defines a set of useful types and functions for typical
 * applications of Multilevel Monte Carlo.
 *
 * @author Mike Giles (Mathematical Institute, University of Oxford)
 * @author Abdul-Lateef Haji-Ali (Mathematical Institute, University of Oxford)
 * @author David Warne (School of Mathematical Sciences, Queensland University of Technology)
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
 * @details Users must define a function of this type to enable the
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
    /** @var eps
     * @brief Relative error tolerance
     * @see mlmc_create_options()
     * @see mlmc_create_options_simple()
     */
    double eps;

    /** @var Lmin
     * @brief Minimum number of levels
     * @see mlmc_set_level_limits()
     */
    unsigned int Lmin;
    /** @var Lmax
     * @brief Maximum number of levels.
     * @see mlmc_set_level_limits()
     */
    unsigned int Lmax;

    /** @var N0
     * @brief Initial number of samples to compute variance estimates.
     * @see mlmc_set_initial_samples_num()
     */
    unsigned long long N0;

    /** @var alpha
     * @brief Weak convergence rate.
     * @note Must be positive
     * @see mlmc_set_exponents()
     */
    double alpha;
    /** @var beta
     * @brief Variance convergence rate.
     * @note Must be positive
     * @see mlmc_set_exponents()
     */
    double beta;
    /** @var gamma
     * @brief Work rate.
     * @note Must be positive
     * @see mlmc_set_exponents()
     */
    double gamma;
    // It will be difficult to support variable length vectors efficiently.
    // I think we should focus on fixed length vectors for now.
    /** @var per_sample
     * @brief Vector length for vector quantities
     * @see mlmc_set_sample_size()
     */
    unsigned int per_sample;

    /** @var fn_mlmc_sample_levels
     * @brief user supplied function to sample all levels.
     * @see mlmc_create_options()
     * @see mlmc_create_option_cli()
     */
    fn_mlmc_sample_levels_t fn_mlmc_sample_levels;
    /** @var fn_mlmc_sample_levels
     * @brief user supplied function to sample a  single level.
     * @note This function will be ignored if fn_mlmc_sample_levels
     * is set.
     * @see mlmc_create_options_simple()
     * @see mlmc_create_option_cli()
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
 * @details contains estimator moments and compuational cost.
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

/**
 * Create mlmc_options based on named program arguments.
 * @param argc Number of arguments, usually first parameter of main()
 * @param argv Array of arguments, usually second parameter of main()
 * @return Pointer to mlmc_options structure
 * @todo Full documentation
 */
mlmc_options*
mlmc_create_option_cli(int argc, char ** argv);

/**
 * Create an MLMC option structure to be used to be used with mlmc_run.
 * This function accepts a pointer to a function of type fn_mlmc_sample_levels
 * All arguemtns
 * @param eps Error tolerance
 * @param fn_mlmc_sample_levels Pointer to function that returns samples on
 *  multiple levels.
 * @return Pointer to mlmc_options structure
 */
mlmc_options*
mlmc_create_option(double eps, fn_mlmc_sample_levels_t fn_mlmc_sample_levels);

/**
 * Create an MLMC option structure to be used to be used with mlmc_run.
 * @param eps Error tolerance
 * @param fn_mlmc_sample_level  Pointer to function that returns samples on
 *  a single level.
 * @return Pointer to mlmc_options structure
 */
mlmc_options*
mlmc_create_option_simple(double eps, fn_mlmc_sample_level_t fn_mlmc_sample_level);

/**
 * Sets the lower and upper bounds on the number of levels in the MLMC algorithm
 * @param opt pointer to mlmc_option, created with mlmc_create_options or
 * mlmc_create_options_simple
 * @param Lmin Minimum number of levels
 * @param Lmax Maximum number of levels
 * @note default Lmin is 2 and Lmax is 32
 */
void
mlmc_set_level_limits(mlmc_options* opt, unsigned int Lmin, unsigned int LMax);

/**
 * Sets the initial number of samples that are used to compute variance estimates
 * @param opt pointer to mlmc_option, created with mlmc_create_options or
 * mlmc_create_options_simple
 * @param N0 Initial number of samples
 * @todo Consider making this function accept an array of initial samples number
 * corresponding to each level (0-extrapolating the other levels)
 * @note default is 100
 */
void
mlmc_set_initial_samples_num(mlmc_options* opt, unsigned int N0);

/**
 * Sets the size of each sample that is to be expected when calling the
 * provided fn_mlmc_sample_levels_t or fn_mlmc_sample_levels_t that was given
 * when calling mlmc_create_option or mlmc_create_option_simple, respectively.
 * @param opt pointer to mlmc_option, created with mlmc_create_options or
 * mlmc_create_options_simple
 * @param per_sample Sample size
 * @note default is 1
 */
void
mlmc_set_sample_size(mlmc_options* opt, unsigned int per_sample);

/**
 * Set the exponents that are used in the MLMC algorithm.
 * @param opt pointer to mlmc_option, created with mlmc_create_options or
 * mlmc_create_options_simple
 * @param alpha weak convergence rate. If negative, this exponent will be
 * fit based on computed data.
 * @param beta variance convergence rate. If negative, this exponent will be
 * fit based on computed data.
 * @param gamma work rate. If negative, this exponent will be
 * fit based on computed data.
 * @note default alpha, beta and gamma are all -1.
 */
void
mlmc_set_exponents(mlmc_options* opt,
                   double alpha, double beta, double gamma,
                   int fit_alpha);

/**
 * Frees the memory allocated by mlmc_create_options or
 * mlmc_create_options_simple
 * @param alpha weak convergence rate. If negative, this exponent will be
 * fit based on computed data.
 * @param beta variance convergence rate. If negative, this exponent will be
 * fit based on computed data.
 * @param gamma work rate. If negative, this exponent will be
 * fit based on computed data.
 * @note default alpha, beta and gamma are all -1.
 */
void
mlmc_free_options(mlmc_options* opt);

/**
 * Runs the full MLMC algorithm based on given options.
 * @param opt Pointer to mlmc_options
 * @return Pointer to mlmc_output that contains estimator moments and compuational cost.
 * @see mlmc_create_options
 * @see mlmc_create_options_simple
 * @see mlmc_create_options_cli
 */
mlmc_output*
mlmc_run(const mlmc_options* opt); // Options can be NULL for default


/**
 * Frees the memory allocated by mlmc_run()
 * @param out Pointer to mlmc_output as returned by mlmc_run()
 */
void
mlmc_free_output(mlmc_output* out);

/**
 * Computes the expectation and its error estimate
 * @param out Pointer to mlmc_output as returned by mlmc_run()
 * @param estimate Estimate of expectation
 * @param discretization_error Estimated discretization error
 * @param stat_error Estimated statistical error
 * @param size size of arrays, estimate, discretization_error and stat_error.
 * @return Number of elements that were actually filled in arrays,
 * estimate, discretization_error and stat_error. This would usually be the same
 * as the sample size set by mlmc_set_sample_size()
 */
unsigned int
mlmc_compute_estimate(const mlmc_output* out,
                      double *estimate,
                      double* discretization_error,
                      double *stat_error,
                      unsigned int size);

#endif   // __MLMC_H__

#include <iostream> // std::cout
#include <cmath> // std::log
#include <complex> // std::complex
#include <numeric> // std::accumulate
#include <vector> // std::vector

#include <secdecutil/integrators/cuba.hpp> // secdecutil::cuba::Divonne
#include <secdecutil/series.hpp> // secdecutil::Series
#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation
#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply

%(includes)s

typedef double real_t;
typedef std::complex<real_t> complex_t;
template<typename T> using nested_series_t = %(nested_series_type)s;

template<typename integrand_return_t, typename integrand_t, typename ...other_types>
nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>> compute
(
    std::vector<nested_series_t<integrand_t>> integrands,
    nested_series_t<integrand_return_t> prefactor,
    secdecutil::Integrator<integrand_return_t,other_types...> integrator
)
{
    // add integrands of sectors (together flag)
    const nested_series_t<integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * prefactor;
}

#define COMPUTE(RESULT, NAME, DEFORMATION_PARAMETERS_MAXIMUM, INTEGRATOR) \
\
NAME::nested_series_t<secdecutil::UncorrelatedDeviation<NAME::integrand_return_t>> RESULT = \
compute \
( \
    NAME::make_integrands( \
                                    real_parameters, \
                                    complex_parameters, \
                                    10000, /*number of presamples*/ \
                                    DEFORMATION_PARAMETERS_MAXIMUM \
                                ), \
    NAME::prefactor(real_parameters,complex_parameters), \
    INTEGRATOR \
)

int main()
{
    const std::vector<real_t> real_parameters = {};
    const std::vector<complex_t> complex_parameters = {};

    auto divonne = secdecutil::cuba::Divonne<complex_t>();
    divonne.flags = 2; // verbose output
    divonne.epsrel = 1e-8;
    divonne.epsabs = 1e-8;
    divonne.maxeval = 1e6;
    divonne.border = 1e-8;

%(computes)s

%(parameter_define)s

// TODO - generate correct regulator nesting
//    const secdecutil::Series<real_t> eps = {1,1,{1},false,"eps"};

%(sum_of_integrals)s

%(parameter_undef)s

    std::cout << result << std::endl;
    
    return 0;
}

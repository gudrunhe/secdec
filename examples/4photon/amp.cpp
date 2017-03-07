#include <iostream> // std::cout
#include <cmath> // std::log
#include <complex> // std::complex
#include <numeric> // std::accumulate
#include <vector> // std::vector

#include <secdecutil/integrators/cuba.hpp> // secdecutil::cuba::Vegas
#include <secdecutil/series.hpp> // secdecutil::Series
#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation
#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply

#include "yyyy_bubble/yyyy_bubble.hpp"
#include "yyyy_box6Dim/yyyy_box6Dim.hpp"


/*
 * pySecDec Master Integrals
 */

// one loop bubble
yyyy_bubble::nested_series_t<secdecutil::UncorrelatedDeviation<yyyy_bubble::integrand_return_t>> bubble(yyyy_bubble::real_t uORt)
{
    using namespace yyyy_bubble;

    const std::vector<real_t> real_parameters{uORt};
    const std::vector<complex_t> complex_parameters{};

    // optimize contour
    const std::vector<nested_series_t<yyyy_bubble::integrand_t>> integrands = yyyy_bubble::make_integrands(real_parameters, complex_parameters
        // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
        // optionally set here as additional arguments.
        );

    // add integrands of sectors (together flag)
    const yyyy_bubble::nested_series_t<yyyy_bubble::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e7;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * yyyy_bubble::prefactor(real_parameters, complex_parameters);
}

// one loop box in 6 dimensions
yyyy_box6Dim::nested_series_t<secdecutil::UncorrelatedDeviation<yyyy_box6Dim::integrand_return_t>> box6Dim(yyyy_box6Dim::real_t t, yyyy_box6Dim::real_t u)
{
    using namespace yyyy_box6Dim;

    const std::vector<real_t> real_parameters{t,u};
    const std::vector<complex_t> complex_parameters{};

    // optimize contour
    const std::vector<nested_series_t<yyyy_box6Dim::integrand_t>> integrands = yyyy_box6Dim::make_integrands(real_parameters, complex_parameters
        // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
        // optionally set here as additional arguments.
        );

    // add integrands of sectors (together flag)
    const yyyy_box6Dim::nested_series_t<yyyy_box6Dim::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * yyyy_box6Dim::prefactor(real_parameters, complex_parameters);
}

/*
 * numerical amplitude using pySecDec Master Integrals
 */
secdecutil::Series<secdecutil::UncorrelatedDeviation<std::complex<double>>> yyyy_numerical(double s, double t, double u)
{
    return -8.*( 1. + (t*t + u*u)/s * box6Dim(t,u) +  (t-u)/s*( bubble(u)-bubble(t) ) );
}

/*
 * analytic result
 */
double yyyy_analytic(double s, double t, double u)
{
    constexpr double M_pi = 3.14159265358979323846;
    const double L = std::log(t/u);
    return -8.*( 1. + (t-u)/s*L  + (t*t + u*u)/(2.*s*s)*( M_pi*M_pi + L*L ) );
}

int main()
{
    /*
     * phase-space point
     */
    const double s = .9, t = -.1;
    const double u = -s -t;

    /*
     * compute the amplitude
     */
    const auto numerical_result = yyyy_numerical(s,t,u);
    const auto analytic_result = yyyy_analytic(s,t,u);

    /*
     * print the amplitude
     */
    std::cout << "------------" << std::endl << std::endl;
    std::cout << "amplitude M++--" << std::endl << std::endl;

    std::cout << "numerical result" << std::endl;
    std::cout << numerical_result << std::endl << std::endl;

    std::cout << "analytic result" << std::endl;
    std::cout << analytic_result << " + O(eps)" << std::endl << std::endl;

    return 0;
}

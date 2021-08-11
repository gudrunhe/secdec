#include <iostream> // std::cout
#include <cmath> // std::log
#include <numeric> // std::accumulate
#include <vector> // std::vector

#include <secdecutil/integrators/cuba.hpp> // secdecutil::cuba::Vegas
#include <secdecutil/series.hpp> // secdecutil::Series
#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation
#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply

#include "box2L_invprop/box2L_invprop.hpp"
#include "box2L_contracted_tensor/box2L_contracted_tensor.hpp"

/*
 * the integral from pySecDec using two different input forms
 */

// reference result from pySecDec paper Table 2 arXiv:1703.09692v2
secdecutil::Series<secdecutil::UncorrelatedDeviation<std::complex<double>>> box2L_reference
{
    -4, // order min
     0, // order max
    {
        // {<value>, <uncertainty>}
        {-0.2916, 0.0022}, // eps^-4
        { 0.7410, 0.0076}, // eps^-3
        {-0.3056, 0.0095}, // eps^-2
        {-2.2966, 0.0313}, // eps^-1
        { 1.1460, 0.0504}  // eps^0
    },
    true, // series is truncated
    "eps" // expansion parameter
};

// 'box2L_invprop' - input numerator as inverse propagator
box2L_invprop::nested_series_t<secdecutil::UncorrelatedDeviation<box2L_invprop::integrand_return_t>> compute_box2L_invprop(box2L_invprop::real_t s, box2L_invprop::real_t t)
{
    using namespace box2L_invprop;

    const std::vector<real_t> real_parameters{s,t};
    const std::vector<complex_t> complex_parameters{};

    // Set up Integrator
    secdecutil::integrators::Qmc<
                                    integrand_return_t,
                                    maximal_number_of_integration_variables,
                                    integrators::transforms::Korobov<3>::type,
                                    user_integrand_t
                                > integrator;
    integrator.verbosity = 1;
    integrator.epsrel = 5e-3;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;
    
    // Construct the amplitudes
    std::vector<nested_series_t<sum_t>> unwrapped_amplitudes =
        make_amplitudes(real_parameters, complex_parameters, "box2L_invprop/box2L_invprop_data", integrator);

    // Pack amplitudes into handler
    handler_t<amplitudes_t> amplitudes
    (
        unwrapped_amplitudes,
        integrator.epsrel, integrator.epsabs
        // further optional arguments: maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
    );

        // compute the amplitudes
    const std::vector<nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>> result = amplitudes.evaluate();
    
    return result.at(0);
}

// 'box2L_contracted_tensor' - input numerator as contracted tensor
box2L_contracted_tensor::nested_series_t<secdecutil::UncorrelatedDeviation<box2L_contracted_tensor::integrand_return_t>> compute_box2L_contracted_tensor(box2L_invprop::real_t s, box2L_invprop::real_t t)
{
    using namespace box2L_contracted_tensor;

    const std::vector<real_t> real_parameters{s,t};
    const std::vector<complex_t> complex_parameters{};

    // Set up Integrator
    secdecutil::integrators::Qmc<
                                    integrand_return_t,
                                    maximal_number_of_integration_variables,
                                    integrators::transforms::Korobov<3>::type,
                                    user_integrand_t
                                > integrator;
    integrator.verbosity = 1;
    integrator.epsrel = 5e-3;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;

    // Construct the amplitudes
    std::vector<nested_series_t<sum_t>> unwrapped_amplitudes =
        make_amplitudes(real_parameters, complex_parameters, "box2L_contracted_tensor/box2L_contracted_tensor_data", integrator);

    // Pack amplitudes into handler
    handler_t<amplitudes_t> amplitudes
    (
        unwrapped_amplitudes,
        integrator.epsrel, integrator.epsabs
        // further optional arguments: maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
    );

        // compute the amplitudes
    const std::vector<nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>> result = amplitudes.evaluate();
    
    return result.at(0);
}

int main()
{
    /*
     * phase-space point
     */
    const double s = -3.0, t = -2.0;

    /*
     * compute the two integrals
     */
    const auto box2L_invprop = compute_box2L_invprop(s,t);
    const auto box2L_contracted_tensor = compute_box2L_contracted_tensor(s,t);

    /*
     * print the result
     */
    std::cout << "------------" << std::endl << std::endl;
    std::cout << "box2L" << std::endl << std::endl;

    std::cout << "result from pySecDec paper" << std::endl;
    std::cout << box2L_reference << std::endl << std::endl;

    std::cout << "result inverse propagator" << std::endl;
    std::cout << box2L_invprop << std::endl << std::endl;

    std::cout << "result contracted tensor" << std::endl;
    std::cout << box2L_contracted_tensor << std::endl << std::endl;

    std::cout << "<inverse propagator> - <contracted tensor>" << std::endl;
    std::cout << box2L_invprop - box2L_contracted_tensor << std::endl << std::endl;

    return 0;
}

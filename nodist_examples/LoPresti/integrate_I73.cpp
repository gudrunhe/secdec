#include <cmath> // std::sqrt
#include <iostream> // std::cout
#include <vector> // std::vector

#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation

#include "I73.hpp"

int main()
{
    const double v1 = -3, v2 = -3, v3 = -1, v4 = -1, v5 = -1;
    const double sqrtDelta = std::sqrt(v1*v1 * (v2-v5)*(v2-v5) + (v2*v3 + v4*(v5-v3))*(v2*v3 + v4*(v5-v3)) + 2*v1 * (-v2*v2*v3 + v4*(v3-v5)*v5 + v2*(v3*v4 + (v3 + v4) * v5)));

    // User Specified Phase-space point
    const std::vector<I73::real_t> real_parameters = { -3,-3,-1,-1,-1,sqrtDelta };
    const std::vector<I73::complex_t> complex_parameters = {  };

    // Construct the amplitude. Further options for the individual integrals can be set in the
    // corresponding "<integral_name>_weighted_integral.cpp" file in the "src/" directory
    I73::nested_series_t<I73::sum_t> unwrapped_amplitude =
        I73::make_amplitude(real_parameters, complex_parameters);

    // pack amplitude into handler
    I73::handler_t<I73::nested_series_t> amplitude
    (
        unwrapped_amplitude,
        1e-6, 1e-14 // further optional arguments: epsrel, epsabs, maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
    );
    amplitude.verbose = true;

    // The optional further arguments of the handler are set for all orders.
    // To specify different settings for a particular order, type e.g.: amplitude.expression.at(0).epsrel = 1e-5;

    // optionally set wall clock limits (in seconds)
    // Note: Only the wall clock time spent in "amplitude.evaluate()" is considered for these limits.
    amplitude.soft_wall_clock_limit = 4 * 60 * 60;
    amplitude.hard_wall_clock_limit = 5 * 60 * 60;

    // optionally compute multiple integrals concurrently
    // Note: The integrals themselves may also be computed in parallel irresprective of this option.
    amplitude.number_of_threads = 16;

    // The cuda driver does not automatically remove unnecessary functions from the device memory
    // such that the device may run out of memry after some time. This option controls how many
    // after how many integrals "cudaDeviceReset()" is called to clear the memory. With the default
    // "0", "cudaDeviceReset()" is never called. This option is ignored if compiled without cuda.
    amplitude.reset_cuda_after = 2000;

    // compute the amplitude
    const I73::nested_series_t<secdecutil::UncorrelatedDeviation<I73::integrand_return_t>> result = amplitude.evaluate();

    // print the result
    std::cout << "amplitude = " << result << std::endl;
}

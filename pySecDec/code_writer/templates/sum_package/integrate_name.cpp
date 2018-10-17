#include <iostream> // std::cout
#include <vector> // std::vector

#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation

#include "%(name)s.hpp"

int main()
{
    // User Specified Phase-space point
    const std::vector<%(name)s::real_t> real_parameters = { /* EDIT: insert real parameter values here */ };
    const std::vector<%(name)s::complex_t> complex_parameters = { /* EDIT: insert complex parameter values here */ };

    // Construct the amplitude. Further options for the individual integrals can be set in the
    // corresponding "<integral_name>_weighted_integral.cpp" file in the "src/" directory
    %(name)s::nested_series_t<%(name)s::sum_t> unwrapped_amplitude =
        %(name)s::make_amplitude(real_parameters, complex_parameters);

    // pack amplitude into handler
    %(name)s::handler_t<%(name)s::nested_series_t> amplitude
    (
        unwrapped_amplitude
        // further optional arguments: epsrel, epsabs, maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
    );
    amplitude.verbose = true;

    // The optional further arguments of the handler are set for all orders.
    // To specify different settings for a particular order, type e.g.: amplitude.expression.at(0).epsrel = 1e-5;

    // optionally set wall clock limits (in seconds)
    // Note: Only the wall clock time spent in "amplitude.evaluate()" is considered for these limits.
    // amplitude.soft_wall_clock_limit = 60 *  8;
    // amplitude.hard_wall_clock_limit = 60 * 10;

    // optionally compute multiple integrals concurrently
    // Note: The integrals themselves may also be computed in parallel irresprective of this option.
    // amplitude.number_of_threads = 12;

    // The cuda driver does not automatically remove unneccessary functions from the device memory
    // such that the device may run out of memry after some time. This option controls how many
    // after how many integrals "cudaDeviceReset()" is called to clear the memory. With the default
    // "0", "cudaDeviceReset()" is never called. This option is ignored if compiled without cuda.
    // amplitude.reset_cuda_after = 2000;

    // compute the amplitude
    const %(name)s::nested_series_t<secdecutil::UncorrelatedDeviation<%(name)s::integrand_return_t>> result = amplitude.evaluate();

    // print the result
    std::cout << "amplitude = " << result << std::endl;
}

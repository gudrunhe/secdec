#include "%(name)s.hpp"

#include <vector>
#include <memory> // std::shared_ptr, std::make_shared
#include <string>
#include <sstream>

#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation

#define INTEGRAL_NAME %(name)s

#include <secdecutil/pylink.hpp> // The python-C binding is general and therefore contained in the util
#include <secdecutil/pylink_amplitude.hpp>

using namespace INTEGRAL_NAME;

template<typename T> using amplitudes_t = std::vector<nested_series_t<T>>;

extern "C"
{
    int compute_integral
    (
        std::string * integral_without_prefactor_strptr, std::string * prefactor_strptr, std::string * integral_with_prefactor_strptr, // output
        const secdecutil::Integrator<integrand_return_t,real_t> * integrator, // pointer to the integrator
        const double real_parameters_input[], // real parameters
        const double complex_parameters_input[], // complex parameters serialized as real(x0), imag(x0), real(x1), imag(x1), ...
        const bool together, // integrate sectors together
        const unsigned number_of_presamples,
        const real_t deformation_parameters_maximum,
        const real_t deformation_parameters_minimum,
        const real_t deformation_parameters_decrease_factor,
        const real_t epsrel,
        const real_t epsabs,
        const unsigned long long int maxeval,
        const unsigned long long int mineval,
        const real_t maxincreasefac,
        const real_t min_epsrel,
        const real_t min_epsabs,
        const real_t max_epsrel,
        const real_t max_epsabs,
        const real_t min_decrease_factor,
        const real_t decrease_to_percentage, // of remaining time
        const real_t soft_wall_clock_limit,
        const real_t hard_wall_clock_limit,
        const size_t number_of_threads,
        const size_t reset_cuda_after,
        const int verbose
    )
    {
        size_t i;
        std::stringstream sstream;

        // fix output formatting
        sstream.precision(std::numeric_limits<real_t>::max_digits10); // force enough digits to ensure unique recreation
        sstream << std::scientific; // stringify floats as #.#e#
        
        // read real parameters
        std::vector<real_t> real_parameters(number_of_real_parameters);
        for (i=0 ; i<number_of_real_parameters ; ++i)
            real_parameters[i] = real_parameters_input[i];

        // read complex parameters
        std::vector<complex_t> complex_parameters(number_of_complex_parameters);
        for (i=0 ; i<number_of_complex_parameters ; ++i)
            complex_parameters[i] = complex_t(complex_parameters_input[2*i],complex_parameters_input[2*i + 1]);

        // Construct the amplitudes
        std::vector<nested_series_t<sum_t>> unwrapped_amplitudes = make_amplitudes(real_parameters, complex_parameters, integrator);

        // pack amplitude into handler
        handler_t<amplitudes_t> amplitudes
        (
            unwrapped_amplitudes,
            epsrel, epsabs, maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
        );
        amplitudes.min_decrease_factor = min_decrease_factor;
        amplitudes.decrease_to_percentage = decrease_to_percentage;
        amplitudes.soft_wall_clock_limit = soft_wall_clock_limit;
        amplitudes.hard_wall_clock_limit = hard_wall_clock_limit;
        amplitudes.number_of_threads = number_of_threads;
        amplitudes.reset_cuda_after = reset_cuda_after;
        amplitudes.verbose = verbose;

        // compute the amplitude
        std::vector<nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>> result;
        try {
            result = amplitudes.evaluate();
        } catch (std::exception& e){
            std::cout << "Encountered an exception of type '" << typeid(e).name() << "'" << std::endl;
            std::cout << "  what():  " << e.what() << std::endl;
            return -1;
        }
        
        // populate output strings:
        //   - integral without prefactor
        sstream.str("");
        for(const auto amp : result)
            sstream << amp;
        *integral_without_prefactor_strptr = sstream.str();
        
        //   - prefactor
        sstream.str("");
        sstream << "1";
        *prefactor_strptr = sstream.str();
        
        //   - full result (prefactor*integral)
        sstream.str("");
        for(const auto amp : result)
            sstream << amp;
        *integral_with_prefactor_strptr = sstream.str();
        
        return 0;
    }

}

#ifndef SecDecUtil_pylink_amplitude_hpp_included
#define SecDecUtil_pylink_amplitude_hpp_included

#include <iostream>
#include <limits> // std::numeric_limits
#include <memory> // std::unique_ptr
#include <numeric> // std::accumulate
#include <string>
#include <sstream>
#include <vector>

#include <secdecutil/deep_apply.hpp> // deep_apply
#include <secdecutil/series.hpp> // Series
#include <secdecutil/uncertainties.hpp> // UncorrelatedDeviation

/*
 * First perform replacements, then put item into quotes.
 * With only a single macro, the preprocessor does not
 * consider replacements in "item".
 */
#define EXPAND_STRINGIFY(item) STRINGIFY(item)
#define STRINGIFY(item) #item

extern "C"
{
    /*
     * get integral info
     */
    void get_integral_info(std::string * str_ptr)
    {
        std::stringstream sstream;

        // fix output formatting
        sstream.precision(std::numeric_limits<real_t>::max_digits10); // force enough digits to ensure unique recreation
        sstream << std::scientific; // stringify floats as #.#e#

        sstream << "name = " << EXPAND_STRINGIFY(INTEGRAL_NAME) << std::endl;

        sstream << "number_of_regulators = " << number_of_regulators << std::endl;
        sstream << "names_of_regulators =";
        for ( const auto& name : names_of_regulators )
            sstream << " " << name;
        sstream << std::endl;

        sstream << "number_of_real_parameters = " << number_of_real_parameters << std::endl;
        sstream << "names_of_real_parameters =";
        for ( const auto& name : names_of_real_parameters )
            sstream << " " << name;
        sstream << std::endl;

        sstream << "number_of_complex_parameters = " << number_of_complex_parameters << std::endl;
        sstream << "names_of_complex_parameters =";
        for ( const auto& name : names_of_complex_parameters )
            sstream << " " << name;
        // sstream << std::endl;

        *str_ptr = sstream.str();

    }

    /*
     * function to compute the integral
     */
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
        const real_t wall_clock_limit,
        const size_t number_of_threads,
        const size_t reset_cuda_after,
        const bool verbose,
        const int errormode_enum,
        const char *lib_path
    )
    {
        int i;
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
        if(verbose) std::cerr << "Generating amplitudes (optimising contour if required)" << std::endl;
        std::vector<nested_series_t<sum_t>> unwrapped_amplitudes = make_amplitudes(
                                                                                          real_parameters,
                                                                                          complex_parameters,
                                                                                          std::string(lib_path),
                                                                                          integrator
                                                                                          #if integral_contour_deformation
                                                                                              , number_of_presamples,
                                                                                              deformation_parameters_maximum,
                                                                                              deformation_parameters_minimum,
                                                                                              deformation_parameters_decrease_factor
                                                                                          #endif
                                                                                          );

        // pack amplitude into handler
        if(verbose) std::cerr << "Packing amplitudes into handler" << std::endl;
        handler_t<amplitudes_t> amplitudes
        (
            unwrapped_amplitudes,
            epsrel, epsabs, maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
        );
        amplitudes.min_decrease_factor = min_decrease_factor;
        amplitudes.decrease_to_percentage = decrease_to_percentage;
        amplitudes.wall_clock_limit = wall_clock_limit;
        amplitudes.number_of_threads = number_of_threads;
        amplitudes.reset_cuda_after = reset_cuda_after;
        amplitudes.verbose = verbose;
        amplitudes.errormode = static_cast<handler_t<amplitudes_t>::ErrorMode>(errormode_enum);

        // compute the amplitude
        if(verbose) std::cerr << "Integrating" << std::endl;
        std::vector<nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>> result;
        try {
            result = amplitudes.evaluate();
        } catch (std::exception& e){
            std::cerr << "Encountered an exception of type '" << typeid(e).name() << "'" << std::endl;
            std::cerr << "  what():  " << e.what() << std::endl;
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
            sstream << amp << "\n";
        *integral_with_prefactor_strptr = sstream.str();
        
        return 0;
    }

    /*
     * function to compute the integral using cuda
     */
    #ifdef SECDEC_WITH_CUDA
        int cuda_compute_integral
        (
            std::string * integral_without_prefactor_strptr, std::string * prefactor_strptr, std::string * integral_with_prefactor_strptr, // output
            const secdecutil::Integrator<integrand_return_t,real_t,cuda_together_integrand_t> * together_integrator, // pointer to the integrator for together=true (not used)
            const secdecutil::Integrator<integrand_return_t,real_t,cuda_integrand_t> * separate_integrator, // pointer to the integrator for together=false
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
            const real_t wall_clock_limit,
            const size_t number_of_threads,
            const size_t reset_cuda_after,
            const bool verbose,
            const int errormode_enum,
            const char *lib_path
        )
        {
            int i;
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
            if(verbose) std::cerr << "Generating amplitudes (optimising contour if required)" << std::endl;
            std::vector<nested_series_t<sum_t>> unwrapped_amplitudes = make_amplitudes(
                                                                                              real_parameters,
                                                                                              complex_parameters,
                                                                                              std::string(lib_path),
                                                                                              separate_integrator
                                                                                              #if integral_contour_deformation
                                                                                                  , number_of_presamples,
                                                                                                  deformation_parameters_maximum,
                                                                                                  deformation_parameters_minimum,
                                                                                                  deformation_parameters_decrease_factor
                                                                                              #endif
                                                                                              );

            // pack amplitude into handler
            if(verbose) std::cerr << "Packing amplitudes into handler" << std::endl;
            handler_t<amplitudes_t> amplitudes
            (
                unwrapped_amplitudes,
                epsrel, epsabs, maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
            );
            amplitudes.min_decrease_factor = min_decrease_factor;
            amplitudes.decrease_to_percentage = decrease_to_percentage;
            amplitudes.wall_clock_limit = wall_clock_limit;
            amplitudes.number_of_threads = number_of_threads;
            amplitudes.reset_cuda_after = reset_cuda_after;
            amplitudes.verbose = verbose;

            // compute the amplitude
            if(verbose) std::cerr << "Integrating" << std::endl;
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
                sstream << amp << "\n";
            *integral_with_prefactor_strptr = sstream.str();
            
            return 0;
        }
    #endif

    #undef EXPAND_STRINGIFY
    #undef STRINGIFY

}

#endif

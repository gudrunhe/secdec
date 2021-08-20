#ifndef SecDecUtil_pylink_integral_hpp_included
#define SecDecUtil_pylink_integral_hpp_included

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

        sstream << "number_of_sectors = " << number_of_sectors << std::endl;

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
        sstream << std::endl;

        sstream << "lowest_orders =";
        for ( const auto& lowest_order : lowest_orders )
            sstream << " " << lowest_order;
        sstream << std::endl;

        sstream << "highest_orders =";
        for ( const auto& highest_order : highest_orders )
            sstream << " " << highest_order;
        sstream << std::endl;

        sstream << "lowest_prefactor_orders =";
        for ( const auto& highest_order : lowest_prefactor_orders )
            sstream << " " << highest_order;
        sstream << std::endl;

        sstream << "highest_prefactor_orders =";
        for ( const auto& highest_order : highest_prefactor_orders )
            sstream << " " << highest_order;
        sstream << std::endl;

        sstream << "requested_orders =";
        for ( const auto& requested_order : requested_orders )
            sstream << " " << requested_order;
        sstream << std::endl;

        sstream << "pole_structures =";
        for ( const auto& polestruct : pole_structures )
        {
            for ( const auto& variable_power : polestruct )
                sstream << " " << variable_power;
            sstream << " , ";
        }

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
     
        // following parameters unused (required for compatibility with pylink_amplitude)
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

        // optimize the deformation (if any)
        if(verbose) std::cerr << "Generating integrands (optimising contour if required)" << std::endl;
        const std::vector<nested_series_t<secdec_integrand_t>> sector_integrands =
        make_integrands
        (
            real_parameters, complex_parameters
            #if integral_contour_deformation
                ,number_of_presamples,
                deformation_parameters_maximum,
                deformation_parameters_minimum,
                deformation_parameters_decrease_factor
            #endif
        );

        std::unique_ptr<nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>> result_all;
        try{
            if (together) {
                // add integrands of sectors (together flag)
                if(verbose) std::cerr << "Summing integrands" << std::endl;
                const nested_series_t<secdec_integrand_t> all_sectors = std::accumulate( ++sector_integrands.begin(), sector_integrands.end(), *sector_integrands.begin() );

                // perform the integration
                if(verbose) std::cerr << "Integrating" << std::endl;
                result_all.reset
                (
                    new nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>
                    (
                        secdecutil::deep_apply( all_sectors, integrator->integrate )
                    )
                );
            } else {
                // perform the integration
                if(verbose) std::cerr << "Integrating" << std::endl;
                const std::vector<nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>> integrated_sectors = secdecutil::deep_apply( sector_integrands, integrator->integrate );

                // add integrated sectors
                if(verbose) std::cerr << "Summing integrals" << std::endl;
                result_all.reset
                (
                    new nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>
                    (
                        std::accumulate( ++integrated_sectors.begin(), integrated_sectors.end(), *integrated_sectors.begin() )
                    )
                );
            }
        } catch (std::exception& e){
            std::cerr << "Encountered an exception of type '" << typeid(e).name() << "'" << std::endl;
            std::cerr << "  what():  " << e.what() << std::endl;
            return -1;
        }

        // populate output strings:
        //   - integral without prefactor
        sstream.str("");
        sstream << *result_all;
        *integral_without_prefactor_strptr = sstream.str();

        //   - prefactor
        const nested_series_t<integrand_return_t> evaluated_prefactor = prefactor(real_parameters, complex_parameters);
        sstream.str("");
        sstream << evaluated_prefactor;
        *prefactor_strptr = sstream.str();

        //   - full result (prefactor*integral)
        sstream.str("");
        sstream << evaluated_prefactor * (*result_all);
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
            const secdecutil::Integrator<integrand_return_t,real_t,cuda_together_integrand_t> * together_integrator, // pointer to the integrator if together=true
            const secdecutil::Integrator<integrand_return_t,real_t,cuda_integrand_t> * separate_integrator, // pointer to the integrator if together=false
            const double real_parameters_input[], // real parameters
            const double complex_parameters_input[], // complex parameters serialized as real(x0), imag(x0), real(x1), imag(x1), ...
            const bool together, // integrate sectors together
            const unsigned number_of_presamples,
            const real_t deformation_parameters_maximum,
            const real_t deformation_parameters_minimum,
            const real_t deformation_parameters_decrease_factor,
         
            // following parameters unused (required for compatibility with pylink_amplitude)
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

            // optimize the deformation (if any)
            if(verbose) std::cerr << "Generating integrands (optimising contour if required)" << std::endl;
            const std::vector<nested_series_t<cuda_integrand_t>> sector_integrands =
            make_cuda_integrands
            (
                real_parameters, complex_parameters
                #if integral_contour_deformation
                    ,number_of_presamples,
                    deformation_parameters_maximum,
                    deformation_parameters_minimum,
                    deformation_parameters_decrease_factor
                #endif
            );

            std::unique_ptr<nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>> result_all;
            try{
                if (together) {
                    // add integrands of sectors (together flag)
                    if(verbose) std::cerr << "Summing integrands" << std::endl;
                    const nested_series_t<cuda_together_integrand_t> all_sectors =
                        std::accumulate( ++sector_integrands.begin(), sector_integrands.end(), cuda_together_integrand_t()+*sector_integrands.begin() );

                    // perform the integration
                    if(verbose) std::cerr << "Integrating" << std::endl;
                    result_all.reset
                    (
                        new nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>
                        (
                            secdecutil::deep_apply( all_sectors, together_integrator->integrate )
                        )
                    );
                } else {
                    // perform the integration
                    if(verbose) std::cerr << "Integrating" << std::endl;
                    const std::vector<nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>> integrated_sectors = secdecutil::deep_apply( sector_integrands, separate_integrator->integrate );

                    // add integrated sectors
                    if(verbose) std::cerr << "Summing integrals" << std::endl;
                    result_all.reset
                    (
                        new nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>
                        (
                            std::accumulate( ++integrated_sectors.begin(), integrated_sectors.end(), *integrated_sectors.begin() )
                        )
                    );
                }
            } catch (std::exception& e){
                std::cerr << "Encountered an exception of type '" << typeid(e).name() << "'" << std::endl;
                std::cerr << "  what():  " << e.what() << std::endl;
                return -1;
            }

            // populate output strings:
            //   - integral without prefactor
            sstream.str("");
            sstream << *result_all;
            *integral_without_prefactor_strptr = sstream.str();

            //   - prefactor
            const nested_series_t<integrand_return_t> evaluated_prefactor = prefactor(real_parameters, complex_parameters);
            sstream.str("");
            sstream << evaluated_prefactor;
            *prefactor_strptr = sstream.str();

            //   - full result (prefactor*integral)
            sstream.str("");
            sstream << evaluated_prefactor * (*result_all);
            *integral_with_prefactor_strptr = sstream.str();

            return 0;
        }
    #endif

    #undef EXPAND_STRINGIFY
    #undef STRINGIFY

}

#endif

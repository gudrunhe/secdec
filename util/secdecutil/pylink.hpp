#ifndef SecDecUtil_pylink_hpp_included
#define SecDecUtil_pylink_hpp_included

#include <iostream>
#include <limits> // std::numeric_limits
#include <memory> // std::unique_ptr
#include <numeric> // std::accumulate
#include <string>
#include <sstream>
#include <vector>

#include <secdecutil/deep_apply.hpp> // deep_apply
#include <secdecutil/integrators/cquad.hpp> // CQuad
#include <secdecutil/integrators/qmc.hpp> // Qmc
#include <secdecutil/integrators/cuba.hpp> // Vegas, Suave, Divonne, Cuhre
#include <secdecutil/series.hpp> // Series
#include <secdecutil/uncertainties.hpp> // UncorrelatedDeviation


extern "C"
{

    #if integral_has_complex_parameters || integral_contour_deformation || integral_enforce_complex_return_type
        #define SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX(NOARGS) do {integrator->together = real_complex_together;} while (false)
    #else
        #define SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX(NOARGS)
    #endif

    /*
     * First perform replacements, then put item into quotes.
     * With only a single macro, the preprocessor does not
     * consider replacements in "item".
     */
    #define EXPAND_STRINGIFY(item) STRINGIFY(item)
    #define STRINGIFY(item) #item

    using namespace INTEGRAL_NAME;
    using secdec_integrand_t = INTEGRAL_NAME::integrand_t; // avoid name conflict with cuba

    /*
     * string (de)allocation
     */
    const char * string2charptr(std::string * str)
    {
        return str->c_str();
    }
    std::string * allocate_string()
    {
        return new std::string;
    }
    void free_string(std::string * strptr)
    {
        delete strptr;
    }

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
     * integrator (de)allocation
     */
    secdecutil::Integrator<integrand_return_t,real_t> *
    allocate_MultiIntegrator(
                                 secdecutil::Integrator<integrand_return_t,real_t>* low_dim_integrator,
                                 secdecutil::Integrator<integrand_return_t,real_t>* high_dim_integrator,
                                 int critical_dim
                            )
    {
        auto integrator = new secdecutil::MultiIntegrator<integrand_return_t,real_t>
            (*low_dim_integrator,*high_dim_integrator,critical_dim);
        return integrator;
    }

    secdecutil::Integrator<integrand_return_t,real_t> *
    allocate_gsl_cquad(
                            double epsrel,
                            double epsabs,
                            unsigned int n,
                            bool verbose,
                            double zero_border
                       )
    {
        auto integrator = new secdecutil::gsl::CQuad<integrand_return_t>
            (epsrel,epsabs,n,verbose,zero_border);
        return integrator;
    }

    secdecutil::Integrator<integrand_return_t,real_t> *
    allocate_cuba_Vegas(
                            double epsrel,
                            double epsabs,
                            int flags,
                            int seed,
                            long long int mineval,
                            long long int maxeval,
                            double zero_border,
                            long long int nstart,
                            long long int nincrease,
                            long long int nbatch,
                            bool real_complex_together
                       )
    {
        auto integrator = new secdecutil::cuba::Vegas<integrand_return_t>
            (epsrel,epsabs,flags,seed,mineval,maxeval,zero_border,nstart,nincrease,nbatch);
        SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX();
        return integrator;
    }
    secdecutil::Integrator<integrand_return_t,real_t> *
    allocate_cuba_Suave(
                            double epsrel,
                            double epsabs,
                            int flags,
                            int seed,
                            long long int mineval,
                            long long int maxeval,
                            double zero_border,
                            long long int nnew,
                            long long int nmin,
                            double flatness,
                            bool real_complex_together
                       )
    {
        auto integrator = new secdecutil::cuba::Suave<integrand_return_t>
            (epsrel,epsabs,flags,seed,mineval,maxeval,zero_border,nnew,nmin,flatness);
        SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX();
        return integrator;
    }
    secdecutil::Integrator<integrand_return_t,real_t> *
    allocate_cuba_Divonne(
                            double epsrel,
                            double epsabs,
                            int flags,
                            int seed,
                            long long int mineval,
                            long long int maxeval,
                            double zero_border,
                            int key1,
                            int key2,
                            int key3,
                            int maxpass,
                            double border,
                            double maxchisq,
                            double mindeviation,
                            bool real_complex_together
                         )
    {
        auto integrator = new secdecutil::cuba::Divonne<integrand_return_t>
            (epsrel,epsabs,flags,seed,mineval,maxeval,zero_border,key1,key2,key3,maxpass,
             border,maxchisq,mindeviation);
        SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX();
        return integrator;
    }
    secdecutil::Integrator<integrand_return_t,real_t> *
    allocate_cuba_Cuhre(
                            double epsrel,
                            double epsabs,
                            int flags,
                            long long int mineval,
                            long long int maxeval,
                            double zero_border,
                            int key,
                            bool real_complex_together
                       )
    {
        auto integrator = new secdecutil::cuba::Cuhre<integrand_return_t>
            (epsrel,epsabs,flags,mineval,maxeval,zero_border,key);
        SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX();
        return integrator;
    }
    secdecutil::Integrator<integrand_return_t,real_t> *
    allocate_integrators_Qmc(
                            unsigned long long int minN,
                            unsigned long long int m,
                            unsigned long long int blockSize,
                            long long int seed
                        )
    {
        auto integrator = new secdecutil::integrators::Qmc<integrand_return_t>;
        // If an argument is set to 0 then use the default of the Qmc library
        if ( minN != 0 )
            integrator->integrator.minN = minN;
        if ( m != 0 )
            integrator->integrator.m = m;
        if ( blockSize != 0 )
            integrator->integrator.blockSize = blockSize;
        if ( seed != 0 )
            integrator->integrator.randomGenerator.seed(seed);
        return integrator;
    }
    void free_integrator (secdecutil::Integrator<integrand_return_t,real_t> * integrator)
    {
        delete integrator;
    }

    /*
     * function to compute the integral
     */
    void compute_integral
    (
        std::string * integral_without_prefactor_strptr, std::string * prefactor_strptr, std::string * integral_with_prefactor_strptr, // output
        const secdecutil::Integrator<integrand_return_t,real_t> * integrator, // pointer to the integrator
        const double real_parameters_input[], // real parameters
        const double complex_parameters_input[], // complex parameters serialized as real(x0), imag(x0), real(x1), imag(x1), ...
        const bool together, // integrate sectors together
        const unsigned number_of_presamples,
        const real_t deformation_parameters_maximum,
        const real_t deformation_parameters_minimum,
        const real_t deformation_parameters_decrease_factor
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
            complex_parameters[i] = {complex_parameters_input[2*i],complex_parameters_input[2*i + 1]};

        // optimize the deformation (if any)
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
        if (together) {
            // add integrands of sectors (together flag)
            const nested_series_t<secdec_integrand_t> all_sectors = std::accumulate( ++sector_integrands.begin(), sector_integrands.end(), *sector_integrands.begin() );

            // perform the integration
            result_all.reset
            (
                new nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>
                (
                    secdecutil::deep_apply( all_sectors, integrator->integrate )
                )
            );
        } else {
            // perform the integration
            const std::vector<nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>> integrated_sectors = secdecutil::deep_apply( sector_integrands, integrator->integrate );

            // add integrated sectors
            result_all.reset
            (
                new nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>>
                (
                    std::accumulate( ++integrated_sectors.begin(), integrated_sectors.end(), *integrated_sectors.begin() )
                )
            );
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
    }

    #undef SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX
    #undef EXPAND_STRINGIFY
    #undef STRINGIFY

}

#endif

#ifndef %(name)s_hpp_included
#define %(name)s_hpp_included

#ifdef SECDEC_WITH_CUDA
    #include <thrust/complex.h>
#else
    #include <complex>
#endif
#include <string>
#include <vector>
#include <string>

#include <secdecutil/integrand_container.hpp>
#include <secdecutil/series.hpp>
#include <secdecutil/uncertainties.hpp>
#include <secdecutil/amplitude.hpp> // secdecutil::amplitude::Integral, secdecutil::amplitude::WeightedIntegral, secdecutil::amplitude::WeightedIntegralHandler

namespace %(name)s
{
    // whether or not to use contour deformation
    #define %(name)s_contour_deformation %(contour_deformation)i

    // some information about the integral
    // --{
    const unsigned long long number_of_integrals = %(number_of_integrals)i;
    const unsigned int number_of_amplitudes = %(number_of_amplitudes)i;

    const unsigned int number_of_real_parameters = %(number_of_real_parameters)i;
    const std::vector<std::string> names_of_real_parameters = {%(names_of_real_parameters)s};

    const unsigned int number_of_complex_parameters = %(number_of_complex_parameters)i;
    const std::vector<std::string> names_of_complex_parameters = {%(names_of_complex_parameters)s};

    const unsigned int number_of_regulators = %(number_of_regulators)i;
    const std::vector<std::string> names_of_regulators = {%(names_of_regulators)s};

    const std::vector<int> requested_orders = {%(requested_orders)s};
    // --}
    
    // basic data types
    // --{
    typedef double real_t;
    #ifdef SECDEC_WITH_CUDA
        typedef thrust::complex<real_t> complex_t;
    #else
        typedef std::complex<real_t> complex_t;
    #endif
    
    const unsigned int maximal_number_of_integration_variables = %(number_of_integration_variables)i;

    // all integrals must have the same return type, assume complex
    typedef complex_t integrand_return_t;
    template<typename T> using nested_series_t = %(nested_series_type)s;
    template<typename T> using amplitudes_t = std::vector<nested_series_t<T>>;
    typedef secdecutil::IntegrandContainer<integrand_return_t, real_t const * const, real_t> integrand_t;
    // --}

    // amplitude-related data types
    // --{
    typedef secdecutil::amplitude::Integral<integrand_return_t,real_t> integral_t;
    typedef secdecutil::amplitude::WeightedIntegral<integral_t,integrand_return_t> weighted_integral_t;
    typedef std::vector<weighted_integral_t> sum_t;
    template<template<typename...> class container_t> using handler_t = secdecutil::amplitude::WeightedIntegralHandler<integrand_return_t,real_t,integrand_return_t,container_t>;
    
    #ifdef SECDEC_WITH_CUDA
        #if %(name)s_contour_deformation
            typedef secdecutil::CudaIntegrandContainerWithDeformation
                    <
                        real_t,complex_t,1/*maximal_number_of_functions*/,
                        maximal_number_of_integration_variables,
                        number_of_real_parameters,number_of_complex_parameters,
                        %(name_as_char_pack)s
                    >
                    cuda_integrand_t;
            typedef cuda_integrand_t cuda_together_integrand_t;
        #else
            typedef secdecutil::CudaIntegrandContainerWithoutDeformation
                    <
                        real_t,complex_t,integrand_return_t,
                        1/*maximal_number_of_functions*/,
                        number_of_real_parameters,number_of_complex_parameters,
                        %(name_as_char_pack)s
                    >
                    cuda_integrand_t;
            typedef cuda_integrand_t cuda_together_integrand_t;
        #endif
        typedef cuda_integrand_t user_integrand_t;
    #else
        typedef integrand_t user_integrand_t;
    #endif
    // --}

    // amplitude getter functions
    // --{
    template<typename integrator_t>
    std::vector<nested_series_t<sum_t>> make_amplitudes
    (
        const std::vector<real_t>& real_parameters,
        const std::vector<complex_t>& complex_parameters,
        const std::string& lib_path,
        const integrator_t& integrator
        #if %(name)s_contour_deformation
            ,unsigned number_of_presamples = 100000,
            real_t deformation_parameters_maximum = 1.,
            real_t deformation_parameters_minimum = 1.e-5,
            real_t deformation_parameters_decrease_factor = 0.9
        #endif
    );
    std::vector<nested_series_t<sum_t>> make_amplitudes
    (
        const std::vector<real_t>& real_parameters,
        const std::vector<complex_t>& complex_parameters,
        const std::string& lib_path,
        const secdecutil::Integrator<integrand_return_t,real_t> * integrator
        #if %(name)s_contour_deformation
            ,unsigned number_of_presamples,
            real_t deformation_parameters_maximum,
            real_t deformation_parameters_minimum,
            real_t deformation_parameters_decrease_factor
        #endif
    );
    #ifdef SECDEC_WITH_CUDA
        std::vector<nested_series_t<sum_t>> make_amplitudes //<const secdecutil::Integrator<integrand_return_t,real_t,cuda_integrand_t>*>
        (
            const std::vector<real_t>& real_parameters,
            const std::vector<complex_t>& complex_parameters,
            const std::string& lib_path,
            const secdecutil::Integrator<integrand_return_t,real_t,cuda_integrand_t> * integrator
            #if %(name)s_contour_deformation
                ,unsigned number_of_presamples,
                real_t deformation_parameters_maximum,
                real_t deformation_parameters_minimum,
                real_t deformation_parameters_decrease_factor
            #endif
        );
    #endif
    // --}
};

#ifdef SECDEC_WITH_CUDA
// provide custom common_type (together + integral) = together
namespace std 
{
    template <>
    struct common_type<%(name)s::cuda_together_integrand_t, %(name)s::cuda_integrand_t> 
    {
        using type = %(name)s::cuda_together_integrand_t;
    };
}
#endif

#endif

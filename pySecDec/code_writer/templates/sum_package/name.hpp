#ifndef %(name)s_hpp_included
#define %(name)s_hpp_included

#include <vector> // std::vector

#include <secdecutil/amplitude.hpp> // secdecutil::amplitude::Integral, secdecutil::amplitude::WeightedIntegral, secdecutil::amplitude::WeightedIntegralHandler
#include <secdecutil/series.hpp> // secdecutil::Series

namespace %(name)s
{
    // basic data types
    // --{
    typedef double real_t;
    #ifdef SECDEC_WITH_CUDA
        typedef thrust::complex<real_t> complex_t;
    #else
        typedef std::complex<real_t> complex_t;
    #endif

    // all integrals must have the same return type, assume complex
    typedef complex_t integrand_return_t;
    template<typename T> using nested_series_t = %(nested_series_type)s;
    // --}

    // amplitude-related data types
    // --{
    typedef secdecutil::amplitude::Integral<integrand_return_t,real_t> integral_t;
    typedef secdecutil::amplitude::WeightedIntegral<integral_t,integrand_return_t> weighted_integral_t;
    typedef std::vector<weighted_integral_t> sum_t;
    template<template<typename...> class container_t> using handler_t = secdecutil::amplitude::WeightedIntegralHandler<integrand_return_t,real_t,integrand_return_t,container_t>;
    // --}

    // amplitude getter functions
    // --{
    nested_series_t<sum_t> make_amplitude
    (
        const std::vector<real_t>& real_parameters,
        const std::vector<complex_t>& complex_parameters
    );
    // --}

    // some information about the integral
    // --{
    const unsigned long long number_of_integrals = %(number_of_integrals)i;

    const unsigned int number_of_regulators = %(number_of_regulators)i;
    const std::vector<std::string> names_of_regulators = {%(names_of_regulators)s};

    const std::vector<int> requested_orders = {%(requested_orders)s};
    // --}
};
#endif

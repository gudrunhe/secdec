#ifndef %(name)s_hpp_included
#define %(name)s_hpp_included

#include <complex>
#include <string>
#include <vector>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/sector_container.hpp>
#include <secdecutil/series.hpp>
#include <secdecutil/uncertainties.hpp>

namespace %(name)s
{
    // whether or not to use contour deformation
    #define %(name)s_contour_deformation %(contour_deformation)i

    // whether or not complex parameters are present
    #define %(name)s_has_complex_parameters %(have_complex_parameters)i

    // basic data types
    // --{
    typedef double real_t;
    typedef std::complex<real_t> complex_t;
    #if %(name)s_has_complex_parameters || %(name)s_contour_deformation
        typedef complex_t integrand_return_t;
    #else
        typedef real_t integrand_return_t;
    #endif
    // --}

    const unsigned int number_of_sectors = %(number_of_sectors)i;

    const unsigned int number_of_regulators = %(number_of_regulators)i;
    const std::vector<std::string> names_of_regulators = {%(names_of_regulators)s};

    const unsigned int number_of_real_parameters = %(number_of_real_parameters)i;
    const std::vector<std::string> names_of_real_parameters = {%(names_of_real_parameters)s};

    const unsigned int number_of_complex_parameters = %(number_of_complex_parameters)i;
    const std::vector<std::string> names_of_complex_parameters = {%(names_of_complex_parameters)s};

    const std::vector<int> lowest_orders = {%(lowest_orders)s}; // not including the prefactor
    const std::vector<int> highest_orders = {%(highest_orders)s}; // not including the prefactor
    const std::vector<int> lowest_prefactor_orders = {%(lowest_prefactor_orders)s};
    const std::vector<int> highest_prefactor_orders = {%(highest_prefactor_orders)s};
    const std::vector<int> requested_orders = {%(requested_orders)s};

    extern const std::vector<%(sector_container_type)s> sectors;
    %(prefactor_type)s  prefactor(const std::vector<real_t>& real_parameters, const std::vector<complex_t>& complex_parameters);

    %(make_integrands_return_t)s make_integrands
    (
        const std::vector<real_t>& real_parameters,
        const std::vector<complex_t>& complex_parameters
        #if %(name)s_contour_deformation
            ,unsigned number_of_samples = 100000,
            real_t deformation_parameters_maximum = 1.,
            real_t deformation_parameters_minimum = 1.e-5,
            real_t deformation_parameters_decrease_factor = 0.9,
            real_t deformation_parameters_deformation_offset = 1.e-3
        #endif
    );

    #undef %(name)s_contour_deformation
    #undef %(name)s_has_complex_parameters

};
#endif

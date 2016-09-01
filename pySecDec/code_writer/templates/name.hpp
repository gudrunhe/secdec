#ifndef %(name)s_hpp_included
#define %(name)s_hpp_included

#include <cmath>
#include <complex>
#include <limits>
#include <vector>
#include <secdecutil/deep_apply.hpp>
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
    const unsigned int number_of_regulators = %(number_of_regulators)i; //TODO: names of regulators
    const unsigned int number_of_real_parameters = %(number_of_real_parameters)i; //TODO: names of real_parameters
    const unsigned int number_of_complex_parameters = %(number_of_complex_parameters)i; //TODO: names of complex_parameters
    const std::vector<int> lowest_orders = {%(lowest_orders)s}; // not including the prefactor // TODO: lowest_prefactor_orders
    const std::vector<int> highest_orders = {%(highest_orders)s}; // not including the prefactor // TODO: highest_prefactor_orders
    const std::vector<int> requested_orders = {%(requested_orders)s};
    extern const std::vector<%(sector_container_type)s> sectors;
    %(prefactor_type)s  prefactor(const std::vector<real_t>& real_parameters, const std::vector<complex_t>& complex_parameters);

    auto make_integrands
    (
        const std::vector<real_t>& real_parameters,
        const std::vector<complex_t>& complex_parameters
        #if %(name)s_contour_deformation
            ,unsigned number_of_samples = 100000,
            real_t deformation_parameters_maximum = 1.,
            real_t deformation_parameters_minimum = 1.e-5,
            real_t deformation_parameters_decrease_factor = 0.9
        #endif
    )
    -> decltype
    (
        #if %(name)s_contour_deformation
            secdecutil::deep_apply( sectors, secdecutil::SectorContainerWithDeformation_to_IntegrandContainer(real_parameters, complex_parameters) )
        #else
            secdecutil::deep_apply( sectors, secdecutil::SectorContainerWithoutDeformation_to_IntegrandContainer<integrand_return_t>(real_parameters, complex_parameters) )
        #endif
    );

    #undef %(name)s_contour_deformation
    #undef %(name)s_has_complex_parameters


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO: move everything below to a different file

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    constexpr complex_t i_{0,1}; // the imaginary unit

    // required special functions
    // define your own additionally needed special function here
    // --{

    /*
     * We need the nonstandard continuation on the negative x-axis;
     * i.e. log(-1) = -i*pi --> subtract a tiny imaginary part.
     */
    // use std::log(real_t) but override log(complex_t)
    using std::log;
    inline complex_t log(complex_t arg)
    {
        return std::log( arg - i_*std::numeric_limits<real_t>::min() );
    }

    /*
     * We do not want to use "std::pow(double, int)" because the g++ compiler
     * casts it to "pow(double, double)" which is extremely inefficient but
     * demanded by the c++ standard.
     *
     * Note: Using std::pow and the flags "-O2" and "-ffast-math" with g++ is even faster,
     *       but "-ffast-math" is g++ specific and allows "unsafe math optimizations".
     *       The intel compiler produces code that runs faster when using std::pow and
     *       "-O2" than with this function.
     *       However, the c++ standard requires that the second argument of
     *       "std::pow(double, int)" is casted to double. To comply with the standard, we
     *       decided to implement our own optimized power function rather than relying on
     *       the compiler to perform optimizations possibly disallowed by the c++ standard.
     *       Playing around with "std::pow" and the aforementioned switches is nevertheless
     *       worth a try in practical applications where high performance is needed.
     */
    template <typename T> inline T pow(T base, int exponent)
    {
        if (exponent < 0)
            return 1./pow(base, -exponent);

        else if (exponent == 0)
            return 1.;

        else if (exponent == 1)
            return base;

        else if (exponent == 2)
            return base * base;

        else if (exponent == 3)
            return base * base * base;

        else if (exponent == 4)
        {
            T result = base;
            result *= result;
            result *= result;
            return result;
        }

        else if (exponent == 5)
        {
            T result = base;
            result *= result;
            result *= result;
            return result * base;
        }

        else if (exponent == 6)
        {
            T result = base * base * base;
            return result * result;
        }

        else if (exponent == 7)
        {
            T result = base * base * base;
            return result * result * base;
        }

        else if (exponent == 8)
        {
            T result = base;
            result *= result;
            result *= result;
            return result * result;
        }

        else if (exponent == 16)
        {
            T tmp = base * base;
            tmp *= tmp;
            tmp *= tmp;
            tmp *= tmp;
            return tmp;
        }

        unsigned half_exponent = exponent / 2;
        T out = pow(base, half_exponent);

        out *= out;
        if (2 * half_exponent == exponent) // exponent is even
            return out;
        else // exponent is odd --> need another factor of the base due to integer division above
            return out * base;
    }
    // --}

};
#endif

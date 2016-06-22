#ifndef __SecDec_include_guard_%(name)s_config
#define __SecDec_include_guard_%(name)s_config

#include <array>
#include <cmath>
#include <complex>
#include <exception>
#include <functional>
#include <limits>
#include <vector>
// TODO: which of these headers are needed?

namespace secdecutil {}; // TODO: include external library "secdecutil" instead of this

namespace %(name)s
{
    using namespace secdecutil;

    // whether or not to use contour deformation
    // TODO: set this variable in the Makefile
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

    // function types
    // --{
    // the call signature of an integrand
    typedef integrand_return_t IntegrandFunction(real_t const * const integration_variables, real_t const * const real_parameters, complex_t const * const complex_parameters);

    #if %(name)s_contour_deformation

        // the call signature of an integrand to be deformed
        typedef complex_t DeformableIntegrandFunction(complex_t const * const transformed_integration_variables, real_t const * const real_parameters, complex_t const * const complex_parameters);

        // the return type of the integral transformation (contour deformation)
        struct integral_transformation_t
        {
            std::vector<complex_t> transformed_variables;
            complex_t Jacobian_determinant;
        };

        // the call signature of the integral transformation (contour deformation)
        typedef integral_transformation_t ContourDeformationFunction
        (
            real_t const * const integration_variables,
            real_t const * const deformation_parameters,
            real_t const * const real_parameters,
            complex_t const * const complex_parameters,
        );

    #endif

    // --}

    // this error is thrown if the sign check of the deformation (contour_deformation_polynomial.imag() <= 0) fails
    #if %(name)s_contour_deformation
        struct sign_check_error : public std::runtime_error
        {
            sign_check_error (const std::string& what_arg) : std::runtime_error(what_arg) {};
            sign_check_error (const char*        what_arg) : std::runtime_error(what_arg) {};
        };
    #endif

    // the imaginary unit
    constexpr complex_t i_{0,1};

    // container class to collect the integrand functions
    struct IntegrandContainer
    {
        const int number_of_integration_variables;

        #if %(name)s_contour_deformation

            const std::function<IntegrandFunction> integrand;

        #else

            const std::function<DeformableIntegrandFunction> undeformed_integrand;
            const std::function<DeformableIntegrandFunction> contour_deformation_polynomial;
            const std::function<ContourDeformationFunction> contour_deformation;

            void optimize_deformation_parameters
            (
                const real_t const * const initial_guess,
                real_t const * const real_parameters,
                complex_t const * const complex_parameters,
                const size_t number_of_samples = 1000,
                const double rel_error = 1.e-2,
                const double abs_error = 1.e-7,
                const int maxiter = 1000
            ); // TODO: what arguments are needed?

            inline complex_t deformed_integrand
            (
                real_t const * const integration_variables,
                real_t const * const real_parameters,
                complex_t const * const complex_parameters,
                deformation_parameter_t const * const deformation_parameters
            ) const
            {
                auto deformation = contour_deformation(integration_variables, Mandelstam, mass, deformation_parameters);
                if (F(deformation.transformed_variables.data(), Mandelstam, mass).imag() > 0.)
                    throw sign_check_error("Contour deformation yields the wrong sign of \"contour_deformation_polynomial.imag\". Choose smaller \"deformation_parameters.\"");
                return deformation.Jacobian_determinant * integrand(deformation.transformed_variables.data(), Mandelstam, mass);
            };
            inline complex_t integrand
            (
                integration_variable_t const * const integration_variables,
                Mandelstam_t const * const Mandelstam,
                mass_t const * const mass
            ) const
            {
                auto deformation = contour_deformation(integration_variables, Mandelstam, mass, optimized_deformation_parameters.data());
                if (F(deformation.transformed_variables.data(), Mandelstam, mass).imag() > 0.)
                    throw sign_check_error("Contour deformation yields the wrong sign of \"contour_deformation_polynomial.imag\". Choose a larger \"number_of_samples\" in \"optimize_deformation_parameters.\"");
                return deformation.Jacobian_determinant * integrand(deformation.transformed_variables.data(), Mandelstam, mass);
            };

            // need an explicit constructor due to the private members
            integrand_t
            (
                const int&& number_of_Feynman_parameters,
                const IntegrandFunction&& integrand,
                const ContourDeformationFunction&& contour_deformation,
                const IntegrandFunction&& contour_deformation_polynomial
            ) : number_of_Feynman_parameters(number_of_Feynman_parameters),
                integrand(integrand),
                contour_deformation(contour_deformation),
                contour_deformation_polynomial(contour_deformation_polynomial),
                optimized_deformation_parameters(number_of_Feynman_parameters, 1.0) // fill constructor (fill with ones)
            {};

            inline const std::vector<real_t> get_optimized_deformation_parameters() const
            {
                return optimized_deformation_parameters;
            };

            private:

                std::vector<real_t> optimized_deformation_parameters;

        #endif
    };

    // TODO: move special functions to an other file
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
    template <typename T> inline T pow(T base, unsigned exponent)
    {
        if (exponent == 0)
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

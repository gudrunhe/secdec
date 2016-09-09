#ifndef %(name)s_src_SecDecInternalFunctions_hpp_included
#define %(name)s_src_SecDecInternalFunctions_hpp_included

#include "%(name)s.hpp"

#define %(name)s_contour_deformation %(contour_deformation)i
#if %(name)s_contour_deformation
    #include "contour_deformation_functions.hpp"
#endif
#undef %(name)s_contour_deformation

#include <cmath>
#include <complex>

namespace %(name)s
{
    constexpr complex_t i_{0,1}; // the imaginary unit

    // required functions
    // --{

    /*
     * We need the nonstandard continuation on the negative x-axis;
     * i.e. log(-1) = -i*pi.
     */
    // use std::log(real_t) but override log(complex_t)
    #define %(name)s_contour_deformation %(contour_deformation)i
    #define %(name)s_has_complex_parameters %(have_complex_parameters)i
    #if %(name)s_has_complex_parameters || %(name)s_contour_deformation
        inline complex_t log(complex_t arg)
        {
            if (arg.imag() == 0)
                arg = complex_t(arg.real(),-0.);
            return std::log(arg);
        }
    #else
        inline real_t log(real_t arg)
        {
            return std::log(arg);
        }
    #endif
    #undef %(name)s_contour_deformation
    #undef %(name)s_has_complex_parameters

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
    using std::pow;
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

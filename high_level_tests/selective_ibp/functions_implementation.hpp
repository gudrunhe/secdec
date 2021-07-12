#ifndef userdefined_cpp_integral_src_functions_hpp_included
#define userdefined_cpp_integral_src_functions_hpp_included

#include "userdefined_cpp_integral.hpp"
#include "SecDecInternalFunctions.hpp"

#include <cmath>
#include <complex>

namespace userdefined_cpp_integral
{

    /*
     * Declarations of the `functions` and their required
     * derivatives are declared here. The derivative of a function
     * 'f' by its i-th argument is denoted as 'dfdi'. To implement
     * the functions listed below, you can either add "inline"
     * keywords and define these functions here, or you define the
     * functions in a separate '.cpp' file. If you decide for a
     * separate file, the file name can be arbitrary up to the
     * '.cpp' suffix. Furthermore, the '.cpp' file must be located
     * in this directory ('src/'). If you want to link against
     * an external library (e.g. the gsl), you should add the
     * corresponding compiler and linker flags to the "Makefile.conf"
     * in the top level directory.
     *
     * Note: Not all functions listed here may actually be needed.
     *       This file lists all derivatives that occurred in the
     *       calculation. It is possible that some dropped out due
     *       to algebraic simplifications after this list was
     *       generated.
     */

    template<typename T>
    integrand_return_t func(T arg)
    {
        return 1.0;
    }

    template<typename T>
    integrand_return_t dfuncd0(T arg)
    {
        return 0.0;
    }

    inline
    integrand_return_t HeavisideTheta(real_t x)
    {
        return x > 0.0 ? 1.0 : 0.0;
    }

};
#endif

#ifndef userdefined_cpp_src_functions_hpp_included
#define userdefined_cpp_src_functions_hpp_included

#include "userdefined_cpp.hpp"
#include "SecDecInternalFunctions.hpp"

#include <cmath>
#include <complex>

namespace userdefined_cpp
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

/*
    template<typename T0, typename T1, typename T2>
    integrand_return_t dum2(T0 arg0, T1 arg1, T2 arg2);

    template<typename T0, typename T1, typename T2, typename T3>
    integrand_return_t dum1(T0 arg0, T1 arg1, T2 arg2, T3 arg3);

*/

};
#endif

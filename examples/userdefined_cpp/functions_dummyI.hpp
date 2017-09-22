#ifndef dummyI_src_functions_hpp_included
#define dummyI_src_functions_hpp_included

#include "dummyI.hpp"
#include "SecDecInternalFunctions.hpp"

#include <cmath>
#include <complex>

namespace dummyI
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


    template<typename T0, typename T1, typename T2>
    integrand_return_t dum2(T0 arg0, T1 arg1, T2 arg2)
    {
        return  arg0*arg0 + arg1*arg1 +arg2*arg2 + 4*arg0*arg1 +3*arg0*arg0*arg1*arg1 - sqrt(arg0*arg1*arg2);
    };
    
    template<typename T0, typename T1, typename T2, typename T3>
    integrand_return_t dum1(T0 arg0, T1 arg1, T2 arg2, T3 arg3)
    {
        return arg0*arg0 + pow(arg1,3) +pow(arg2,4) +pow(arg3,5) + 4*arg0*arg1*arg2*arg3 +2 -arg0*arg0 *pow(arg1,3) *pow(arg2,4) *pow(arg3,5);
    };

    template<typename T0, typename T1, typename T2, typename T3>
    integrand_return_t ddum1d0(T0 arg0, T1 arg1, T2 arg2, T3 arg3);

    template<typename T0, typename T1, typename T2>
    integrand_return_t ddum2d0(T0 arg0, T1 arg1, T2 arg2);



};
#endif

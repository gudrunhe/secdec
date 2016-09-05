#ifndef %(name)s_src_functions_hpp_included
#define %(name)s_src_functions_hpp_included

#include "%(name)s.hpp"
#include "SecDecInternalFunctions.hpp"

#include <cmath>
#include <complex>

namespace %(name)s
{

    /*
     * Declarations of the `functions` and their required
     * derivatives are declared here. You can either add "inline"
     * keywords and define these functions here, or you define the
     * functions in a separate '.cpp' file. If you decide for a
     * separate file, the file name can be arbitrary up to the
     * '.cpp' suffix. Furthermore, the '.cpp' file must be located
     * in this directory ('src/'). If you want to link against
     * an external library (e.g. the gsl), you should add the
     * corresponding compiler and linker flags to the "Makefile.conf"
     * in the top level directory.
     */

/*
%(function_declarations)s
*/

};
#endif

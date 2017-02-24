#define INTEGRAL_NAME %(name)s

// whether or not to use contour deformation
#define integral_contour_deformation %(contour_deformation)i

// whether or not complex parameters are present
#define integral_has_complex_parameters %(have_complex_parameters)i

// whether or no the return type should be complex in any case
#define integral_enforce_complex_return_type %(enforce_complex_return_type)i



#include "%(name)s.hpp"



// The python-C binding is general and therefore contained in the util
#include <secdecutil/pylink.hpp>



#undef integral_contour_deformation
#undef integral_has_complex_parameters
#undef integral_enforce_complex_return_type

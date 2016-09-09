#ifndef %(name)s_src_contour_deformation_functions_hpp_included
#define %(name)s_src_contour_deformation_functions_hpp_included

#include <cmath>

namespace %(name)s
{
    // required functions for contour deformation

    /*
     * Numerically stable implementation of x*exp(-mu/x)
     * for x in [0,1].
     */
    inline real_t SecDecInternalXExpMinusMuOverX(real_t mu, real_t x)
    {
        if (x == 0)
            return 0;

        // x*exp(-mu/x) = exp(-mu/x + log(x))
        return std::exp(  -mu/x + std::log(x)  );
    }

    /*
     * Numerically stable implementation of the first
     * derivative of x*exp(-mu/x) by x for x in [0,1].
     */
    inline real_t dSecDecInternalXExpMinusMuOverXd1(real_t mu, real_t x)
    {
        if (x == 0)
            return 0;

        return std::exp(-mu/x) + mu * std::exp(  -mu/x - std::log(x)  );
    }

    /*
     * Numerically stable implementation of exp(-mu/x)
     * for x in [0,1].
     */
    inline real_t SecDecInternalExpMinusMuOverX(real_t mu, real_t x)
    {
        if (x == 0)
            return 0;
        return std::exp(-mu/x);
    }

};
#endif

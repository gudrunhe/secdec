#ifndef %(name)s_src_contour_deformation_functions_hpp_included
#define %(name)s_src_contour_deformation_functions_hpp_included

#include <cmath>
#include <complex>

namespace %(name)s
{
    // required functions for contour deformation

    inline real_t SecDecInternalXExpMinusMuOverX(real_t mu, real_t x)
    {
        if (x == 0)
            return 0;

        return x * std::exp(-mu/x);
    }

    // The first derivative of x*exp(-mu/x) by x.
    inline real_t dSecDecInternalXExpMinusMuOverXd1(real_t mu, real_t x)
    {
        if (x == 0)
            return 0;

        real_t minus_mu_over_x = -mu/x;
        return std::exp(minus_mu_over_x) - minus_mu_over_x * std::exp(minus_mu_over_x);
    }

    inline real_t SecDecInternalExpMinusMuOverX(real_t mu, real_t x)
    {
        if (x == 0)
            return 0;

        return std::exp(-mu/x);
    }

    inline real_t SecDecInternalRealPart(complex_t x)
    {
        return std::real(x);
    }

};
#endif

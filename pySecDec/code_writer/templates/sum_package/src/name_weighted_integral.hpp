#ifndef %(sub_integral_name)s_weighted_integral_hpp_included
#define %(sub_integral_name)s_weighted_integral_hpp_included

#include <vector> // std::vector

#include "%(name)s.hpp"
#include "%(sub_integral_name)s/%(sub_integral_name)s.hpp"

namespace %(name)s
{
    namespace %(sub_integral_name)s
    {
        nested_series_t<sum_t> make_weighted_integral
        (
            const std::vector<real_t>& real_parameters,
            const std::vector<complex_t>& complex_parameters
        );
    }
};
#endif

#include <vector> // std::vector

#include "%(name)s.hpp"

%(weighted_integral_includes)s

namespace %(name)s
{
    nested_series_t<sum_t> make_amplitude
    (
        const std::vector<real_t>& real_parameters,
        const std::vector<complex_t>& complex_parameters
    )
    {
        nested_series_t<sum_t>
        %(weighted_integral_sum_initialization)s
        return amplitude;
    };
};

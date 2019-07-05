#include <vector> // std::vector

#include "%(name)s.hpp"

%(weighted_integral_includes)s

namespace %(name)s
{
    std::vector<nested_series_t<sum_t>> make_amplitudes
    (
        const std::vector<real_t>& real_parameters,
        const std::vector<complex_t>& complex_parameters
    )
    {
        std::vector<nested_series_t<sum_t>> amplitudes;
        amplitudes.reserve(number_of_amplitudes);

        for (unsigned int amp_idx = 0; amp_idx < number_of_amplitudes; ++amp_idx)
        {
            nested_series_t<sum_t>
            %(weighted_integral_sum_initialization)s

            amplitudes.push_back(amplitude);
        }

        return amplitudes;
    };
};

#include <secdecutil/series.hpp>
#include <secdecutil/uncertainties.hpp>
#include <vector>

#include "%(name)s.hpp"

namespace %(name)s
{
    %(prefactor_type)s prefactor(const std::vector<real_t>& real_parameters, const std::vector<complex_t>& complex_parameters)
    {
        %(prefactor_function_body)s
    }
};

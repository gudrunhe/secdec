#include <secdecutil/series.hpp>
#include <secdecutil/uncertainties.hpp>
#include <vector>

#include "%(name)s.hpp"
#include "functions.hpp"

namespace %(name)s
{
    nested_series_t<integrand_return_t> prefactor(const std::vector<real_t>& real_parameters, const std::vector<complex_t>& complex_parameters)
    {
        %(prefactor_function_body)s
    }
};

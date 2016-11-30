#include <secdecutil/deep_apply.hpp>
#include <secdecutil/sector_container.hpp>
#include <secdecutil/series.hpp>
#include <string>
#include <vector>

#include "%(name)s.hpp"
%(sector_includes)s

namespace %(name)s
{
    const std::vector<nested_series_t<sector_container_t>> sectors = {%(sectors_initializer)s};

    #define %(name)s_contour_deformation %(contour_deformation)i
    std::vector<nested_series_t<secdecutil::IntegrandContainer<integrand_return_t, real_t const * const>>> make_integrands
    (
        const std::vector<real_t>& real_parameters,
        const std::vector<complex_t>& complex_parameters
        #if %(name)s_contour_deformation
            ,unsigned number_of_samples,
            real_t deformation_parameters_maximum,
            real_t deformation_parameters_minimum,
            real_t deformation_parameters_decrease_factor
        #endif
    )
    {
        if ( real_parameters.size() != %(name)s::number_of_real_parameters )
            throw std::logic_error(
                                        "Called \"%(name)s::make_integrands\" with " +
                                        std::to_string(real_parameters.size()) + " \"real_parameters\" (" +
                                        std::to_string(%(name)s::number_of_real_parameters) + " expected)."
                                  );

        if ( complex_parameters.size() != %(name)s::number_of_complex_parameters )
            throw std::logic_error(
                                        "Called \"%(name)s::make_integrands\" with " +
                                        std::to_string(complex_parameters.size()) + " \"complex_parameters\" (" +
                                        std::to_string(%(name)s::number_of_complex_parameters) + " expected)."
                                  );

        #if %(name)s_contour_deformation
            return secdecutil::deep_apply
            (
                sectors,
                secdecutil::SectorContainerWithDeformation_to_IntegrandContainer
                    (
                        real_parameters,
                        complex_parameters,
                        number_of_samples,
                        deformation_parameters_maximum,
                        deformation_parameters_minimum,
                        deformation_parameters_decrease_factor
                    )
            );
        #else
            return secdecutil::deep_apply( sectors, secdecutil::SectorContainerWithoutDeformation_to_IntegrandContainer<integrand_return_t>(real_parameters, complex_parameters) );
        #endif
    };
    #undef %(name)s_contour_deformation
};

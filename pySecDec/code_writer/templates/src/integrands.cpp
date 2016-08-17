#include <secdecutil/sector_container.hpp>
#include <secdecutil/series.hpp>
#include <vector>

#include "%(name)s.hpp"
%(sector_includes)s

namespace %(name)s
{
    const std::vector<%(sector_container_type)s> sectors = {%(sectors_initializer)s};

    #define %(name)s_contour_deformation %(contour_deformation)i
    auto make_integrands
    (
        const std::vector<real_t>& real_parameters,
        const std::vector<complex_t>& complex_parameters
        #if %(name)s_contour_deformation
            ,unsigned number_of_samples,
            real_t deformation_parameters_initial_guess
        #endif
    )
    -> decltype
    (
        #if %(name)s_contour_deformation
            secdecutil::deep_apply( sectors, secdecutil::SectorContainerWithDeformation_to_IntegrandContainer(real_parameters, complex_parameters) )
        #else
            secdecutil::deep_apply( sectors, secdecutil::SectorContainerWithoutDeformation_to_IntegrandContainer<%(name)s::integrand_return_t>(real_parameters, complex_parameters) )
        #endif
    )
    {
        #if %(name)s_contour_deformation
            return secdecutil::deep_apply
            (
                sectors,
                secdecutil::SectorContainerWithDeformation_to_IntegrandContainer
                    (
                        real_parameters,
                        complex_parameters,
                        number_of_samples,
                        deformation_parameters_initial_guess
                    )
            );
        #else
            return secdecutil::deep_apply( sectors, secdecutil::SectorContainerWithoutDeformation_to_IntegrandContainer<integrand_return_t>(real_parameters, complex_parameters) );
        #endif
    };
    #undef %(name)s_contour_deformation
};

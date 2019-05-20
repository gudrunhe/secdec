#include <cassert> // assert
#include <fstream> // std::ifstream
#include <memory> // std::shared_ptr, std::make_shared
#include <numeric> // std::accumulate
#include <vector> // std::vector

#include <secdecutil/amplitude.hpp> // secdecutil::amplitude::Integral, secdecutil::amplitude::CubaIntegral, secdecutil::amplitude::QmcIntegral
#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply
#include <secdecutil/ginac_coefficient_parser.hpp> // secdecutil::ginac::read_coefficient
#include <secdecutil/integrators/cuba.hpp> // secdecutil::cuba::Vegas, secdecutil::cuba::Suave, secdecutil::cuba::Cuhre, secdecutil::cuba::Divonne
#include <secdecutil/integrators/qmc.hpp> // secdecutil::integrators::Qmc
#include <secdecutil/series.hpp> // secdecutil::Series

#include "%(name)s.hpp"
#include "%(sub_integral_name)s/%(sub_integral_name)s.hpp"
#include "%(sub_integral_name)s_weighted_integral.hpp"
#include "config_%(name)s.hpp"


namespace %(sub_integral_name)s
{
    const std::vector<int> lowest_coefficient_orders{%(lowest_coefficient_orders)s};

    static std::vector<int> compute_required_orders()
    {
        assert(%(name)s::requested_orders.size() == %(sub_integral_name)s::requested_orders.size());
        size_t number_of_regulators = %(name)s::requested_orders.size();
        std::vector<int> orders; orders.reserve( %(name)s::requested_orders.size() );
        for(size_t i = 0; i < number_of_regulators; ++i)
            orders.push_back(
                                  %(name)s::requested_orders.at(i) + 1
                                - %(sub_integral_name)s::lowest_coefficient_orders.at(i)
                                - %(sub_integral_name)s::lowest_orders.at(i)
                                - %(sub_integral_name)s::lowest_prefactor_orders.at(i)
                            );
        return orders;
    }

    nested_series_t<complex_t> coefficient(const std::vector<real_t>& real_parameters, const std::vector<complex_t>& complex_parameters)
    {
        std::ifstream coeffile("lib/%(sub_integral_name)s_coefficient.txt");
        assert( coeffile.is_open() );
        return secdecutil::ginac::read_coefficient<nested_series_t>
               (
                    coeffile, compute_required_orders(),
                    names_of_regulators, names_of_real_parameters, names_of_complex_parameters,
                    real_parameters, complex_parameters
               );
    }
};

namespace %(name)s
{
    namespace %(sub_integral_name)s
    {
        template<bool with_cuda>
        struct WithCuda
        {
            using integrand_t = ::%(sub_integral_name)s::integrand_t;
            constexpr static auto make_raw_integrands = ::%(sub_integral_name)s::make_integrands;

            // with contour deformation
            static std::vector<nested_series_t<integrand_t>> make_integrands
            (
                const std::vector<real_t>& real_parameters,
                const std::vector<complex_t>& complex_parameters,
                std::vector<nested_series_t<integrand_t>> (*make_integrands_with_contour_deformation)
                    (
                        const std::vector<real_t>& real_parameters,
                        const std::vector<complex_t>& complex_parameters,
                        unsigned number_of_presamples,
                        real_t deformation_parameters_maximum,
                        real_t deformation_parameters_minimum,
                        real_t deformation_parameters_decrease_factor
                    )
            )
            {
                return make_integrands_with_contour_deformation
                (
                    real_parameters,
                    complex_parameters,
                    Options::ContourDeformation::number_of_presamples,
                    Options::ContourDeformation::deformation_parameters_maximum,
                    Options::ContourDeformation::deformation_parameters_minimum,
                    Options::ContourDeformation::deformation_parameters_decrease_factor
                );
            };

            // without contour deformation
            static std::vector<nested_series_t<integrand_t>> make_integrands
            (
                const std::vector<real_t>& real_parameters,
                const std::vector<complex_t>& complex_parameters,
                std::vector<nested_series_t<integrand_t>> (*make_integrands_without_contour_deformation)
                    (
                        const std::vector<real_t>& real_parameters,
                        const std::vector<complex_t>& complex_parameters
                    )
            )
            {
                return make_integrands_without_contour_deformation
                (
                    real_parameters,
                    complex_parameters
                );
            };
        };

        #ifdef SECDEC_WITH_CUDA
            template<>
            struct WithCuda<true>
            {
                using integrand_t = ::%(sub_integral_name)s::cuda_integrand_t;
                constexpr static auto make_raw_integrands = ::%(sub_integral_name)s::make_cuda_integrands;

                // with contour deformation
                static std::vector<nested_series_t<integrand_t>> make_integrands
                (
                    const std::vector<real_t>& real_parameters,
                    const std::vector<complex_t>& complex_parameters,
                    std::vector<nested_series_t<integrand_t>> (*make_integrands_with_contour_deformation)
                        (
                            const std::vector<real_t>& real_parameters,
                            const std::vector<complex_t>& complex_parameters,
                            unsigned number_of_presamples,
                            real_t deformation_parameters_maximum,
                            real_t deformation_parameters_minimum,
                            real_t deformation_parameters_decrease_factor
                        )
                )
                {
                    return make_integrands_with_contour_deformation
                    (
                        real_parameters,
                        complex_parameters,
                        Options::ContourDeformation::number_of_presamples,
                        Options::ContourDeformation::deformation_parameters_maximum,
                        Options::ContourDeformation::deformation_parameters_minimum,
                        Options::ContourDeformation::deformation_parameters_decrease_factor
                    );
                };

                // without contour deformation
                static std::vector<nested_series_t<integrand_t>> make_integrands
                (
                    const std::vector<real_t>& real_parameters,
                    const std::vector<complex_t>& complex_parameters,
                    std::vector<nested_series_t<integrand_t>> (*make_integrands_without_contour_deformation)
                        (
                            const std::vector<real_t>& real_parameters,
                            const std::vector<complex_t>& complex_parameters
                        )
                )
                {
                    return make_integrands_without_contour_deformation
                    (
                        real_parameters,
                        complex_parameters
                    );
                };
            };
        #endif

        nested_series_t<sum_t> make_weighted_integral
        (
            const std::vector<real_t>& real_parameters,
            const std::vector<complex_t>& complex_parameters
        )
        {
            using types = WithCuda<Options::cuda_compliant_integrator>;
            using integrand_t = types::integrand_t;
            using Integral = Options::Integral<
                                                   integrand_t,
                                                   ::%(sub_integral_name)s::maximal_number_of_integration_variables
                                              >;
            using integrator_t = Integral::integrator_t;

            const std::vector<nested_series_t<integrand_t>> raw_integrands =
            types::make_integrands
            (
                real_parameters,
                complex_parameters,
                types::make_raw_integrands
            );

            const std::shared_ptr<integrator_t> integrator = Integral::configure_integrator();

            const std::function<sum_t(const integrand_t& integrand)> convert_integrands =
                [ integrator ] (const integrand_t& integrand) -> sum_t
                {
                    return { /* constructor of std::vector */
                                { /* constructor of WeightedIntegral */
                                    std::make_shared<Integral::integral_t>(integrator,integrand)
                                }
                        };
                };

            const std::vector<nested_series_t<sum_t>> integrals = deep_apply(raw_integrands, convert_integrands);
            nested_series_t<sum_t> amplitude = std::accumulate(++integrals.begin(), integrals.end(), *integrals.begin() );
            amplitude *= ::%(sub_integral_name)s::prefactor(real_parameters,complex_parameters)
                       * ::%(sub_integral_name)s::coefficient(real_parameters,complex_parameters);
            return amplitude;
        }
    };
};

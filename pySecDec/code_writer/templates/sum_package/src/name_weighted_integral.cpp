#include <cassert> // assert
#include <fstream> // std::ifstream
#include <memory> // std::shared_ptr, std::make_shared
#include <numeric> // std::accumulate
#include <string> // std::string, std::to_string
#include <vector> // std::vector

#include <secdecutil/amplitude.hpp> // secdecutil::amplitude::Integral, secdecutil::amplitude::CubaIntegral, secdecutil::amplitude::QmcIntegral
#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply
#include <secdecutil/coefficient_parser.hpp> // secdecutil::exparse::read_coefficient
#include <secdecutil/integrators/cquad.hpp> // secdecutil::gsl::CQuad
#include <secdecutil/integrators/cuba.hpp> // secdecutil::cuba::Vegas, secdecutil::cuba::Suave, secdecutil::cuba::Cuhre, secdecutil::cuba::Divonne
#include <secdecutil/integrators/qmc.hpp> // secdecutil::integrators::Qmc
#include <secdecutil/series.hpp> // secdecutil::Series

#include "%(name)s.hpp"
#include "%(sub_integral_name)s/%(sub_integral_name)s.hpp"
#include "%(sub_integral_name)s_weighted_integral.hpp"

#define INTEGRAL_NAME %(name)s
#ifdef SECDEC_WITH_CUDA
    #define INTEGRAND_TYPE cuda_integrand_t
#else
    #define INTEGRAND_TYPE integrand_t
#endif

namespace %(sub_integral_name)s
{
    const std::vector<std::vector<int>> lowest_coefficient_orders{%(lowest_coefficient_orders)s};

    static std::vector<int> compute_required_orders(const unsigned int amp_idx)
    {
        assert(%(name)s::requested_orders.size() == %(sub_integral_name)s::requested_orders.size());
        size_t number_of_regulators = %(name)s::requested_orders.size();
        std::vector<int> orders; orders.reserve( %(name)s::requested_orders.size() );
        for(size_t i = 0; i < number_of_regulators; ++i)
        {
            int order = (
                              %(name)s::requested_orders.at(i) + 1
                            - %(sub_integral_name)s::lowest_coefficient_orders.at(amp_idx).at(i)
                            - %(sub_integral_name)s::lowest_orders.at(i)
                            - %(sub_integral_name)s::lowest_prefactor_orders.at(i)
                        );
            orders.push_back(order > 0 ? order : 1);
        }
        return orders;
    }

    nested_series_t<complex_t> coefficient(const std::vector<real_t>& real_parameters, const std::vector<complex_t>& complex_parameters, const unsigned int amp_idx, const std::string& lib_path)
    {
        std::string coeff_filename(lib_path + "/%(sub_integral_name)s_coefficient" + std::to_string(amp_idx) + ".txt");
        return secdecutil::exparse::read_coefficient<nested_series_t>
               (
                    coeff_filename, compute_required_orders(amp_idx),
                    names_of_regulators, names_of_real_parameters, names_of_complex_parameters,
                    real_parameters, complex_parameters
               );
    }
};

namespace %(name)s
{
    namespace %(sub_integral_name)s
    {
        typedef secdecutil::MultiIntegrator<INTEGRAL_NAME::integrand_return_t,INTEGRAL_NAME::real_t,INTEGRAL_NAME::integrand_t> multiintegrator_t;

        template<bool with_cuda>
        struct WithCuda
        {
            using integrand_t = ::%(sub_integral_name)s::integrand_t;

            // Note: we define make_integrands with %(name)s_contour_deformation
            // but call ::sub_integral_name::make_integrands with %(sub_integral_name)s_contour_deformation
            // i.e we drop contour deformation parameters not relevant for this integral
            static std::vector<nested_series_t<integrand_t>> make_integrands
            (
                const std::vector<real_t>& real_parameters,
                const std::vector<complex_t>& complex_parameters
                #if %(name)s_contour_deformation
                    ,unsigned number_of_presamples,
                    real_t deformation_parameters_maximum,
                    real_t deformation_parameters_minimum,
                    real_t deformation_parameters_decrease_factor
                #endif
            )
            {
                    return ::%(sub_integral_name)s::make_integrands
                    (
                        real_parameters,
                        complex_parameters
                        #if %(sub_integral_name)s_contour_deformation
                            ,number_of_presamples,
                            deformation_parameters_maximum,
                            deformation_parameters_minimum,
                            deformation_parameters_decrease_factor
                        #endif
                    );
            };
        };

        #ifdef SECDEC_WITH_CUDA
            template<>
            struct WithCuda<true>
            {
                using integrand_t = ::%(sub_integral_name)s::cuda_integrand_t;

                // Note: we define make_integrands with %(name)s_contour_deformation
                // but call ::sub_integral_name::make_integrands with %(sub_integral_name)s_contour_deformation
                // i.e we drop contour deformation parameters not relevant for this integral
                static std::vector<nested_series_t<integrand_t>> make_integrands
                (
                    const std::vector<real_t>& real_parameters,
                    const std::vector<complex_t>& complex_parameters
                    #if %(name)s_contour_deformation
                        ,unsigned number_of_presamples,
                        real_t deformation_parameters_maximum,
                        real_t deformation_parameters_minimum,
                        real_t deformation_parameters_decrease_factor
                    #endif
                )
                {
                        return ::%(sub_integral_name)s::make_cuda_integrands
                        (
                            real_parameters,
                            complex_parameters
                            #if %(sub_integral_name)s_contour_deformation
                                ,number_of_presamples,
                                deformation_parameters_maximum,
                                deformation_parameters_minimum,
                                deformation_parameters_decrease_factor
                            #endif
                        );
                };
            };
        #endif
        
        // secdecutil::cuba::Vegas, secdecutil::cuba::Suave, secdecutil::cuba::Cuhre, secdecutil::cuba::Divonne
        template<typename integrand_return_t, typename real_t, typename integrator_t, typename integrand_t>
        struct AmplitudeIntegral
        {
            using amplitude_integrator_t = integrator_t;
            using amplitude_integral_t = secdecutil::amplitude::CubaIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>;
        };
        
        // secdecutil::gsl::CQuad
        template<typename integrand_return_t, typename real_t, typename integrand_t>
        struct AmplitudeIntegral<integrand_return_t, real_t, secdecutil::gsl::CQuad<INTEGRAL_NAME::integrand_return_t>, integrand_t>
        {
            using amplitude_integrator_t = secdecutil::gsl::CQuad<INTEGRAL_NAME::integrand_return_t>;
            using amplitude_integral_t = secdecutil::amplitude::CQuadIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>;
        };
        
        // secdecutil::MultiIntegrator
        template<typename integrand_return_t, typename real_t, typename integrand_t>
        struct AmplitudeIntegral<integrand_return_t, real_t, multiintegrator_t, integrand_t>
        {
            using amplitude_integrator_t = multiintegrator_t;
            using amplitude_integral_t = secdecutil::amplitude::MultiIntegratorIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>;
        };

        // Note: for the CUDA-compatible integrators we take an integrator with INTEGRAL_NAME::(cuda_)integrand_t
        // and return one with ::%(sub_integral_name)s::(cuda_)integrand_t instead
        #define INSTANTIATE_AMPLITUDE_INTEGRAL_NONE_QMC() \
            template<typename integrand_return_t, typename real_t, typename integrand_t> \
            struct AmplitudeIntegral<integrand_return_t, real_t, secdecutil::integrators::Qmc< \
                 INTEGRAL_NAME::integrand_return_t, \
                 INTEGRAL_NAME::maximal_number_of_integration_variables, \
                 ::integrators::transforms::None::type, \
                 INTEGRAL_NAME::INTEGRAND_TYPE, \
                 secdecutil::integrators::void_template \
            > \
            , integrand_t> \
            { \
                using amplitude_integrator_t = secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::None::type, \
                     ::%(sub_integral_name)s::INTEGRAND_TYPE, \
                     secdecutil::integrators::void_template \
                >; \
                using amplitude_integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>; \
            }; \
            template<typename integrand_return_t, typename real_t, typename integrand_t> \
            struct AmplitudeIntegral<integrand_return_t, real_t, secdecutil::integrators::Qmc< \
                 INTEGRAL_NAME::integrand_return_t, \
                 INTEGRAL_NAME::maximal_number_of_integration_variables, \
                 ::integrators::transforms::None::type, \
                 INTEGRAL_NAME::INTEGRAND_TYPE, \
                 ::integrators::fitfunctions::None::type \
            > \
            , integrand_t> \
            { \
                using amplitude_integrator_t = secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::None::type, \
                     ::%(sub_integral_name)s::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::None::type \
                >; \
                using amplitude_integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>; \
            }; \
            template<typename integrand_return_t, typename real_t, typename integrand_t> \
            struct AmplitudeIntegral<integrand_return_t, real_t, secdecutil::integrators::Qmc< \
                 INTEGRAL_NAME::integrand_return_t, \
                 INTEGRAL_NAME::maximal_number_of_integration_variables, \
                 ::integrators::transforms::None::type, \
                 INTEGRAL_NAME::INTEGRAND_TYPE, \
                 ::integrators::fitfunctions::PolySingular::type \
            > \
            , integrand_t> \
            { \
                using amplitude_integrator_t = secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::None::type, \
                     ::%(sub_integral_name)s::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::PolySingular::type \
                >; \
                using amplitude_integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>; \
            };

        #define INSTANTIATE_AMPLITUDE_INTEGRAL_BAKER_QMC() \
            template<typename integrand_return_t, typename real_t, typename integrand_t> \
            struct AmplitudeIntegral<integrand_return_t, real_t, secdecutil::integrators::Qmc< \
                 INTEGRAL_NAME::integrand_return_t, \
                 INTEGRAL_NAME::maximal_number_of_integration_variables, \
                 ::integrators::transforms::Baker::type, \
                 INTEGRAL_NAME::INTEGRAND_TYPE, \
                 secdecutil::integrators::void_template \
            > \
            , integrand_t> \
            { \
                using amplitude_integrator_t = secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Baker::type, \
                     ::%(sub_integral_name)s::INTEGRAND_TYPE, \
                     secdecutil::integrators::void_template \
                >; \
                using amplitude_integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>; \
            }; \
            template<typename integrand_return_t, typename real_t, typename integrand_t> \
            struct AmplitudeIntegral<integrand_return_t, real_t, secdecutil::integrators::Qmc< \
                 INTEGRAL_NAME::integrand_return_t, \
                 INTEGRAL_NAME::maximal_number_of_integration_variables, \
                 ::integrators::transforms::Baker::type, \
                 INTEGRAL_NAME::INTEGRAND_TYPE, \
                 ::integrators::fitfunctions::None::type \
            > \
            , integrand_t> \
            { \
                using amplitude_integrator_t = secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Baker::type, \
                     ::%(sub_integral_name)s::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::None::type \
                >; \
                using amplitude_integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>; \
            }; \
            template<typename integrand_return_t, typename real_t, typename integrand_t> \
            struct AmplitudeIntegral<integrand_return_t, real_t, secdecutil::integrators::Qmc< \
                 INTEGRAL_NAME::integrand_return_t, \
                 INTEGRAL_NAME::maximal_number_of_integration_variables, \
                 ::integrators::transforms::Baker::type, \
                 INTEGRAL_NAME::INTEGRAND_TYPE, \
                 ::integrators::fitfunctions::PolySingular::type \
            > \
            , integrand_t> \
            { \
                using amplitude_integrator_t = secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Baker::type, \
                     ::%(sub_integral_name)s::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::PolySingular::type \
                >; \
                using amplitude_integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>; \
            };
        
        #define INSTANTIATE_AMPLITUDE_INTEGRAL_KOROBOV_QMC(KOROBOVDEGREE1,KOROBOVDEGREE2) \
            template<typename integrand_return_t, typename real_t, typename integrand_t> \
            struct AmplitudeIntegral<integrand_return_t, real_t, secdecutil::integrators::Qmc< \
                 INTEGRAL_NAME::integrand_return_t, \
                 INTEGRAL_NAME::maximal_number_of_integration_variables, \
                 ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                 INTEGRAL_NAME::INTEGRAND_TYPE, \
                 secdecutil::integrators::void_template \
            > \
            , integrand_t> \
            { \
                using amplitude_integrator_t = secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                     ::%(sub_integral_name)s::INTEGRAND_TYPE, \
                     secdecutil::integrators::void_template \
                >; \
                using amplitude_integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>; \
            }; \
            template<typename integrand_return_t, typename real_t, typename integrand_t> \
            struct AmplitudeIntegral<integrand_return_t, real_t, secdecutil::integrators::Qmc< \
                 INTEGRAL_NAME::integrand_return_t, \
                 INTEGRAL_NAME::maximal_number_of_integration_variables, \
                 ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                 INTEGRAL_NAME::INTEGRAND_TYPE, \
                 ::integrators::fitfunctions::None::type \
            > \
            , integrand_t> \
            { \
                using amplitude_integrator_t = secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                     ::%(sub_integral_name)s::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::None::type \
                >; \
                using amplitude_integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>; \
            }; \
            template<typename integrand_return_t, typename real_t, typename integrand_t> \
            struct AmplitudeIntegral<integrand_return_t, real_t, secdecutil::integrators::Qmc< \
                 INTEGRAL_NAME::integrand_return_t, \
                 INTEGRAL_NAME::maximal_number_of_integration_variables, \
                 ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                 INTEGRAL_NAME::INTEGRAND_TYPE, \
                 ::integrators::fitfunctions::PolySingular::type \
            > \
            , integrand_t> \
            { \
                using amplitude_integrator_t = secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                     ::%(sub_integral_name)s::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::PolySingular::type \
                >; \
                using amplitude_integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>; \
            };
        
        #define INSTANTIATE_AMPLITUDE_INTEGRAL_SIDI_QMC(SIDIDEGREE) \
            template<typename integrand_return_t, typename real_t, typename integrand_t> \
            struct AmplitudeIntegral<integrand_return_t, real_t, secdecutil::integrators::Qmc< \
                 INTEGRAL_NAME::integrand_return_t, \
                 INTEGRAL_NAME::maximal_number_of_integration_variables, \
                 ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                 INTEGRAL_NAME::INTEGRAND_TYPE, \
                 secdecutil::integrators::void_template \
            > \
            , integrand_t> \
            { \
                using amplitude_integrator_t = secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                     ::%(sub_integral_name)s::INTEGRAND_TYPE, \
                     secdecutil::integrators::void_template \
                >; \
                using amplitude_integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>; \
            }; \
            template<typename integrand_return_t, typename real_t, typename integrand_t> \
            struct AmplitudeIntegral<integrand_return_t, real_t, secdecutil::integrators::Qmc< \
                 INTEGRAL_NAME::integrand_return_t, \
                 INTEGRAL_NAME::maximal_number_of_integration_variables, \
                 ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                 INTEGRAL_NAME::INTEGRAND_TYPE, \
                 ::integrators::fitfunctions::None::type \
            > \
            , integrand_t> \
            { \
                using amplitude_integrator_t = secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                     ::%(sub_integral_name)s::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::None::type \
                >; \
                using amplitude_integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>; \
            }; \
            template<typename integrand_return_t, typename real_t, typename integrand_t> \
            struct AmplitudeIntegral<integrand_return_t, real_t, secdecutil::integrators::Qmc< \
                 INTEGRAL_NAME::integrand_return_t, \
                 INTEGRAL_NAME::maximal_number_of_integration_variables, \
                 ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                 INTEGRAL_NAME::INTEGRAND_TYPE, \
                 ::integrators::fitfunctions::PolySingular::type \
            > \
            , integrand_t> \
            { \
                using amplitude_integrator_t = secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                     ::%(sub_integral_name)s::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::PolySingular::type \
                >; \
                using amplitude_integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,amplitude_integrator_t,integrand_t>; \
            };
        
        // secdecutil::integrators::Qmc
        %(pylink_qmc_instantiate_amplitude_integral)s
        
        // Note: we define make_integrands with %(name)s_contour_deformation
        // but call ::sub_integral_name::make_integrands with %(sub_integral_name)s_contour_deformation
        // i.e we drop contour deformation parameters not relevant for this integral
        template<typename integrator_t>
        std::vector<nested_series_t<sum_t>> make_integral
        (
            const std::vector<real_t>& real_parameters,
            const std::vector<complex_t>& complex_parameters,
            const integrator_t& integrator
            #if %(name)s_contour_deformation
                ,unsigned number_of_presamples,
                real_t deformation_parameters_maximum,
                real_t deformation_parameters_minimum,
                real_t deformation_parameters_decrease_factor
            #endif
        )
        {
            using types = WithCuda<integrator_t::cuda_compliant_integrator>;
            using integrand_t = typename types::integrand_t;
            using integrand_return_t = %(name)s::integrand_return_t;
            using amplitude_integrator_t = typename AmplitudeIntegral<integrand_return_t, real_t, integrator_t, integrand_t>::amplitude_integrator_t;
            using amplitude_integral_t = typename AmplitudeIntegral<integrand_return_t, real_t, integrator_t, integrand_t>::amplitude_integral_t;

            const std::vector<nested_series_t<integrand_t>> raw_integrands =
            types::make_integrands
            (
                real_parameters,
                complex_parameters
                #if %(name)s_contour_deformation
                    ,number_of_presamples,
                    deformation_parameters_maximum,
                    deformation_parameters_minimum,
                    deformation_parameters_decrease_factor
                #endif
            );

            // Instantiate an amplitude_integrator_t from integrator, store this instance in a shared pointer
            const std::shared_ptr<amplitude_integrator_t> integrator_ptr = std::make_shared<amplitude_integrator_t>(integrator);

            const std::function<sum_t(const integrand_t& integrand)> convert_integrands =
                [ integrator_ptr ] (const integrand_t& integrand) -> sum_t
                {
                    const std::shared_ptr<amplitude_integral_t> integral_ptr = std::make_shared<amplitude_integral_t>(integrator_ptr, integrand);
                    integral_ptr->display_name = ::%(sub_integral_name)s::package_name + "_" + integrand.display_name;
                    return { /* constructor of std::vector */
                                { /* constructor of WeightedIntegral */
                                    integral_ptr
                                }
                        };
                };

            return deep_apply(raw_integrands, convert_integrands);
        };

        nested_series_t<sum_t> make_weighted_integral
        (
            const std::vector<real_t>& real_parameters,
            const std::vector<complex_t>& complex_parameters,
            const std::vector<nested_series_t<sum_t>>& integrals,
            const unsigned int amp_idx,
            const std::string& lib_path
        )
        {
            nested_series_t<sum_t> amplitude = std::accumulate(++integrals.begin(), integrals.end(), *integrals.begin() );
            amplitude *= ::%(sub_integral_name)s::prefactor(real_parameters,complex_parameters)
                       * ::%(sub_integral_name)s::coefficient(real_parameters,complex_parameters,amp_idx,lib_path);
            return amplitude;
        }
        
        #if %(name)s_contour_deformation
        
            #define INSTANTIATE_MAKE_INTEGRAL(INTEGRATOR) \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                        const INTEGRATOR&, unsigned, real_t, real_t, real_t);

            #define INSTANTIATE_MAKE_INTEGRAL_NONE_QMC() \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::None::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         secdecutil::integrators::void_template \
                    >&, \
                    unsigned, real_t, real_t, real_t); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::None::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         ::integrators::fitfunctions::None::type \
                    >&, \
                    unsigned, real_t, real_t, real_t); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::None::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&, \
                    unsigned, real_t, real_t, real_t);
        
            #define INSTANTIATE_MAKE_INTEGRAL_BAKER_QMC() \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::Baker::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         secdecutil::integrators::void_template \
                    >&, \
                    unsigned, real_t, real_t, real_t); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::Baker::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         ::integrators::fitfunctions::None::type \
                    >&, \
                    unsigned, real_t, real_t, real_t); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Baker::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&, \
                    unsigned, real_t, real_t, real_t);
        
            #define INSTANTIATE_MAKE_INTEGRAL_KOROBOV_QMC(KOROBOVDEGREE1,KOROBOVDEGREE2) \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         secdecutil::integrators::void_template \
                    >&, \
                    unsigned, real_t, real_t, real_t); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         ::integrators::fitfunctions::None::type \
                    >&, \
                    unsigned, real_t, real_t, real_t); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&, \
                    unsigned, real_t, real_t, real_t);
        
            #define INSTANTIATE_MAKE_INTEGRAL_SIDI_QMC(SIDIDEGREE) \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         secdecutil::integrators::void_template \
                    >&, \
                    unsigned, real_t, real_t, real_t); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         ::integrators::fitfunctions::None::type \
                    >&, \
                    unsigned, real_t, real_t, real_t); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&, \
                    unsigned, real_t, real_t, real_t);
        
        #else
        
            #define INSTANTIATE_MAKE_INTEGRAL(INTEGRATOR) \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, const INTEGRATOR&);
  
            #define INSTANTIATE_MAKE_INTEGRAL_NONE_QMC() \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::None::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         secdecutil::integrators::void_template \
                    >&); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::None::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         ::integrators::fitfunctions::None::type \
                    >&); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::None::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&);
        
            #define INSTANTIATE_MAKE_INTEGRAL_BAKER_QMC() \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::Baker::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         secdecutil::integrators::void_template \
                    >&); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::Baker::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         ::integrators::fitfunctions::None::type \
                    >&); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Baker::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&);
        
            #define INSTANTIATE_MAKE_INTEGRAL_KOROBOV_QMC(KOROBOVDEGREE1,KOROBOVDEGREE2) \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         secdecutil::integrators::void_template \
                    >&); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         ::integrators::fitfunctions::None::type \
                    >&); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&);
        
            #define INSTANTIATE_MAKE_INTEGRAL_SIDI_QMC(SIDIDEGREE) \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         secdecutil::integrators::void_template \
                    >&); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                         INTEGRAL_NAME::integrand_return_t, \
                         INTEGRAL_NAME::maximal_number_of_integration_variables, \
                         ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                         INTEGRAL_NAME::INTEGRAND_TYPE, \
                         ::integrators::fitfunctions::None::type \
                    >&); \
                template std::vector<nested_series_t<sum_t>> make_integral(const std::vector<real_t>&, const std::vector<complex_t>&, \
                    const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&);
        #endif
        
        INSTANTIATE_MAKE_INTEGRAL(secdecutil::gsl::CQuad<INTEGRAL_NAME::integrand_return_t>)

        // secdecutil::cuba::Vegas, secdecutil::cuba::Suave, secdecutil::cuba::Cuhre, secdecutil::cuba::Divonne
        INSTANTIATE_MAKE_INTEGRAL(secdecutil::cuba::Vegas<INTEGRAL_NAME::integrand_return_t>)
        INSTANTIATE_MAKE_INTEGRAL(secdecutil::cuba::Suave<INTEGRAL_NAME::integrand_return_t>)
        INSTANTIATE_MAKE_INTEGRAL(secdecutil::cuba::Cuhre<INTEGRAL_NAME::integrand_return_t>)
        INSTANTIATE_MAKE_INTEGRAL(secdecutil::cuba::Divonne<INTEGRAL_NAME::integrand_return_t>)
        
        // secdecutil::MultiIntegrator
        INSTANTIATE_MAKE_INTEGRAL(multiintegrator_t)
        
        // secdecutil::integrators::Qmc
        %(pylink_qmc_instantiate_make_integral)s
        
        #undef INSTANTIATE_AMPLITUDE_INTEGRAL_NONE_QMC
        #undef INSTANTIATE_AMPLITUDE_INTEGRAL_BAKER_QMC
        #undef INSTANTIATE_AMPLITUDE_INTEGRAL_KOROBOV_QMC
        #undef INSTANTIATE_AMPLITUDE_INTEGRAL_SIDI_QMC
        #undef INSTANTIATE_MAKE_AMPLITUDES
        #undef INSTANTIATE_MAKE_INTEGRAL_NONE_QMC
        #undef INSTANTIATE_MAKE_INTEGRAL_BAKER_QMC
        #undef INSTANTIATE_MAKE_INTEGRAL_KOROBOV_QMC
        #undef INSTANTIATE_MAKE_INTEGRAL_SIDI_QMC
        
    };
};

#undef INTEGRAL_NAME
#undef INTEGRAND_TYPE

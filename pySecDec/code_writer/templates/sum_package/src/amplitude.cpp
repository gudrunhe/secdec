#include <vector> // std::vector
#include <stdexcept> // std::invalid_argument
#include <string> // std::string
#include <typeinfo> // typeid

#include <secdecutil/integrators/cquad.hpp> // secdecutil::gsl::CQuad
#include <secdecutil/integrators/cuba.hpp> // secdecutil::cuba::Vegas, secdecutil::cuba::Suave, secdecutil::cuba::Cuhre, secdecutil::cuba::Divonne
#include <secdecutil/integrators/qmc.hpp> // secdecutil::integrators::Qmc

#include "%(name)s.hpp"

%(weighted_integral_includes)s

#define EXPAND_STRINGIFY(item) STRINGIFY(item)
#define STRINGIFY(item) #item

#define INTEGRAL_NAME %(name)s
#ifdef SECDEC_WITH_CUDA
    #define INTEGRAND_TYPE cuda_integrand_t
#else
    #define INTEGRAND_TYPE integrand_t
#endif

namespace %(name)s
{
    typedef secdecutil::MultiIntegrator<INTEGRAL_NAME::integrand_return_t,INTEGRAL_NAME::real_t,INTEGRAL_NAME::integrand_t> multiintegrator_t;
    
    template<typename integrator_t>
    std::vector<nested_series_t<sum_t>> make_amplitudes
    (
        const std::vector<real_t>& real_parameters,
        const std::vector<complex_t>& complex_parameters,
        const std::string& lib_path,
        const integrator_t& integrator
        #if %(name)s_contour_deformation
            ,unsigned number_of_presamples,
            real_t deformation_parameters_maximum,
            real_t deformation_parameters_minimum,
            real_t deformation_parameters_decrease_factor
        #endif
    )
    {
        // Construct Integrals
        #if %(name)s_contour_deformation
            %(integral_initialization_with_contour_deformation)s
        #else
            %(integral_initialization)s
        #endif
        
        // Construct Amplitudes
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
    
    // Specialisation for pylink interface
    #if %(name)s_contour_deformation
    
        #define DYNAMIC_CAST_INTEGRATOR_NONE_QMC() \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::None::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                secdecutil::integrators::void_template \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::None::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        secdecutil::integrators::void_template \
                    >&>(*integrator) \
                    ,number_of_presamples, \
                    deformation_parameters_maximum, \
                    deformation_parameters_minimum, \
                    deformation_parameters_decrease_factor \
                ); \
            }\
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::None::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::None::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::None::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::None::type \
                    >&>(*integrator) \
                    ,number_of_presamples, \
                    deformation_parameters_maximum, \
                    deformation_parameters_minimum, \
                    deformation_parameters_decrease_factor \
                ); \
            } \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::None::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::PolySingular::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::None::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&>(*integrator) \
                    ,number_of_presamples, \
                    deformation_parameters_maximum, \
                    deformation_parameters_minimum, \
                    deformation_parameters_decrease_factor \
                ); \
            }
    
        #define DYNAMIC_CAST_INTEGRATOR_BAKER_QMC() \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Baker::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                secdecutil::integrators::void_template \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Baker::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        secdecutil::integrators::void_template \
                    >&>(*integrator) \
                    ,number_of_presamples, \
                    deformation_parameters_maximum, \
                    deformation_parameters_minimum, \
                    deformation_parameters_decrease_factor \
                ); \
            }\
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Baker::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::None::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Baker::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::None::type \
                    >&>(*integrator) \
                    ,number_of_presamples, \
                    deformation_parameters_maximum, \
                    deformation_parameters_minimum, \
                    deformation_parameters_decrease_factor \
                ); \
            } \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Baker::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::PolySingular::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Baker::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&>(*integrator) \
                    ,number_of_presamples, \
                    deformation_parameters_maximum, \
                    deformation_parameters_minimum, \
                    deformation_parameters_decrease_factor \
                ); \
            }
    
        #define DYNAMIC_CAST_INTEGRATOR_KOROBOV_QMC(KOROBOVDEGREE1,KOROBOVDEGREE2) \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                secdecutil::integrators::void_template \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        secdecutil::integrators::void_template \
                    >&>(*integrator) \
                    ,number_of_presamples, \
                    deformation_parameters_maximum, \
                    deformation_parameters_minimum, \
                    deformation_parameters_decrease_factor \
                ); \
            }\
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::None::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::None::type \
                    >&>(*integrator) \
                    ,number_of_presamples, \
                    deformation_parameters_maximum, \
                    deformation_parameters_minimum, \
                    deformation_parameters_decrease_factor \
                ); \
            } \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::PolySingular::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&>(*integrator) \
                    ,number_of_presamples, \
                    deformation_parameters_maximum, \
                    deformation_parameters_minimum, \
                    deformation_parameters_decrease_factor \
                ); \
            }
    
        #define DYNAMIC_CAST_INTEGRATOR_SIDI_QMC(SIDIDEGREE) \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                secdecutil::integrators::void_template \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        secdecutil::integrators::void_template \
                    >&>(*integrator) \
                    ,number_of_presamples, \
                    deformation_parameters_maximum, \
                    deformation_parameters_minimum, \
                    deformation_parameters_decrease_factor \
                ); \
            }\
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::None::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::None::type \
                    >&>(*integrator) \
                    ,number_of_presamples, \
                    deformation_parameters_maximum, \
                    deformation_parameters_minimum, \
                    deformation_parameters_decrease_factor \
                ); \
            } \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::PolySingular::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&>(*integrator) \
                    ,number_of_presamples, \
                    deformation_parameters_maximum, \
                    deformation_parameters_minimum, \
                    deformation_parameters_decrease_factor \
                ); \
            }
    
        #define DYNAMIC_CAST_INTEGRATOR(INTEGRATOR) \
            if(dynamic_cast<const INTEGRATOR*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const INTEGRATOR&>(*integrator) \
                    ,number_of_presamples, \
                    deformation_parameters_maximum, \
                    deformation_parameters_minimum, \
                    deformation_parameters_decrease_factor \
                ); \
            }
    
    #else
    
        #define DYNAMIC_CAST_INTEGRATOR_NONE_QMC() \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::None::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                secdecutil::integrators::void_template \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::None::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        secdecutil::integrators::void_template \
                    >&>(*integrator) \
                ); \
            }\
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::None::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::None::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::None::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::None::type \
                    >&>(*integrator) \
                ); \
            } \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::None::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::PolySingular::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::None::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&>(*integrator) \
                ); \
            }

        #define DYNAMIC_CAST_INTEGRATOR_BAKER_QMC() \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Baker::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                secdecutil::integrators::void_template \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Baker::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        secdecutil::integrators::void_template \
                    >&>(*integrator) \
                ); \
            }\
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Baker::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::None::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Baker::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::None::type \
                    >&>(*integrator) \
                ); \
            } \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Baker::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::PolySingular::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Baker::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&>(*integrator) \
                ); \
            }
    
        #define DYNAMIC_CAST_INTEGRATOR_KOROBOV_QMC(KOROBOVDEGREE1,KOROBOVDEGREE2) \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                secdecutil::integrators::void_template \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        secdecutil::integrators::void_template \
                    >&>(*integrator) \
                ); \
            }\
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::None::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::None::type \
                    >&>(*integrator) \
                ); \
            } \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::PolySingular::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&>(*integrator) \
                ); \
            }
    
        #define DYNAMIC_CAST_INTEGRATOR_SIDI_QMC(SIDIDEGREE) \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                secdecutil::integrators::void_template \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        secdecutil::integrators::void_template \
                    >&>(*integrator) \
                ); \
            }\
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::None::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::None::type \
                    >&>(*integrator) \
                ); \
            } \
            if(dynamic_cast<const secdecutil::integrators::Qmc< \
                INTEGRAL_NAME::integrand_return_t, \
                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                INTEGRAL_NAME::INTEGRAND_TYPE, \
                ::integrators::fitfunctions::PolySingular::type \
            >*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const secdecutil::integrators::Qmc< \
                        INTEGRAL_NAME::integrand_return_t, \
                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                        INTEGRAL_NAME::INTEGRAND_TYPE, \
                        ::integrators::fitfunctions::PolySingular::type \
                    >&>(*integrator) \
                ); \
            }
    
        #define DYNAMIC_CAST_INTEGRATOR(INTEGRATOR) \
            if(dynamic_cast<const INTEGRATOR*>(integrator)) \
            { \
                return make_amplitudes( \
                    real_parameters, \
                    complex_parameters, \
                    lib_path, \
                    dynamic_cast<const INTEGRATOR&>(*integrator) \
                ); \
            }
    
    #endif
    
    std::vector<nested_series_t<sum_t>> make_amplitudes
    (
        const std::vector<real_t>& real_parameters,
        const std::vector<complex_t>& complex_parameters,
        const std::string& lib_path,
        const secdecutil::Integrator<integrand_return_t,real_t> * integrator
        #if %(name)s_contour_deformation
            ,unsigned number_of_presamples,
            real_t deformation_parameters_maximum,
            real_t deformation_parameters_minimum,
            real_t deformation_parameters_decrease_factor
        #endif
    )
    {
        DYNAMIC_CAST_INTEGRATOR(secdecutil::gsl::CQuad<INTEGRAL_NAME::integrand_return_t>)
        
        // secdecutil::cuba::Vegas, secdecutil::cuba::Suave, secdecutil::cuba::Cuhre, secdecutil::cuba::Divonne
        DYNAMIC_CAST_INTEGRATOR(secdecutil::cuba::Vegas<INTEGRAL_NAME::integrand_return_t>)
        DYNAMIC_CAST_INTEGRATOR(secdecutil::cuba::Suave<INTEGRAL_NAME::integrand_return_t>)
        DYNAMIC_CAST_INTEGRATOR(secdecutil::cuba::Cuhre<INTEGRAL_NAME::integrand_return_t>)
        DYNAMIC_CAST_INTEGRATOR(secdecutil::cuba::Divonne<INTEGRAL_NAME::integrand_return_t>)
        
        // secdecutil::MultiIntegrator
        DYNAMIC_CAST_INTEGRATOR(multiintegrator_t)
        
        // secdecutil::integrators::Qmc
        %(pylink_qmc_dynamic_cast_integrator)s
        
        // None of the above dynamic_casts succeeded, throw and give up
        throw std::invalid_argument("Trying to call \"" EXPAND_STRINGIFY(INTEGRAL_NAME) "::make_amplitudes\" with unknown \"secdecutil::Integrator\" derived type: " + std::string(typeid(*integrator).name()));
    }

    #ifdef SECDEC_WITH_CUDA
        std::vector<nested_series_t<sum_t>> make_amplitudes
        (
            const std::vector<real_t>& real_parameters,
            const std::vector<complex_t>& complex_parameters,
            const std::string& lib_path,
            const secdecutil::Integrator<integrand_return_t,real_t,cuda_integrand_t> * integrator
            #if %(name)s_contour_deformation
                ,unsigned number_of_presamples,
                real_t deformation_parameters_maximum,
                real_t deformation_parameters_minimum,
                real_t deformation_parameters_decrease_factor
            #endif
        )
        {
            // secdecutil::integrators::Qmc
            %(pylink_qmc_dynamic_cast_integrator)s
            
            // None of the above dynamic_casts succeeded, throw and give up
            throw std::invalid_argument("Trying to call \"" EXPAND_STRINGIFY(INTEGRAL_NAME) "::make_amplitudes\" with unknown \"secdecutil::Integrator\" derived type: " + std::string(typeid(*integrator).name()));
        }
    #endif
    
    #if %(name)s_contour_deformation
    
        #define INSTANTIATE_MAKE_AMPLITUDES(INTEGRATOR) \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                    const INTEGRATOR&, unsigned, real_t, real_t, real_t);

        #define INSTANTIATE_MAKE_AMPLITUDES_NONE_QMC() \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::None::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     secdecutil::integrators::void_template \
                >&, \
                unsigned, real_t, real_t, real_t); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::None::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::None::type \
                >&, \
                unsigned, real_t, real_t, real_t); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                    INTEGRAL_NAME::integrand_return_t, \
                    INTEGRAL_NAME::maximal_number_of_integration_variables, \
                    ::integrators::transforms::None::type, \
                    INTEGRAL_NAME::INTEGRAND_TYPE, \
                    ::integrators::fitfunctions::PolySingular::type \
                >&, \
                unsigned, real_t, real_t, real_t);

        #define INSTANTIATE_MAKE_AMPLITUDES_BAKER_QMC() \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Baker::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     secdecutil::integrators::void_template \
                >&, \
                unsigned, real_t, real_t, real_t); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Baker::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::None::type \
                >&, \
                unsigned, real_t, real_t, real_t); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                    INTEGRAL_NAME::integrand_return_t, \
                    INTEGRAL_NAME::maximal_number_of_integration_variables, \
                    ::integrators::transforms::Baker::type, \
                    INTEGRAL_NAME::INTEGRAND_TYPE, \
                    ::integrators::fitfunctions::PolySingular::type \
                >&, \
                unsigned, real_t, real_t, real_t);
    
        #define INSTANTIATE_MAKE_AMPLITUDES_KOROBOV_QMC(KOROBOVDEGREE1,KOROBOVDEGREE2) \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     secdecutil::integrators::void_template \
                >&, \
                unsigned, real_t, real_t, real_t); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::None::type \
                >&, \
                unsigned, real_t, real_t, real_t); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                    INTEGRAL_NAME::integrand_return_t, \
                    INTEGRAL_NAME::maximal_number_of_integration_variables, \
                    ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                    INTEGRAL_NAME::INTEGRAND_TYPE, \
                    ::integrators::fitfunctions::PolySingular::type \
                >&, \
                unsigned, real_t, real_t, real_t);
    
        #define INSTANTIATE_MAKE_AMPLITUDES_SIDI_QMC(SIDIDEGREE) \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     secdecutil::integrators::void_template \
                >&, \
                unsigned, real_t, real_t, real_t); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::None::type \
                >&, \
                unsigned, real_t, real_t, real_t); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                    INTEGRAL_NAME::integrand_return_t, \
                    INTEGRAL_NAME::maximal_number_of_integration_variables, \
                    ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                    INTEGRAL_NAME::INTEGRAND_TYPE, \
                    ::integrators::fitfunctions::PolySingular::type \
                >&, \
                unsigned, real_t, real_t, real_t);
    
    #else
    
        #define INSTANTIATE_MAKE_AMPLITUDES(INTEGRATOR) \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, const INTEGRATOR&);
    
        #define INSTANTIATE_MAKE_AMPLITUDES_NONE_QMC() \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::None::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     secdecutil::integrators::void_template \
                >&); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::None::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::None::type \
                >&); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                    INTEGRAL_NAME::integrand_return_t, \
                    INTEGRAL_NAME::maximal_number_of_integration_variables, \
                    ::integrators::transforms::None::type, \
                    INTEGRAL_NAME::INTEGRAND_TYPE, \
                    ::integrators::fitfunctions::PolySingular::type \
                >&);
    
        #define INSTANTIATE_MAKE_AMPLITUDES_BAKER_QMC() \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Baker::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     secdecutil::integrators::void_template \
                >&); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Baker::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::None::type \
                >&); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                    INTEGRAL_NAME::integrand_return_t, \
                    INTEGRAL_NAME::maximal_number_of_integration_variables, \
                    ::integrators::transforms::Baker::type, \
                    INTEGRAL_NAME::INTEGRAND_TYPE, \
                    ::integrators::fitfunctions::PolySingular::type \
                >&);
    
        #define INSTANTIATE_MAKE_AMPLITUDES_KOROBOV_QMC(KOROBOVDEGREE1,KOROBOVDEGREE2) \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     secdecutil::integrators::void_template \
                >&); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::None::type \
                >&); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                    INTEGRAL_NAME::integrand_return_t, \
                    INTEGRAL_NAME::maximal_number_of_integration_variables, \
                    ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                    INTEGRAL_NAME::INTEGRAND_TYPE, \
                    ::integrators::fitfunctions::PolySingular::type \
                >&);
    
        #define INSTANTIATE_MAKE_AMPLITUDES_SIDI_QMC(SIDIDEGREE) \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     secdecutil::integrators::void_template \
                >&); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                     INTEGRAL_NAME::integrand_return_t, \
                     INTEGRAL_NAME::maximal_number_of_integration_variables, \
                     ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                     INTEGRAL_NAME::INTEGRAND_TYPE, \
                     ::integrators::fitfunctions::None::type \
                >&); \
            template std::vector<nested_series_t<sum_t>> make_amplitudes(const std::vector<real_t>&, const std::vector<complex_t>&, const std::string&, \
                const secdecutil::integrators::Qmc< \
                    INTEGRAL_NAME::integrand_return_t, \
                    INTEGRAL_NAME::maximal_number_of_integration_variables, \
                    ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                    INTEGRAL_NAME::INTEGRAND_TYPE, \
                    ::integrators::fitfunctions::PolySingular::type \
                >&);
    #endif
    
    INSTANTIATE_MAKE_AMPLITUDES(secdecutil::gsl::CQuad<INTEGRAL_NAME::integrand_return_t>)

    // secdecutil::cuba::Vegas, secdecutil::cuba::Suave, secdecutil::cuba::Cuhre, secdecutil::cuba::Divonne
    INSTANTIATE_MAKE_AMPLITUDES(secdecutil::cuba::Vegas<INTEGRAL_NAME::integrand_return_t>)
    INSTANTIATE_MAKE_AMPLITUDES(secdecutil::cuba::Suave<INTEGRAL_NAME::integrand_return_t>)
    INSTANTIATE_MAKE_AMPLITUDES(secdecutil::cuba::Cuhre<INTEGRAL_NAME::integrand_return_t>)
    INSTANTIATE_MAKE_AMPLITUDES(secdecutil::cuba::Divonne<INTEGRAL_NAME::integrand_return_t>)
    
    // secdecutil::MultiIntegrator
    INSTANTIATE_MAKE_AMPLITUDES(multiintegrator_t)
    
    // secdecutil::integrators::Qmc
    %(pylink_qmc_instantiate_make_amplitudes)s
    
    #undef INTEGRAL_NAME
    #undef INTEGRAND_TYPE
    #undef DYNAMIC_CAST_INTEGRATOR_NONE_QMC
    #undef DYNAMIC_CAST_INTEGRATOR_BAKER_QMC
    #undef DYNAMIC_CAST_INTEGRATOR_KOROBOV_QMC
    #undef DYNAMIC_CAST_INTEGRATOR_SIDI_QMC
    #undef DYNAMIC_CAST_INTEGRATOR
    #undef INSTANTIATE_MAKE_AMPLITUDES
    #undef INSTANTIATE_MAKE_AMPLITUDES_NONE_QMC
    #undef INSTANTIATE_MAKE_AMPLITUDES_BAKER_QMC
    #undef INSTANTIATE_MAKE_AMPLITUDES_KOROBOV_QMC
    #undef INSTANTIATE_MAKE_AMPLITUDES_SIDI_QMC
};

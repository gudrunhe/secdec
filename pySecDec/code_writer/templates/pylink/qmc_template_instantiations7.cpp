#include "%(name)s.hpp"
#include <secdecutil/integrators/qmc.hpp> // Qmc

#define INTEGRAL_NAME %(name)s
#define %(name)s_number_of_sectors %(number_of_sectors)i

#ifdef SECDEC_WITH_CUDA
    #define INSTANTIATE_SIDI_QMC_SEPARATE(SIDIDEGREE) \
        template class secdecutil::integrators::Qmc< \
                                                        INTEGRAL_NAME::integrand_return_t, \
                                                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                        INTEGRAL_NAME::cuda_integrand_t, \
                                                        secdecutil::integrators::void_template \
                                                   >; \
        template class secdecutil::integrators::Qmc< \
                                                        INTEGRAL_NAME::integrand_return_t, \
                                                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                        INTEGRAL_NAME::cuda_integrand_t, \
                                                        ::integrators::fitfunctions::None::type \
                                                   >; \
        template class secdecutil::integrators::Qmc< \
                                                        INTEGRAL_NAME::integrand_return_t, \
                                                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                        INTEGRAL_NAME::cuda_integrand_t, \
                                                        ::integrators::fitfunctions::PolySingular::type \
                                                   >;
    #if %(name)s_number_of_sectors != 1
        #define INSTANTIATE_SIDI_QMC(SIDIDEGREE) \
            INSTANTIATE_SIDI_QMC_SEPARATE(SIDIDEGREE) \
            template class secdecutil::integrators::Qmc< \
                                                            INTEGRAL_NAME::integrand_return_t, \
                                                            INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                            ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                            INTEGRAL_NAME::cuda_together_integrand_t, \
                                                            secdecutil::integrators::void_template \
                                                       >; \
            template class secdecutil::integrators::Qmc< \
                                                            INTEGRAL_NAME::integrand_return_t, \
                                                            INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                            ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                            INTEGRAL_NAME::cuda_together_integrand_t, \
                                                            ::integrators::fitfunctions::None::type \
                                                       >; \
            template class secdecutil::integrators::Qmc< \
                                                            INTEGRAL_NAME::integrand_return_t, \
                                                            INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                            ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                            INTEGRAL_NAME::cuda_together_integrand_t, \
                                                            ::integrators::fitfunctions::PolySingular::type \
                                                       >;
    #else
        #define INSTANTIATE_SIDI_QMC(SIDIDEGREE) INSTANTIATE_SIDI_QMC_SEPARATE(SIDIDEGREE)
    #endif
#else
    #define INSTANTIATE_SIDI_QMC(SIDIDEGREE) \
        template class secdecutil::integrators::Qmc< \
                                                        INTEGRAL_NAME::integrand_return_t, \
                                                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                        INTEGRAL_NAME::integrand_t, \
                                                        secdecutil::integrators::void_template \
                                                   >; \
        template class secdecutil::integrators::Qmc< \
                                                        INTEGRAL_NAME::integrand_return_t, \
                                                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                        INTEGRAL_NAME::integrand_t, \
                                                        ::integrators::fitfunctions::None::type \
                                                   >; \
        template class secdecutil::integrators::Qmc< \
                                                        INTEGRAL_NAME::integrand_return_t, \
                                                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                        ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                        INTEGRAL_NAME::integrand_t, \
                                                        ::integrators::fitfunctions::PolySingular::type \
                                                   >;
#endif
INSTANTIATE_SIDI_QMC(1) INSTANTIATE_SIDI_QMC(2) INSTANTIATE_SIDI_QMC(3) INSTANTIATE_SIDI_QMC(4) INSTANTIATE_SIDI_QMC(5) INSTANTIATE_SIDI_QMC(6)
#undef INSTANTIATE_SIDI_QMC_SEPARATE
#undef INSTANTIATE_SIDI_QMC
#undef %(name)s_number_of_sectors

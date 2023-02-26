#include "%(name)s.hpp"
#include <secdecutil/integrators/qmc.hpp> // Qmc

#define INTEGRAL_NAME %(name)s

#ifdef SECDEC_WITH_CUDA
    #define INSTANTIATE_NONE_QMC() \
        template class secdecutil::integrators::Qmc< \
                                                                INTEGRAL_NAME::integrand_return_t, \
                                                                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                                ::integrators::transforms::None::type, \
                                                                INTEGRAL_NAME::cuda_integrand_t \
                                                        >; \
        template class secdecutil::integrators::Qmc< \
                                                                INTEGRAL_NAME::integrand_return_t, \
                                                                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                                ::integrators::transforms::None::type, \
                                                                INTEGRAL_NAME::cuda_integrand_t, \
                                                                ::integrators::fitfunctions::None::type \
                                                        >; \
        template class secdecutil::integrators::Qmc< \
                                                                INTEGRAL_NAME::integrand_return_t, \
                                                                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                                ::integrators::transforms::None::type, \
                                                                INTEGRAL_NAME::cuda_integrand_t, \
                                                                ::integrators::fitfunctions::PolySingular::type \
                                                        >;
    #define INSTANTIATE_BAKER_QMC() \
        template class secdecutil::integrators::Qmc< \
                                                            INTEGRAL_NAME::integrand_return_t, \
                                                            INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                            ::integrators::transforms::Baker::type, \
                                                            INTEGRAL_NAME::cuda_integrand_t \
                                                        >; \
        template class secdecutil::integrators::Qmc< \
                                                            INTEGRAL_NAME::integrand_return_t, \
                                                            INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                            ::integrators::transforms::Baker::type, \
                                                            INTEGRAL_NAME::cuda_integrand_t, \
                                                            ::integrators::fitfunctions::None::type \
                                                        >; \
        template class secdecutil::integrators::Qmc< \
                                                            INTEGRAL_NAME::integrand_return_t, \
                                                            INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                            ::integrators::transforms::Baker::type, \
                                                            INTEGRAL_NAME::cuda_integrand_t, \
                                                            ::integrators::fitfunctions::PolySingular::type \
                                                        >;
    #define INSTANTIATE_KOROBOV_QMC(KOROBOVDEGREE1,KOROBOVDEGREE2) \
            template class secdecutil::integrators::Qmc< \
                                                            INTEGRAL_NAME::integrand_return_t, \
                                                            INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                            ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                            INTEGRAL_NAME::cuda_integrand_t, \
                                                            secdecutil::integrators::void_template \
                                                       >; \
            template class secdecutil::integrators::Qmc< \
                                                            INTEGRAL_NAME::integrand_return_t, \
                                                            INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                            ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                            INTEGRAL_NAME::cuda_integrand_t, \
                                                            ::integrators::fitfunctions::None::type \
                                                       >; \
            template class secdecutil::integrators::Qmc< \
                                                            INTEGRAL_NAME::integrand_return_t, \
                                                            INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                            ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                            INTEGRAL_NAME::cuda_integrand_t, \
                                                            ::integrators::fitfunctions::PolySingular::type \
                                                       >;
#else
    #define INSTANTIATE_NONE_QMC() \
        template class secdecutil::integrators::Qmc< \
                                                                INTEGRAL_NAME::integrand_return_t, \
                                                                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                                ::integrators::transforms::None::type, \
                                                                INTEGRAL_NAME::integrand_t \
                                                        >; \
        template class secdecutil::integrators::Qmc< \
                                                                INTEGRAL_NAME::integrand_return_t, \
                                                                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                                ::integrators::transforms::None::type, \
                                                                INTEGRAL_NAME::integrand_t, \
                                                                ::integrators::fitfunctions::None::type \
                                                        >; \
        template class secdecutil::integrators::Qmc< \
                                                                INTEGRAL_NAME::integrand_return_t, \
                                                                INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                                ::integrators::transforms::None::type, \
                                                                INTEGRAL_NAME::integrand_t, \
                                                                ::integrators::fitfunctions::PolySingular::type \
                                                        >;
    #define INSTANTIATE_BAKER_QMC() \
        template class secdecutil::integrators::Qmc< \
                                                        INTEGRAL_NAME::integrand_return_t, \
                                                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                        ::integrators::transforms::Baker::type, \
                                                        INTEGRAL_NAME::integrand_t \
                                                    >; \
        template class secdecutil::integrators::Qmc< \
                                                        INTEGRAL_NAME::integrand_return_t, \
                                                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                        ::integrators::transforms::Baker::type, \
                                                        INTEGRAL_NAME::integrand_t, \
                                                        ::integrators::fitfunctions::None::type \
                                                    >; \
        template class secdecutil::integrators::Qmc< \
                                                        INTEGRAL_NAME::integrand_return_t, \
                                                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                        ::integrators::transforms::Baker::type, \
                                                        INTEGRAL_NAME::integrand_t, \
                                                        ::integrators::fitfunctions::PolySingular::type \
                                                    >;
    #define INSTANTIATE_KOROBOV_QMC(KOROBOVDEGREE1,KOROBOVDEGREE2) \
        template class secdecutil::integrators::Qmc< \
                                                        INTEGRAL_NAME::integrand_return_t, \
                                                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                        ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                        INTEGRAL_NAME::integrand_t, \
                                                        secdecutil::integrators::void_template \
                                                   >; \
        template class secdecutil::integrators::Qmc< \
                                                        INTEGRAL_NAME::integrand_return_t, \
                                                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                        ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                        INTEGRAL_NAME::integrand_t, \
                                                        ::integrators::fitfunctions::None::type \
                                                   >; \
        template class secdecutil::integrators::Qmc< \
                                                        INTEGRAL_NAME::integrand_return_t, \
                                                        INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                        ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                        INTEGRAL_NAME::integrand_t, \
                                                        ::integrators::fitfunctions::PolySingular::type \
                                                    >;
#endif

#ifdef SECDEC_WITH_CUDA
    #define INSTANTIATE_SIDI_QMC(SIDIDEGREE) \
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

%(pylink_qmc_transforms)s

#undef INSTANTIATE_NONE_QMC_SEPARATE
#undef INSTANTIATE_BAKER_QMC_SEPARATE
#undef INSTANTIATE_KOROBOV_QMC_SEPARATE
#undef INSTANTIATE_NONE_QMC
#undef INSTANTIATE_BAKER_QMC
#undef INSTANTIATE_KOROBOV_QMC
#undef INSTANTIATE_SIDI_QMC_SEPARATE
#undef INSTANTIATE_SIDI_QMC

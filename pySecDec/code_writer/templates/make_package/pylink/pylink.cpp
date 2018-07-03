#include "%(name)s.hpp"
#include <secdecutil/integrators/qmc.hpp> // Qmc



#define INTEGRAL_NAME %(name)s
#define %(name)s_number_of_sectors %(number_of_sectors)i

// whether or not to use contour deformation
#define integral_contour_deformation %(contour_deformation)i

// whether or not complex parameters are present
#define integral_has_complex_parameters %(have_complex_parameters)i

// whether or no the return type should be complex in any case
#define integral_enforce_complex_return_type %(enforce_complex_return_type)i

// delegate some template instatiations to separate translation units
#ifdef SECDEC_WITH_CUDA
    #define EXTERN_KOROBOV_QMC_SEPARATE(KOROBOVDEGREE1,KOROBOVDEGREE2) \
        extern template class secdecutil::integrators::Qmc< \
                                                               INTEGRAL_NAME::integrand_return_t, \
                                                               INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                               ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                               INTEGRAL_NAME::cuda_integrand_t, \
                                                               secdecutil::integrators::void_template \
                                                          >; \
        extern template class secdecutil::integrators::Qmc< \
                                                               INTEGRAL_NAME::integrand_return_t, \
                                                               INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                               ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                               INTEGRAL_NAME::cuda_integrand_t, \
                                                               ::integrators::fitfunctions::None::type \
                                                          >; \
        extern template class secdecutil::integrators::Qmc< \
                                                               INTEGRAL_NAME::integrand_return_t, \
                                                               INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                               ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                               INTEGRAL_NAME::cuda_integrand_t, \
                                                               ::integrators::fitfunctions::PolySingular::type \
                                                          >;
    #define EXTERN_SIDI_QMC_SEPARATE(SIDIDEGREE) \
        extern template class secdecutil::integrators::Qmc< \
                                                               INTEGRAL_NAME::integrand_return_t, \
                                                               INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                               ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                               INTEGRAL_NAME::cuda_integrand_t, \
                                                               secdecutil::integrators::void_template \
                                                          >; \
        extern template class secdecutil::integrators::Qmc< \
                                                               INTEGRAL_NAME::integrand_return_t, \
                                                               INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                               ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                               INTEGRAL_NAME::cuda_integrand_t, \
                                                               ::integrators::fitfunctions::None::type \
                                                          >; \
        extern template class secdecutil::integrators::Qmc< \
                                                               INTEGRAL_NAME::integrand_return_t, \
                                                               INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                               ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                               INTEGRAL_NAME::cuda_integrand_t, \
                                                               ::integrators::fitfunctions::PolySingular::type \
                                                          >;
    #if %(name)s_number_of_sectors != 1
        #define EXTERN_KOROBOV_QMC(KOROBOVDEGREE1,KOROBOVDEGREE2) \
            EXTERN_KOROBOV_QMC_SEPARATE(KOROBOVDEGREE1,KOROBOVDEGREE2) \
            extern template class secdecutil::integrators::Qmc< \
                                                                   INTEGRAL_NAME::integrand_return_t, \
                                                                   INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                                   ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                                   INTEGRAL_NAME::cuda_together_integrand_t, \
                                                                   secdecutil::integrators::void_template \
                                                              >; \
            extern template class secdecutil::integrators::Qmc< \
                                                                   INTEGRAL_NAME::integrand_return_t, \
                                                                   INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                                   ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                                   INTEGRAL_NAME::cuda_together_integrand_t, \
                                                                   ::integrators::fitfunctions::None::type \
                                                              >; \
            extern template class secdecutil::integrators::Qmc< \
                                                                   INTEGRAL_NAME::integrand_return_t, \
                                                                   INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                                   ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                                   INTEGRAL_NAME::cuda_together_integrand_t, \
                                                                   ::integrators::fitfunctions::PolySingular::type \
                                                              >;
        #define EXTERN_SIDI_QMC(SIDIDEGREE) \
            EXTERN_SIDI_QMC_SEPARATE(SIDIDEGREE) \
            extern template class secdecutil::integrators::Qmc< \
                                                                   INTEGRAL_NAME::integrand_return_t, \
                                                                   INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                                   ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                                   INTEGRAL_NAME::cuda_together_integrand_t, \
                                                                   secdecutil::integrators::void_template \
                                                              >; \
            extern template class secdecutil::integrators::Qmc< \
                                                                   INTEGRAL_NAME::integrand_return_t, \
                                                                   INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                                   ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                                   INTEGRAL_NAME::cuda_together_integrand_t, \
                                                                   ::integrators::fitfunctions::None::type \
                                                              >; \
            extern template class secdecutil::integrators::Qmc< \
                                                                   INTEGRAL_NAME::integrand_return_t, \
                                                                   INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                                   ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                                   INTEGRAL_NAME::cuda_together_integrand_t, \
                                                                   ::integrators::fitfunctions::PolySingular::type \
                                                              >;
    #else
        #define EXTERN_KOROBOV_QMC(KOROBOVDEGREE1,KOROBOVDEGREE2) EXTERN_KOROBOV_QMC_SEPARATE(KOROBOVDEGREE1,KOROBOVDEGREE2)
        #define EXTERN_SIDI_QMC(SIDIDEGREE) EXTERN_SIDI_QMC_SEPARATE(SIDIDEGREE)
    #endif
#else
    #define EXTERN_KOROBOV_QMC(KOROBOVDEGREE1,KOROBOVDEGREE2) \
        extern template class secdecutil::integrators::Qmc< \
                                                               INTEGRAL_NAME::integrand_return_t, \
                                                               INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                               ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                               INTEGRAL_NAME::integrand_t, \
                                                               secdecutil::integrators::void_template \
                                                          >; \
        extern template class secdecutil::integrators::Qmc< \
                                                               INTEGRAL_NAME::integrand_return_t, \
                                                               INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                               ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                               INTEGRAL_NAME::integrand_t, \
                                                               ::integrators::fitfunctions::None::type \
                                                          >; \
        extern template class secdecutil::integrators::Qmc< \
                                                               INTEGRAL_NAME::integrand_return_t, \
                                                               INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                               ::integrators::transforms::Korobov<KOROBOVDEGREE1,KOROBOVDEGREE2>::type, \
                                                               INTEGRAL_NAME::integrand_t, \
                                                               ::integrators::fitfunctions::PolySingular::type \
                                                          >;
    #define EXTERN_SIDI_QMC(SIDIDEGREE) \
        extern template class secdecutil::integrators::Qmc< \
                                                               INTEGRAL_NAME::integrand_return_t, \
                                                               INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                               ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                               INTEGRAL_NAME::integrand_t, \
                                                               secdecutil::integrators::void_template \
                                                          >; \
        extern template class secdecutil::integrators::Qmc< \
                                                               INTEGRAL_NAME::integrand_return_t, \
                                                               INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                               ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                               INTEGRAL_NAME::integrand_t, \
                                                               ::integrators::fitfunctions::None::type \
                                                          >; \
        extern template class secdecutil::integrators::Qmc< \
                                                               INTEGRAL_NAME::integrand_return_t, \
                                                               INTEGRAL_NAME::maximal_number_of_integration_variables, \
                                                               ::integrators::transforms::Sidi<SIDIDEGREE>::type, \
                                                               INTEGRAL_NAME::integrand_t, \
                                                               ::integrators::fitfunctions::PolySingular::type \
                                                          >;
#endif
EXTERN_KOROBOV_QMC(1,1) EXTERN_KOROBOV_QMC(1,2) EXTERN_KOROBOV_QMC(1,3) EXTERN_KOROBOV_QMC(1,4) EXTERN_KOROBOV_QMC(1,5) EXTERN_KOROBOV_QMC(1,6)
EXTERN_KOROBOV_QMC(2,1) EXTERN_KOROBOV_QMC(2,2) EXTERN_KOROBOV_QMC(2,3) EXTERN_KOROBOV_QMC(2,4) EXTERN_KOROBOV_QMC(2,5) EXTERN_KOROBOV_QMC(2,6)
EXTERN_KOROBOV_QMC(3,1) EXTERN_KOROBOV_QMC(3,2) EXTERN_KOROBOV_QMC(3,3) EXTERN_KOROBOV_QMC(3,4) EXTERN_KOROBOV_QMC(3,5) EXTERN_KOROBOV_QMC(3,6)
EXTERN_KOROBOV_QMC(4,1) EXTERN_KOROBOV_QMC(4,2) EXTERN_KOROBOV_QMC(4,3) EXTERN_KOROBOV_QMC(4,4) EXTERN_KOROBOV_QMC(4,5) EXTERN_KOROBOV_QMC(4,6)
EXTERN_KOROBOV_QMC(5,1) EXTERN_KOROBOV_QMC(5,2) EXTERN_KOROBOV_QMC(5,3) EXTERN_KOROBOV_QMC(5,4) EXTERN_KOROBOV_QMC(5,5) EXTERN_KOROBOV_QMC(5,6)
EXTERN_KOROBOV_QMC(6,1) EXTERN_KOROBOV_QMC(6,2) EXTERN_KOROBOV_QMC(6,3) EXTERN_KOROBOV_QMC(6,4) EXTERN_KOROBOV_QMC(6,5) EXTERN_KOROBOV_QMC(6,6)
EXTERN_SIDI_QMC(1) EXTERN_SIDI_QMC(2) EXTERN_SIDI_QMC(3) EXTERN_SIDI_QMC(4) EXTERN_SIDI_QMC(5) EXTERN_SIDI_QMC(6)
#undef EXTERN_KOROBOV_QMC
#undef EXTERN_SIDI_QMC
#undef EXTERN_KOROBOV_QMC_SEPARATE
#undef EXTERN_SIDI_QMC_SEPARATE
#undef %(name)s_number_of_sectors



// The python-C binding is general and therefore contained in the util
#include <secdecutil/pylink.hpp>



#undef integral_contour_deformation
#undef integral_has_complex_parameters
#undef integral_enforce_complex_return_type

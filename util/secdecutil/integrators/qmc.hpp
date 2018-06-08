#ifndef SecDecUtil_qmc_hpp_included
#define SecDecUtil_qmc_hpp_included

#ifdef SECDEC_WITH_CUDA
    #include <thrust/complex.h>
#endif
#include <complex>
#include <memory>
#include <stdexcept>
#include <qmc.hpp>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/uncertainties.hpp>
#include <secdecutil/integrators/integrator.hpp>

namespace secdecutil
{
    namespace integrators {

        template<typename return_t, typename container_t = secdecutil::IntegrandContainer<return_t, return_t const * const>, typename transform_t = void>
        struct Qmc : Integrator<return_t,return_t,container_t>, public ::integrators::Qmc<return_t,return_t>
        {
        protected:

            std::function<secdecutil::UncorrelatedDeviation<return_t>(const container_t&)>
            get_integrate()
            {
                std::function<secdecutil::UncorrelatedDeviation<return_t>(const container_t&)> integrate_function = [this] (const container_t& integrand_container)
                {
                    transform_t transform;
                    ::integrators::result<return_t> result = ::integrators::Qmc<return_t,return_t>::integrate
                    (
                        integrand_container,
                        integrand_container.number_of_integration_variables>0 ? integrand_container.number_of_integration_variables : 1, // ensure dim > 0
                        transform
                    );
                    return secdecutil::UncorrelatedDeviation<return_t> { result.integral, result.error };
                };
                return integrate_function;
            };

        public:

            using Integrator<return_t,return_t,container_t>::integrate;

        };

        template<typename return_t, typename container_t>
        struct Qmc<return_t,container_t,void> : Integrator<return_t,return_t,container_t>, public ::integrators::Qmc<return_t,return_t>
        {
        protected:

            std::function<secdecutil::UncorrelatedDeviation<return_t>(const container_t&)>
            get_integrate()
            {
                std::function<secdecutil::UncorrelatedDeviation<return_t>(const container_t&)> integrate_function = [this] (const container_t& integrand_container)
                {
                    ::integrators::result<return_t> result = ::integrators::Qmc<return_t,return_t>::integrate
                    (
                        integrand_container,
                        integrand_container.number_of_integration_variables>0 ? integrand_container.number_of_integration_variables : 1 // ensure dim > 0
                    );
                    return secdecutil::UncorrelatedDeviation<return_t> { result.integral, result.error };
                };
                return integrate_function;
            };

        public:

            using Integrator<return_t,return_t,container_t>::integrate;

        };

    #define COMPLEX_QMC_BODY_WITH_TRANSFORM(complex_template) \
        protected: \
 \
            std::unique_ptr<Integrator<return_t, return_t>> get_real_integrator() \
            { \
                throw std::runtime_error("Separate integration of real and imaginary part has no advantage for this non-adaptive integrator. Please use \"together = true\"."); \
            }; \
 \
            std::function<secdecutil::UncorrelatedDeviation<complex_template<return_t>> \
            (const container_t&)> \
            get_together_integrate() \
            { \
                std::function<secdecutil::UncorrelatedDeviation<complex_template<return_t>>(const container_t&)> integrate_function = [this] (const container_t& integrand_container) \
                { \
                    transform_t transform; \
                    ::integrators::result<complex_template<return_t>> result = ::integrators::Qmc<complex_template<return_t>,return_t>::integrate \
                    ( \
                        integrand_container, \
                        integrand_container.number_of_integration_variables>0 ? integrand_container.number_of_integration_variables : 1, /* ensure dim > 0 */ \
                        transform \
                    ); \
                    return secdecutil::UncorrelatedDeviation<complex_template<return_t>> { result.integral, result.error }; \
                }; \
                return integrate_function; \
            }; \
 \
        public: \
 \
            using Integrator<complex_template<return_t>, return_t, container_t>::integrate; \
 \
            Qmc() { this->together = true;  };

        #define COMPLEX_QMC_BODY_WITHOUT_TRANSFORM(complex_template) \
        protected: \
 \
            std::unique_ptr<Integrator<return_t, return_t>> get_real_integrator() \
            { \
                throw std::runtime_error("Separate integration of real and imaginary part has no advantage for this non-adaptive integrator. Please use \"together = true\"."); \
            }; \
 \
            std::function<secdecutil::UncorrelatedDeviation<complex_template<return_t>> \
            (const container_t&)> \
            get_together_integrate() \
            {\
                std::function<secdecutil::UncorrelatedDeviation<complex_template<return_t>>(const container_t&)> integrate_function = [this] (const container_t& integrand_container) \
                { \
                    ::integrators::result<complex_template<return_t>> result; \
                    result = ::integrators::Qmc<complex_template<return_t>,return_t>::integrate( \
                        integrand_container, \
                        integrand_container.number_of_integration_variables>0 ? integrand_container.number_of_integration_variables : 1 /* ensure dim > 0 */ \
                    ); \
                    return secdecutil::UncorrelatedDeviation<complex_template<return_t>> { result.integral, result.error }; \
                }; \
                return integrate_function; \
            }; \
 \
        public: \
 \
            using Integrator<complex_template<return_t>, return_t, container_t>::integrate; \
 \
            Qmc() { this->together = true;  };

        template<typename return_t, typename container_t, typename transform_t>
        struct Qmc<std::complex<return_t>, container_t, transform_t> : Integrator<std::complex<return_t>, return_t, container_t>, public ::integrators::Qmc<std::complex<return_t>,return_t>
        {
            COMPLEX_QMC_BODY_WITH_TRANSFORM(std::complex)
        };

        template<typename return_t, typename container_t>
        struct Qmc<std::complex<return_t>, container_t> : Integrator<std::complex<return_t>, return_t, container_t>, public ::integrators::Qmc<std::complex<return_t>,return_t>
        {
            COMPLEX_QMC_BODY_WITHOUT_TRANSFORM(std::complex)
        };

        template<typename return_t>
        struct Qmc<std::complex<return_t>> : Integrator<std::complex<return_t>, return_t>, public ::integrators::Qmc<std::complex<return_t>,return_t>
        {
            using container_t = secdecutil::IntegrandContainer<std::complex<return_t>, return_t const * const>;
            COMPLEX_QMC_BODY_WITHOUT_TRANSFORM(std::complex)
        };

        #ifdef SECDEC_WITH_CUDA
            template<typename return_t, typename container_t, typename transform_t>
            struct Qmc<thrust::complex<return_t>,container_t,transform_t> : Integrator<thrust::complex<return_t>, return_t, container_t>, public ::integrators::Qmc<thrust::complex<return_t>,return_t>
            {
                COMPLEX_QMC_BODY_WITH_TRANSFORM(thrust::complex)
            };

            template<typename return_t, typename container_t>
            struct Qmc<thrust::complex<return_t>,container_t> : Integrator<thrust::complex<return_t>, return_t, container_t>, public ::integrators::Qmc<thrust::complex<return_t>,return_t>
            {
                COMPLEX_QMC_BODY_WITHOUT_TRANSFORM(thrust::complex)
            };
        #endif

        #undef COMPLEX_QMC_BODY_WITH_TRANSFORM
        #undef COMPLEX_QMC_BODY_WITHOUT_TRANSFORM

    };

}

#endif

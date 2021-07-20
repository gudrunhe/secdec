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

        using ::integrators::U;

        template<typename,typename,U> using void_template = void;

        template<typename functor_t>
        struct QmcContainer : functor_t
        {
            QmcContainer(const functor_t& functor) : functor_t(functor) {};
            const unsigned number_of_integration_variables = functor_t::number_of_integration_variables >= 1 ? functor_t::number_of_integration_variables : 1; // ensure dim >= 1
        };

        template<typename T> struct remove_complex { using type = T; };
        template<typename T> struct remove_complex<std::complex<T>> { using type = T; };
        #ifdef SECDEC_WITH_CUDA
            template<typename T> struct remove_complex<thrust::complex<T>> { using type = T; };
        #endif

        // base template
        template<
                    typename return_t, U maxdim, template<typename,typename,U> class transform_t,
                    typename container_t = secdecutil::IntegrandContainer<return_t, typename remove_complex<return_t>::type const * const>,
                    template<typename,typename,U> class fitfunction_t = void_template
                >
        struct Qmc : Integrator<return_t,return_t,container_t>, public ::integrators::Qmc<return_t,return_t,maxdim,transform_t,fitfunction_t>
        {
        protected:

            std::function<secdecutil::UncorrelatedDeviation<return_t>(const container_t&)> get_integrate(); // define outside of the class to prevent inline

        public:

            using Integrator<return_t,return_t,container_t>::integrate;
            static constexpr bool cuda_compliant_integrator = true;

            Qmc(){
                // concat relevant generating vectors by default
                std::map<U,std::vector<U>> default_generatingvectors;
                std::map<U,std::vector<U>> next_map;
                if(maxdim<=6)
                {
                    next_map = ::integrators::generatingvectors::cbcpt_dn2_6();
                    default_generatingvectors.insert(next_map.begin(), next_map.end());
                    next_map = ::integrators::generatingvectors::cbcpt_cfftw1_6();
                    default_generatingvectors.insert(next_map.begin(), next_map.end());
                }
                if(maxdim<=10)
                {
                    next_map = ::integrators::generatingvectors::cbcpt_cfftw2_10();
                    default_generatingvectors.insert(next_map.begin(), next_map.end());
                }
                if(maxdim<=100)
                {
                    next_map = ::integrators::generatingvectors::cbcpt_dn1_100();
                    default_generatingvectors.insert(next_map.begin(), next_map.end());
                }
                this->generatingvectors = default_generatingvectors;
            };

            template<typename original_container_t>
            Qmc( const Qmc<return_t,maxdim,transform_t,original_container_t,fitfunction_t>& original) : Integrator<return_t,return_t,container_t>(), ::integrators::Qmc<return_t,return_t,maxdim,transform_t,fitfunction_t>(original)
            {};

            Qmc( const Qmc<return_t,maxdim,transform_t,container_t,fitfunction_t>& original) : Integrator<return_t,return_t,container_t>(), ::integrators::Qmc<return_t,return_t,maxdim,transform_t,fitfunction_t>(original)
            {};
        };
        template<typename return_t, U maxdim, template<typename,typename,U> class transform_t, typename container_t, template<typename,typename,U> class fitfunction_t>
        std::function<secdecutil::UncorrelatedDeviation<return_t>(const container_t&)> Qmc<return_t,maxdim,transform_t,container_t,fitfunction_t>::get_integrate()
        {
            std::function<secdecutil::UncorrelatedDeviation<return_t>(const container_t&)> integrate_function = [this] (const container_t& integrand_container)
            {
                QmcContainer<container_t> qmc_integrand_container(integrand_container);
                ::integrators::result<return_t> result = ::integrators::Qmc<return_t,return_t,maxdim,transform_t,fitfunction_t>::integrate(qmc_integrand_container);
                integrand_container.process_errors();
                return secdecutil::UncorrelatedDeviation<return_t> { result.integral, result.error };
            };
            return integrate_function;
        };

        // default fitfunction
        template<typename return_t, U maxdim, template<typename,typename,U> class transform_t, typename container_t>
        struct Qmc<return_t,maxdim,transform_t,container_t> : Integrator<return_t,return_t,container_t>, public ::integrators::Qmc<return_t,return_t,maxdim,transform_t>
        {
        protected:

            std::function<secdecutil::UncorrelatedDeviation<return_t>(const container_t&)> get_integrate(); // define outside of the class to prevent inline

        public:

            using Integrator<return_t,return_t,container_t>::integrate;
            static constexpr bool cuda_compliant_integrator = true;

            Qmc(){
                // concat relevant generating vectors by default
                std::map<U,std::vector<U>> default_generatingvectors;
                std::map<U,std::vector<U>> next_map;
                if(maxdim<=6)
                {
                    next_map = ::integrators::generatingvectors::cbcpt_dn2_6();
                    default_generatingvectors.insert(next_map.begin(), next_map.end());
                    next_map = ::integrators::generatingvectors::cbcpt_cfftw1_6();
                    default_generatingvectors.insert(next_map.begin(), next_map.end());
                }
                if(maxdim<=10)
                {
                    next_map = ::integrators::generatingvectors::cbcpt_cfftw2_10();
                    default_generatingvectors.insert(next_map.begin(), next_map.end());
                }
                if(maxdim<=100)
                {
                    next_map = ::integrators::generatingvectors::cbcpt_dn1_100();
                    default_generatingvectors.insert(next_map.begin(), next_map.end());
                }
                this->generatingvectors = default_generatingvectors;
            };

            template<typename original_container_t>
            Qmc( const Qmc<return_t,maxdim,transform_t,original_container_t>& original) : Integrator<return_t,return_t,container_t>(), ::integrators::Qmc<return_t,return_t,maxdim,transform_t>(original)
            {};

            Qmc( const Qmc<return_t,maxdim,transform_t,container_t>& original) : Integrator<return_t,return_t,container_t>(), ::integrators::Qmc<return_t,return_t,maxdim,transform_t>(original)
            {};

        };
        template<typename return_t, U maxdim, template<typename,typename,U> class transform_t, typename container_t>
        std::function<secdecutil::UncorrelatedDeviation<return_t>(const container_t&)> Qmc<return_t,maxdim,transform_t,container_t>::get_integrate()
        {
            std::function<secdecutil::UncorrelatedDeviation<return_t>(const container_t&)> integrate_function = [this] (const container_t& integrand_container)
            {
                QmcContainer<container_t> qmc_integrand_container(integrand_container);
                ::integrators::result<return_t> result = ::integrators::Qmc<return_t,return_t,maxdim,transform_t>::integrate(qmc_integrand_container);
                integrand_container.process_errors();
                return secdecutil::UncorrelatedDeviation<return_t> { result.integral, result.error };
            };
            return integrate_function;
        };


        /*
         * complex specilizations
         */

        #define COMPLEX_QMC_BODY(complex_template) \
            protected: \
                /* define outside of the class to prevent inline */ \
                std::function<secdecutil::UncorrelatedDeviation<complex_template<return_t>>(const container_t&)> get_together_integrate(); \
            public: \
                using Integrator<complex_template<return_t>, return_t, container_t>::integrate; \
                static constexpr bool cuda_compliant_integrator = true; \
                Qmc() { \
                    this->together = true; \
                    std::map<U,std::vector<U>> default_generatingvectors; \
                    std::map<U,std::vector<U>> next_map; \
                    if(maxdim<=6) \
                    { \
                        next_map = ::integrators::generatingvectors::cbcpt_dn2_6(); \
                        default_generatingvectors.insert(next_map.begin(), next_map.end()); \
                        next_map = ::integrators::generatingvectors::cbcpt_cfftw1_6(); \
                        default_generatingvectors.insert(next_map.begin(), next_map.end()); \
                    } \
                    if(maxdim<=10) \
                    { \
                        next_map = ::integrators::generatingvectors::cbcpt_cfftw2_10(); \
                        default_generatingvectors.insert(next_map.begin(), next_map.end()); \
                    } \
                    if(maxdim<=100) \
                    { \
                        next_map = ::integrators::generatingvectors::cbcpt_dn1_100(); \
                        default_generatingvectors.insert(next_map.begin(), next_map.end()); \
                    } \
                    this->generatingvectors = default_generatingvectors; \
                }; \
                template<typename original_container_t> \
                Qmc( const Qmc<complex_template<return_t>,maxdim,transform_t,original_container_t,fitfunction_t>& original) : \
                    Integrator<complex_template<return_t>,return_t,container_t>(), \
                    ::integrators::Qmc<complex_template<return_t>,return_t,maxdim,transform_t,fitfunction_t>(original) \
                    { this->together = true; }; \
                Qmc( const Qmc<complex_template<return_t>,maxdim,transform_t,container_t,fitfunction_t>& original) : \
                    Integrator<complex_template<return_t>,return_t,container_t>(), \
                    ::integrators::Qmc<complex_template<return_t>,return_t,maxdim,transform_t,fitfunction_t>(original) \
                    { this->together = true; };

        #define COMPLEX_QMC_BODY_WITHOUT_FITFUNCTION(complex_template) \
            protected: \
                /* define outside of the class to prevent inline */ \
                std::function<secdecutil::UncorrelatedDeviation<complex_template<return_t>>(const container_t&)> get_together_integrate(); \
            public: \
                using Integrator<complex_template<return_t>, return_t, container_t>::integrate; \
                static constexpr bool cuda_compliant_integrator = true; \
                Qmc() { \
                    this->together = true; \
                    std::map<U,std::vector<U>> default_generatingvectors; \
                    std::map<U,std::vector<U>> next_map; \
                    if(maxdim<=6) \
                    { \
                        next_map = ::integrators::generatingvectors::cbcpt_dn2_6(); \
                        default_generatingvectors.insert(next_map.begin(), next_map.end()); \
                        next_map = ::integrators::generatingvectors::cbcpt_cfftw1_6(); \
                        default_generatingvectors.insert(next_map.begin(), next_map.end()); \
                    } \
                    if(maxdim<=10) \
                    { \
                        next_map = ::integrators::generatingvectors::cbcpt_cfftw2_10(); \
                        default_generatingvectors.insert(next_map.begin(), next_map.end()); \
                    } \
                    if(maxdim<=100) \
                    { \
                        next_map = ::integrators::generatingvectors::cbcpt_dn1_100(); \
                        default_generatingvectors.insert(next_map.begin(), next_map.end()); \
                    } \
                    this->generatingvectors = default_generatingvectors; \
                }; \
                template<typename original_container_t> \
                Qmc( const Qmc<complex_template<return_t>,maxdim,transform_t,original_container_t>& original) : \
                    Integrator<complex_template<return_t>,return_t,container_t>(), \
                    ::integrators::Qmc<complex_template<return_t>,return_t,maxdim,transform_t>(original) \
                    { this->together = true; }; \
                Qmc( const Qmc<complex_template<return_t>,maxdim,transform_t,container_t>& original) : \
                    Integrator<complex_template<return_t>,return_t,container_t>(), \
                    ::integrators::Qmc<complex_template<return_t>,return_t,maxdim,transform_t>(original) \
                    { this->together = true; };

        #define COMPLEX_QMC_GET_TOGETHER_INTEGRATE_WITH_FITFUNCTION(complex_template) \
            template<typename return_t, U maxdim, template<typename,typename,U> class transform_t, typename container_t, template<typename,typename,U> class fitfunction_t> \
            std::function<secdecutil::UncorrelatedDeviation<complex_template<return_t>>(const container_t&)> \
            Qmc<complex_template<return_t>,maxdim,transform_t,container_t,fitfunction_t>::get_together_integrate() \
            { \
                std::function<secdecutil::UncorrelatedDeviation<complex_template<return_t>>(const container_t&)> integrate_function = [this] (const container_t& integrand_container) \
                { \
                    QmcContainer<container_t> qmc_integrand_container(integrand_container); \
                    ::integrators::result<complex_template<return_t>> result = ::integrators::Qmc<complex_template<return_t>,return_t,maxdim,transform_t,fitfunction_t>::integrate(qmc_integrand_container); \
                    integrand_container.process_errors(); \
                    return secdecutil::UncorrelatedDeviation<complex_template<return_t>> { result.integral, result.error }; \
                }; \
                return integrate_function; \
            }; \

        #define COMPLEX_QMC_GET_TOGETHER_INTEGRATE_WITHOUT_FITFUNCTION(complex_template) \
            template<typename return_t, U maxdim, template<typename,typename,U> class transform_t, typename container_t> \
            std::function<secdecutil::UncorrelatedDeviation<complex_template<return_t>>(const container_t&)> Qmc<complex_template<return_t>,maxdim,transform_t,container_t>::get_together_integrate() \
            { \
                std::function<secdecutil::UncorrelatedDeviation<complex_template<return_t>>(const container_t&)> integrate_function = [this] (const container_t& integrand_container) \
                { \
                    QmcContainer<container_t> qmc_integrand_container(integrand_container); \
                    ::integrators::result<complex_template<return_t>> result = ::integrators::Qmc<complex_template<return_t>,return_t,maxdim,transform_t>::integrate(qmc_integrand_container); \
                    integrand_container.process_errors(); \
                    return secdecutil::UncorrelatedDeviation<complex_template<return_t>> { result.integral, result.error }; \
                }; \
                return integrate_function; \
            }; \

        template<typename return_t, U maxdim, template<typename,typename,U> class transform_t, typename container_t, template<typename,typename,U> class fitfunction_t>
        struct Qmc<std::complex<return_t>,maxdim,transform_t,container_t,fitfunction_t> : Integrator<std::complex<return_t>,return_t,container_t>,
            public ::integrators::Qmc<std::complex<return_t>,return_t,maxdim,transform_t,fitfunction_t>
        {
            COMPLEX_QMC_BODY(std::complex)
        };
        COMPLEX_QMC_GET_TOGETHER_INTEGRATE_WITH_FITFUNCTION(std::complex)

        template<typename return_t, U maxdim, template<typename,typename,U> class transform_t, typename container_t>
        struct Qmc<std::complex<return_t>,maxdim,transform_t,container_t> : Integrator<std::complex<return_t>,return_t,container_t>,
            public ::integrators::Qmc<std::complex<return_t>,return_t,maxdim,transform_t>
        {
            COMPLEX_QMC_BODY_WITHOUT_FITFUNCTION(std::complex)
        };
        COMPLEX_QMC_GET_TOGETHER_INTEGRATE_WITHOUT_FITFUNCTION(std::complex)

        #ifdef SECDEC_WITH_CUDA
            template<typename return_t, U maxdim, template<typename,typename,U> class transform_t, typename container_t, template<typename,typename,U> class fitfunction_t>
            struct Qmc<thrust::complex<return_t>,maxdim,transform_t,container_t,fitfunction_t> : Integrator<thrust::complex<return_t>,return_t,container_t>,
                public ::integrators::Qmc<thrust::complex<return_t>,return_t,maxdim,transform_t,fitfunction_t>
            {
                COMPLEX_QMC_BODY(thrust::complex)
            };
            COMPLEX_QMC_GET_TOGETHER_INTEGRATE_WITH_FITFUNCTION(thrust::complex)

            template<typename return_t, U maxdim, template<typename,typename,U> class transform_t, typename container_t>
            struct Qmc<thrust::complex<return_t>,maxdim,transform_t,container_t> : Integrator<thrust::complex<return_t>, return_t, container_t>,
                public ::integrators::Qmc<thrust::complex<return_t>,return_t,maxdim,transform_t>
            {
                COMPLEX_QMC_BODY_WITHOUT_FITFUNCTION(thrust::complex)
            };
            COMPLEX_QMC_GET_TOGETHER_INTEGRATE_WITHOUT_FITFUNCTION(thrust::complex)
        #endif

        #undef COMPLEX_QMC_GET_TOGETHER_INTEGRATE
        #undef COMPLEX_QMC_BODY_WITH_FITFUNCTION
        #undef COMPLEX_QMC_BODY_WITHOUT_FITFUNCTION

    };

}

#endif

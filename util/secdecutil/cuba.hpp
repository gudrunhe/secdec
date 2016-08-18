#ifndef SecDecUtil_cuba_hpp_included
#define SecDecUtil_cuba_hpp_included

/*
 * This file implements a convenience wrapper around
 * the CUBA integrator library.
 */

#include <complex>
#include <cuba.h>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/uncertainties.hpp>

namespace secdecutil
{
    namespace cuba
    {
        #define VEGAS_BODY \
            constexpr static long long int nvec = 1; \
            T epsrel; \
            T epsabs; \
            int flags; \
            int seed; \
            long long int mineval; \
            long long int maxeval; \
            long long int nstart; \
            long long int nincrease; \
            long long int nbatch; \
            constexpr static int gridno = 0; \
            constexpr static char * statefile = nullptr; \
            constexpr static void* spin = nullptr; \
            \
            Vegas \
            ( \
                T epsrel = 1e-2, \
                T epsabs = 1e-7, \
                int flags = 0, \
                int seed = 0, \
                long long int mineval = 0, \
                long long int maxeval = 100000, \
                long long int nstart = 1000, \
                long long int nincrease = 500, \
                long long int nbatch = 1000 \
            ) : \
                epsrel(epsrel),epsabs(epsabs), \
                flags(flags),seed(seed),mineval(mineval),maxeval(maxeval), \
                nstart(nstart),nincrease(nincrease),nbatch(nbatch) \
            {};

        #define VEGAS_INTEGRATE_BODY \
        /* Cuba output values */ \
        std::array<T, ncomp> integral; \
        std::array<T, ncomp> error; \
        std::array<T, ncomp> prob; \
        int fail; \
        long long int neval; \
        \
        /* Cuba call */ \
        ::llVegas \
        ( \
            integrand_container.number_of_integration_variables, \
            ncomp, \
            cuba_integrand_prototype, \
            reinterpret_cast<void*>(&integrand_container), \
            nvec, \
            epsrel, \
            epsabs, \
            flags, \
            seed, \
            mineval, \
            maxeval, \
            nstart, \
            nincrease, \
            nbatch, \
            gridno, \
            statefile, \
            spin, \
            &neval, \
            &fail, \
            integral.data(), \
            error.data(), \
            prob.data() \
        );

        template<typename T>
        struct Vegas
        {
            static int cuba_integrand_prototype(const int *ndim, const T integration_variables[], const int *ncomp, T result[], void *userdata)
            {
                auto& integrand_container = *( reinterpret_cast<secdecutil::IntegrandContainer<T, T const * const> *>(userdata) );
                result[0] = integrand_container.integrand(integration_variables);
                return 0;
            };

            constexpr static int ncomp = 1;
            VEGAS_BODY

            std::function<secdecutil::GaussianUncertainty<T>(secdecutil::IntegrandContainer<T, T const * const>)>
            integrate =
            [ this ] (secdecutil::IntegrandContainer<T, T const * const> integrand_container)
            {
                VEGAS_INTEGRATE_BODY
                return secdecutil::GaussianUncertainty<T>(integral.at(0),error.at(0));
            };
        };
        template<typename T>
        struct Vegas<std::complex<T>>
        {
            static int cuba_integrand_prototype(const int *ndim, const T integration_variables[], const int *ncomp, T result[], void *userdata)
            {
                auto& integrand_container = *( reinterpret_cast<secdecutil::IntegrandContainer<std::complex<T>, T const * const> *>(userdata) );
                std::complex<T> evaluated_integrand = integrand_container.integrand(integration_variables);
                result[0] = evaluated_integrand.real();
                result[1] = evaluated_integrand.imag();
                return 0;
            };

            constexpr static int ncomp = 2;
            VEGAS_BODY

            std::function<secdecutil::GaussianUncertainty<std::complex<T>>(secdecutil::IntegrandContainer<std::complex<T>, T const * const>)>
            integrate =
            [ this ] (secdecutil::IntegrandContainer<std::complex<T>, T const * const> integrand_container)
            {
                VEGAS_INTEGRATE_BODY
                return secdecutil::GaussianUncertainty<std::complex<T>>({integral.at(0),integral.at(1)},{error.at(0),error.at(1)});
            };
        };
        #undef VEGAS_BODY
        #undef VEGAS_INTEGRATE_BODY
    };
};
#endif

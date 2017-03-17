#ifndef SecDecUtil_cuba_hpp_included
#define SecDecUtil_cuba_hpp_included

/*
 * This file implements a convenience wrapper around
 * the CUBA integrator library.
 */

#include <array>
#include <cmath>
#include <complex>
#include <cuba.h>
#include <iostream>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/integrators/integrator.hpp>
#include <secdecutil/uncertainties.hpp>
#include <vector>

namespace secdecutil
{
    namespace cuba
    {
     /*
      * Base class for the c++ wrapper of the Cuba library.
      * All Cuba integrators should inherit from this class.
      *
      * Cuba always uses the float type "cubareal" as defined in "cuba.h",
      * which can be different from the type "T" used by the IntegrandContainer.
      */

      #define CUBA_STRUCT_BODY \
        int ndim; \
        void * userdata; \
        virtual void call_cuba() = 0; \
        std::array<cubareal,ncomp> integral; \
        std::array<cubareal,ncomp> error; \
        std::array<cubareal,ncomp> prob;

      #define CUBA_INTEGRATE_BODY \
        ndim = integrand_container.number_of_integration_variables; \
        /* must have at least one integration variable */ \
        if (ndim == 0) \
            ndim = 1; \
        /* nasty removal of constness --> restored in cuba_integrand_prototype */ \
        userdata = const_cast<void*>( reinterpret_cast<const void*>(&integrand_container) ); \
        call_cuba(); \
        if (flags & 3) \
        { \
              std::cout << std::endl << "Cuba chi2-probability:" << std::endl; \
              for ( int i=0; i<ncomp; i++) \
                std::cout << "[" << i+1 << "] " << prob[i] << std::endl; \
              std::cout << std::endl << std::endl; \
        }

      // real version for general type
      template<typename T>
      struct CubaIntegrator : Integrator<T,T>
      {
      protected:
        static int cuba_integrand_prototype(const int *ndim, const cubareal integration_variables[], const int *ncomp, cubareal result[], void *userdata)
        {
          auto& integrand_container = *( reinterpret_cast<const secdecutil::IntegrandContainer<T, T const * const> *>(userdata) );

          /* "integration_variables" is an array of "cubareal", but the integrand expects type T
           * --> copy them into a vector of type T using the iterator contructor.
           * During the copy operation the type is converted implicitly.
           */
          const std::vector<T> converted_integration_variables(integration_variables, integration_variables + (*ndim)*sizeof(cubareal));

          // initialize result with NaN --> result will be NaN if integrand throws an error
          result[0] = std::nan("");

          result[0] = integrand_container.integrand(converted_integration_variables.data()); // implicit conversion of result from type T to cubareal

          return 0;
        };
        static const int ncomp = 1;
        CUBA_STRUCT_BODY

        std::function<secdecutil::UncorrelatedDeviation<T>
          (const secdecutil::IntegrandContainer<T, T const * const>&)> get_integrate()
        {
          return [this] (const secdecutil::IntegrandContainer<T, T const * const>& integrand_container)
            {
              CUBA_INTEGRATE_BODY
              return secdecutil::UncorrelatedDeviation<T>(integral.at(0),error.at(0));
            };
        };
      public:
        int flags;
      };

      // real version specialized for cubareal
      template<>
      struct CubaIntegrator<cubareal> : Integrator<cubareal,cubareal>
      {
      protected:
        static int cuba_integrand_prototype(const int *ndim, const cubareal integration_variables[], const int *ncomp, cubareal result[], void *userdata)
        {
          auto& integrand_container = *( reinterpret_cast<const secdecutil::IntegrandContainer<cubareal, cubareal const * const> *>(userdata) );

          // initialize result with NaN --> result will be NaN if integrand throws an error
          result[0] = std::nan("");

          result[0] = integrand_container.integrand(integration_variables); // pass array "integration_variables" directly

          return 0;
        };
        static const int ncomp = 1;
        CUBA_STRUCT_BODY

        std::function<secdecutil::UncorrelatedDeviation<cubareal>
          (const secdecutil::IntegrandContainer<cubareal, cubareal const * const>&)> get_integrate()
        {
          return [this] (const secdecutil::IntegrandContainer<cubareal, cubareal const * const>& integrand_container)
            {
              CUBA_INTEGRATE_BODY
              return secdecutil::UncorrelatedDeviation<cubareal>(integral.at(0),error.at(0));
            };
        };
      public:
        int flags;
      };


      // complex version
      template<typename T>
      struct CubaIntegrator<std::complex<T>> : Integrator<std::complex<T>,T>
      {
      protected:
        static int cuba_integrand_prototype(const int *ndim, const cubareal integration_variables[], const int *ncomp, cubareal result[], void *userdata)
        {
          auto& integrand_container = *( reinterpret_cast<const secdecutil::IntegrandContainer<std::complex<T>, T const * const> *>(userdata) );

          /* "integration_variables" is an array of "cubareal", but the integrand expects type T
           * --> copy them into a vector of type T using the iterator contructor.
           * During the copy operation the type is converted implicitly.
           */
          const std::vector<T> converted_integration_variables(integration_variables, integration_variables + (*ndim)*sizeof(cubareal));

          // initialize result with NaN --> result will be NaN if integrand throws an error
          result[0] = result[1] = std::nan("");

          std::complex<T> evaluated_integrand = integrand_container.integrand(converted_integration_variables.data());

          result[0] = evaluated_integrand.real(); // implicit conversion of result from type T to cubareal
          result[1] = evaluated_integrand.imag(); // implicit conversion of result from type T to cubareal

          return 0;
        };
        static const int ncomp = 2;
        CUBA_STRUCT_BODY

        std::function<secdecutil::UncorrelatedDeviation<std::complex<T>>
          (const secdecutil::IntegrandContainer<std::complex<T>, T const * const>&)>
          get_together_integrate()
          {
            return [this] (const secdecutil::IntegrandContainer<std::complex<T>, T const * const>& integrand_container) {
              CUBA_INTEGRATE_BODY
              return secdecutil::UncorrelatedDeviation<std::complex<T>>({integral.at(0),integral.at(1)},{error.at(0),error.at(1)});
            };
        };
      public:
        int flags;
      };

      // complex version specialized for cubareal
      template<>
      struct CubaIntegrator<std::complex<cubareal>> : Integrator<std::complex<cubareal>,cubareal>
      {
      protected:
        static int cuba_integrand_prototype(const int *ndim, const cubareal integration_variables[], const int *ncomp, cubareal result[], void *userdata)
        {
          auto& integrand_container = *( reinterpret_cast<const secdecutil::IntegrandContainer<std::complex<cubareal>, cubareal const * const> *>(userdata) );

          // initialize result with NaN --> result will be NaN if integrand throws an error
          result[0] = result[1] = std::nan("");

          std::complex<cubareal> evaluated_integrand = integrand_container.integrand(integration_variables); // pass array "integration_variables" directly

          result[0] = evaluated_integrand.real();
          result[1] = evaluated_integrand.imag();

          return 0;
        };
        static const int ncomp = 2;
        CUBA_STRUCT_BODY

        std::function<secdecutil::UncorrelatedDeviation<std::complex<cubareal>>
          (const secdecutil::IntegrandContainer<std::complex<cubareal>, cubareal const * const>&)>
          get_together_integrate()
          {
            return [this] (const secdecutil::IntegrandContainer<std::complex<cubareal>, cubareal const * const>& integrand_container) {
              CUBA_INTEGRATE_BODY
              return secdecutil::UncorrelatedDeviation<std::complex<cubareal>>({integral.at(0),integral.at(1)},{error.at(0),error.at(1)});
            };
        };
      public:
        int flags;
      };

      #undef CUBA_STRUCT_BODY
      #undef CUBA_INTEGRATE_BODY

      ////////////////////////////////////////// Vegas //////////////////////////////////////////

        #define VEGAS_STRUCT_BODY \
            constexpr static long long int nvec = 1; \
            cubareal epsrel; \
            cubareal epsabs; \
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
                cubareal epsrel = 1e-2, \
                cubareal epsabs = 1e-7, \
                int flags = 0, \
                int seed = 0, \
                long long int mineval = 0, \
                long long int maxeval = 1e6, \
                long long int nstart = 1000, \
                long long int nincrease = 500, \
                long long int nbatch = 1000 \
            ) : \
                epsrel(epsrel),epsabs(epsabs), \
                seed(seed),mineval(mineval),maxeval(maxeval), \
                nstart(nstart),nincrease(nincrease),nbatch(nbatch) \
            { \
                this->flags = flags; \
            };

        #define VEGAS_INTEGRATE_BODY \
        /* Cuba output values */ \
        int fail; \
        long long int neval; \
        \
        /* Cuba call */ \
        ::llVegas \
        ( \
            this->ndim, \
            this->ncomp, \
            this->cuba_integrand_prototype, \
            this->userdata, \
            nvec, \
            epsrel, \
            epsabs, \
            this->flags, \
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
            this->integral.data(), \
            this->error.data(), \
            this->prob.data() \
        );

      template <typename T>
      struct Vegas : CubaIntegrator<T> {
        VEGAS_STRUCT_BODY
        void call_cuba(){
          VEGAS_INTEGRATE_BODY
        };
      };

      template <typename T>
      struct Vegas<std::complex<T>> : CubaIntegrator<std::complex<T>> {
        std::unique_ptr<Integrator<T,T>> get_real_integrator(){
          return std::unique_ptr<Integrator<T,T>>(
                                                     new Vegas<T>(
                                                                     epsrel,epsabs,
                                                                     this->flags,seed,mineval,maxeval,
                                                                     nstart,nincrease,nbatch
                                                                 )
                                                 );
        };
        VEGAS_STRUCT_BODY
        void call_cuba(){
          VEGAS_INTEGRATE_BODY
        };
      };

      #undef VEGAS_STRUCT_BODY
      #undef VEGAS_INTEGRATE_BODY


      ////////////////////////////////////////// Suave //////////////////////////////////////////

        #define SUAVE_STRUCT_BODY \
            constexpr static long long int nvec = 1; \
            cubareal epsrel; \
            cubareal epsabs; \
            int seed; \
            long long int mineval; \
            long long int maxeval; \
            long long int nnew; \
            long long int nmin; \
            cubareal flatness; \
            constexpr static char * statefile = nullptr; \
            constexpr static void* spin = nullptr; \
            \
            Suave \
            ( \
                cubareal epsrel = 1e-2, \
                cubareal epsabs = 1e-7, \
                int flags = 0, \
                int seed = 0, \
                long long int mineval = 0, \
                long long int maxeval = 1e6, \
                long long int nnew = 1000, \
                long long int nmin = 10, \
                cubareal flatness = 25. \
            ) : \
                epsrel(epsrel),epsabs(epsabs), \
                seed(seed),mineval(mineval),maxeval(maxeval), \
                nnew(nnew),nmin(nmin),flatness(flatness)        \
            { \
                this->flags = flags; \
            };

        #define SUAVE_INTEGRATE_BODY \
        /* Cuba output values */ \
        int fail; \
        long long int neval; \
        int nregions; \
        \
        /* Cuba call */ \
        ::llSuave \
        ( \
            this->ndim, \
            this->ncomp, \
            this->cuba_integrand_prototype, \
            this->userdata, \
            nvec, \
            epsrel, \
            epsabs, \
            this->flags, \
            seed, \
            mineval, \
            maxeval, \
            nnew, \
            nmin, \
            flatness, \
            statefile, \
            spin, \
            &nregions, \
            &neval, \
            &fail, \
            this->integral.data(), \
            this->error.data(), \
            this->prob.data() \
        );

      template <typename T>
      struct Suave : CubaIntegrator<T> {
        SUAVE_STRUCT_BODY
        void call_cuba(){
          SUAVE_INTEGRATE_BODY
        }
      };
      template <typename T>
      struct Suave<std::complex<T>> : CubaIntegrator<std::complex<T>> {
        std::unique_ptr<Integrator<T,T>> get_real_integrator(){
          return std::unique_ptr<Integrator<T,T>>(
                                                     new Suave<T>(
                                                                     epsrel,epsabs,
                                                                     this->flags,seed,mineval,maxeval,
                                                                     nnew,nmin,flatness
                                                                 )
                                                 );
        };
        SUAVE_STRUCT_BODY
        void call_cuba(){
          SUAVE_INTEGRATE_BODY
        }
      };

      #undef SUAVE_STRUCT_BODY
      #undef SUAVE_INTEGRATE_BODY

      ////////////////////////////////////////// Divonne //////////////////////////////////////////

        #define DIVONNE_STRUCT_BODY \
            constexpr static long long int nvec = 1; \
            cubareal epsrel; \
            cubareal epsabs; \
            int seed; \
            long long int mineval; \
            long long int maxeval; \
            int key1; \
            int key2; \
            int key3; \
            int maxpass; \
            cubareal border; \
            cubareal maxchisq; \
            cubareal mindeviation; \
            constexpr static long long int ngiven = 0; \
            constexpr static int ldxgiven = 0; \
            constexpr static cubareal * xgiven = nullptr; \
            constexpr static long long int nextra = 0; \
            constexpr static peakfinder_t peakfinder = nullptr; \
            constexpr static char * statefile = nullptr; \
            constexpr static void * spin = nullptr; \
            \
            Divonne \
            ( \
                cubareal epsrel = 1e-2, \
                cubareal epsabs = 1e-7, \
                int flags = 0, \
                int seed = 0, \
                long long int mineval = 0, \
                long long int maxeval = 1e6, \
                int key1 = 2000, \
                int key2 = 1, \
                int key3 = 1, \
                int maxpass = 4, \
                cubareal border = 0., \
                cubareal maxchisq = 1., \
                cubareal mindeviation = .15 \
            ) : \
                epsrel(epsrel),epsabs(epsabs), \
                seed(seed),mineval(mineval),maxeval(maxeval), \
                key1(key1), key2(key2), key3(key3), maxpass(maxpass), \
                border(border), maxchisq(maxchisq), mindeviation(mindeviation) \
            { \
                this->flags = flags; \
            };

        #define DIVONNE_INTEGRATE_BODY \
        /* Cuba output values */ \
        int nregions; \
        long long int neval; \
        int fail; \
        \
        /* Cuba call */ \
        ::llDivonne \
        ( \
            this->ndim, \
            this->ncomp, \
            this->cuba_integrand_prototype, \
            this->userdata, \
            nvec, \
            epsrel, \
            epsabs, \
            this->flags, \
            seed, \
            mineval, \
            maxeval, \
            key1, \
            key2, \
            key3, \
            maxpass, \
            border, \
            maxchisq, \
            mindeviation, \
            ngiven, \
            ldxgiven, \
            xgiven, \
            nextra, \
            peakfinder, \
            statefile, \
            spin, \
            &nregions, \
            &neval, \
            &fail, \
            this->integral.data(), \
            this->error.data(), \
            this->prob.data() \
        );

      template <typename T>
      struct Divonne : CubaIntegrator<T> {
        DIVONNE_STRUCT_BODY
        void call_cuba(){
          DIVONNE_INTEGRATE_BODY
        }
      };
      template <typename T>
      struct Divonne<std::complex<T>> : CubaIntegrator<std::complex<T>> {
        std::unique_ptr<Integrator<T,T>> get_real_integrator(){
          return std::unique_ptr<Integrator<T,T>>(
                                                     new Divonne<T>(
                                                                       epsrel,epsabs,
                                                                       this->flags,seed,mineval,maxeval,
                                                                       key1, key2, key3, maxpass,
                                                                       border, maxchisq, mindeviation
                                                                   )
                                                 );
        };
        DIVONNE_STRUCT_BODY
        void call_cuba(){
          DIVONNE_INTEGRATE_BODY
        }
      };

      #undef DIVONNE_STRUCT_BODY
      #undef DIVONNE_INTEGRATE_BODY

      ////////////////////////////////////////// Cuhre //////////////////////////////////////////

        #define CUHRE_STRUCT_BODY \
            constexpr static long long int nvec = 1; \
            cubareal epsrel; \
            cubareal epsabs; \
            long long int mineval; \
            long long int maxeval; \
            int key; \
            constexpr static char * statefile = nullptr; \
            constexpr static void* spin = nullptr; \
            \
            Cuhre \
            ( \
                cubareal epsrel = 1e-2, \
                cubareal epsabs = 1e-7, \
                int flags = 0, \
                long long int mineval = 0, \
                long long int maxeval = 1e6, \
                int key = 0 \
            ) : \
                epsrel(epsrel),epsabs(epsabs), \
                mineval(mineval),maxeval(maxeval),key(key) \
            { \
                this->flags = flags; \
            };

        #define CUHRE_INTEGRATE_BODY \
        /* Cuba output values */ \
        int nregions; \
        long long int neval; \
        int fail; \
        \
        /* Cuba call */ \
        ::llCuhre \
        ( \
            this->ndim, \
            this->ncomp, \
            this->cuba_integrand_prototype, \
            this->userdata, \
            nvec, \
            epsrel, \
            epsabs, \
            this->flags, \
            mineval, \
            maxeval, \
            key, \
            statefile, \
            spin, \
            &nregions, \
            &neval, \
            &fail, \
            this->integral.data(), \
            this->error.data(), \
            this->prob.data() \
        );


      template <typename T>
      struct Cuhre : CubaIntegrator<T> {
        CUHRE_STRUCT_BODY
        void call_cuba(){
          CUHRE_INTEGRATE_BODY
        }
      };
      template <typename T>
      struct Cuhre<std::complex<T>> : CubaIntegrator<std::complex<T>> {
        std::unique_ptr<Integrator<T,T>> get_real_integrator(){
          return std::unique_ptr<Integrator<T,T>>(
                                                     new Cuhre<T>(
                                                                     epsrel,epsabs,
                                                                     this->flags,mineval,maxeval,key
                                                                 )
                                                 );
        };
        CUHRE_STRUCT_BODY
        void call_cuba(){
          CUHRE_INTEGRATE_BODY
        }
      };

      #undef CUHRE_STRUCT_BODY
      #undef CUHRE_INTEGRATE_BODY
    }
}

#endif

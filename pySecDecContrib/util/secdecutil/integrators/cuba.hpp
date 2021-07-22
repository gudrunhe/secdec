#ifndef SecDecUtil_cuba_hpp_included
#define SecDecUtil_cuba_hpp_included

/*
 * This file implements a convenience wrapper around
 * the CUBA integrator library.
 */

#include <array>
#include <cmath>
#ifdef SECDEC_WITH_CUDA
    #include <thrust/complex.h>
#endif
#include <complex>
#include <cuba.h>
#include <iostream>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/integrators/integrator.hpp>
#include <secdecutil/uncertainties.hpp>
#include <vector>
#include <algorithm> // std::max
#include <string> 

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
        void * userdata = reinterpret_cast<void*>( &typed_userdata ); \
        virtual void call_cuba() = 0; \
        std::array<cubareal,ncomp> integral; \
        std::array<cubareal,ncomp> error; \
        std::array<cubareal,ncomp> prob; 

      #define CUBA_INTEGRATE_BODY \
        ndim = integrand_container.number_of_integration_variables; \
        /* must have at least two integration variables, otherwise Cuhre and Divonne do not work */ \
        if (ndim <= 1) \
            ndim = 2; \
        typed_userdata.integrand_container = &integrand_container; \
        if (flags & 3 and zero_border != 0) \
        { \
              std::cerr << "integrating with zero_border = " << zero_border << std::endl; \
        } \
        call_cuba(); \
        if (flags & 3) \
        { \
              std::cerr << std::endl << "Cuba chi2-probability:" << std::endl; \
              for ( int i=0; i<ncomp; i++) \
                std::cerr << "[" << i+1 << "] " << prob[i] << std::endl; \
              std::cerr << std::endl << std::endl; \
        }

      // real version for general type
      template<typename T>
      struct CubaIntegrator : Integrator<T,T>
      {
      protected:
        template<bool have_zero_border>
        static int cuba_integrand_prototype(const int *ndim, const cubareal integration_variables[], const int *ncomp, cubareal result[], void *userdata)
        {
          auto& typed_userdata = *( reinterpret_cast<const userdata_t *>(userdata) );

          /* "integration_variables" is an array of "cubareal", but the integrand expects type T
           * --> copy them into a vector of type T using the iterator contructor.
           * During the copy operation the type is converted implicitly.
           * Implement "zero_border".
           */
          T bordered_integration_variables[*ndim];
          for (int i = 0 ; i < *ndim ; ++i)
          {
              if (have_zero_border)
                  bordered_integration_variables[i] = integration_variables[i] < typed_userdata.zero_border ? typed_userdata.zero_border : integration_variables[i];
              else
                  bordered_integration_variables[i] = integration_variables[i];
          }

          // initialize result with NaN --> result will be NaN if integrand throws an error
          result[0] = std::nan("");

          // implicit conversion of result from type T to cubareal
          result[0] = (*typed_userdata.integrand_container)(bordered_integration_variables);

          return 0;
        }
        static const int ncomp = 1;
        struct userdata_t
        {
            const secdecutil::IntegrandContainer<T, T const * const> * integrand_container;
            const cubareal& zero_border;
        } typed_userdata{nullptr,zero_border};
        CUBA_STRUCT_BODY

        std::function<secdecutil::UncorrelatedDeviation<T>
          (const secdecutil::IntegrandContainer<T, T const * const>&)> get_integrate()
        {
          return [this] (const secdecutil::IntegrandContainer<T, T const * const>& integrand_container)
            {
              CUBA_INTEGRATE_BODY
              integrand_container.process_errors();
              return secdecutil::UncorrelatedDeviation<T>(integral.at(0),error.at(0));
            };
        };
      public:
        int flags;
        cubareal zero_border;
        static constexpr bool cuda_compliant_integrator = false;
      };

      // real version specialized for cubareal
      template<>
      struct CubaIntegrator<cubareal> : Integrator<cubareal,cubareal>
      {
      protected:
        template<bool have_zero_border>
        static int cuba_integrand_prototype(const int *ndim, const cubareal integration_variables[], const int *ncomp, cubareal result[], void *userdata)
        {
          auto& typed_userdata = *( reinterpret_cast<const userdata_t *>(userdata) );

          // initialize result with NaN --> result will be NaN if integrand throws an error
          result[0] = std::nan("");

          // Implement "zero_border".
          if (have_zero_border)
          {
              cubareal bordered_integration_variables[*ndim];
              for (int i = 0 ; i < *ndim ; ++i)
                  bordered_integration_variables[i] = integration_variables[i] < typed_userdata.zero_border ? typed_userdata.zero_border : integration_variables[i];
              result[0] = (*typed_userdata.integrand_container)(bordered_integration_variables);
              return 0;
          } else {
              result[0] = (*typed_userdata.integrand_container)(integration_variables); // pass array "integration_variables" directly
              return 0;
          }
        }
        static const int ncomp = 1;
        struct userdata_t
        {
            const secdecutil::IntegrandContainer<cubareal, cubareal const * const> * integrand_container;
            const cubareal& zero_border;
        } typed_userdata{nullptr,zero_border};
        CUBA_STRUCT_BODY

        std::function<secdecutil::UncorrelatedDeviation<cubareal>
          (const secdecutil::IntegrandContainer<cubareal, cubareal const * const>&)> get_integrate()
        {
          return [this] (const secdecutil::IntegrandContainer<cubareal, cubareal const * const>& integrand_container)
            {
              CUBA_INTEGRATE_BODY
              integrand_container.process_errors();
              return secdecutil::UncorrelatedDeviation<cubareal>(integral.at(0),error.at(0));
            };
        };
      public:
        int flags;
        cubareal zero_border;
        static constexpr bool cuda_compliant_integrator = false;
      };


      // complex version
      #define COMPLEX_CUBA_INTEGRATOR(complex_template) \
      template<typename T> \
      struct CubaIntegrator<complex_template<T>> : Integrator<complex_template<T>,T> \
      { \
      protected: \
        template<bool have_zero_border> \
        static int cuba_integrand_prototype(const int *ndim, const cubareal integration_variables[], const int *ncomp, cubareal result[], void *userdata) \
        { \
          auto& typed_userdata = *( reinterpret_cast<const userdata_t *>(userdata) ); \
 \
          /* "integration_variables" is an array of "cubareal", but the integrand expects type T \
           * --> copy them into a vector of type T using the iterator contructor. \
           * During the copy operation the type is converted implicitly. \
           * Implement "zero_border". \
           */ \
          T bordered_integration_variables[*ndim]; \
          for (int i = 0 ; i < *ndim ; ++i) \
          { \
              if (have_zero_border) \
                  bordered_integration_variables[i] = integration_variables[i] < typed_userdata.zero_border ? typed_userdata.zero_border : integration_variables[i]; \
              else \
                  bordered_integration_variables[i] = integration_variables[i]; \
          } \
 \
          /* initialize result with NaN --> result will be NaN if integrand throws an error */ \
          result[0] = result[1] = std::nan(""); \
 \
          complex_template<T> evaluated_integrand = (*typed_userdata.integrand_container)(bordered_integration_variables); \
 \
          /* implicit conversion of result from type T to cubareal */ \
          result[0] = evaluated_integrand.real(); \
          result[1] = evaluated_integrand.imag(); \
 \
          return 0; \
        } \
        static const int ncomp = 2; \
        struct userdata_t \
        { \
            const secdecutil::IntegrandContainer<complex_template<T>, T const * const> * integrand_container; \
            const cubareal& zero_border; \
        } typed_userdata{nullptr,zero_border}; \
        CUBA_STRUCT_BODY \
 \
        std::function<secdecutil::UncorrelatedDeviation<complex_template<T>> \
          (const secdecutil::IntegrandContainer<complex_template<T>, T const * const>&)> \
          get_together_integrate() \
          { \
            return [this] (const secdecutil::IntegrandContainer<complex_template<T>, T const * const>& integrand_container) { \
              CUBA_INTEGRATE_BODY \
              integrand_container.process_errors(); \
              return secdecutil::UncorrelatedDeviation<complex_template<T>>({integral.at(0),integral.at(1)},{error.at(0),error.at(1)}); \
            }; \
        }; \
      public: \
        int flags; \
        cubareal zero_border; \
        static constexpr bool cuda_compliant_integrator = false; \
      };
      COMPLEX_CUBA_INTEGRATOR(std::complex)
      #ifdef SECDEC_WITH_CUDA
        COMPLEX_CUBA_INTEGRATOR(thrust::complex)
      #endif
      #undef COMPLEX_CUBA_INTEGRATOR

      // complex version specialized for cubareal
      #define COMPLEX_CUBAREAL_INTEGRATOR(complex_template) \
      template<> \
      struct CubaIntegrator<complex_template<cubareal>> : Integrator<complex_template<cubareal>,cubareal> \
      { \
      protected: \
        template<bool have_zero_border> \
        static int cuba_integrand_prototype(const int *ndim, const cubareal integration_variables[], const int *ncomp, cubareal result[], void *userdata) \
        { \
          auto& typed_userdata = *( reinterpret_cast<const userdata_t *>(userdata) ); \
 \
          /* initialize result with NaN --> result will be NaN if integrand throws an error */ \
          result[0] = result[1] = std::nan(""); \
 \
          /* Implement "zero_border". */ \
          if (have_zero_border) \
          { \
              cubareal bordered_integration_variables[*ndim]; \
              for (int i = 0 ; i < *ndim ; ++i) \
                  bordered_integration_variables[i] = integration_variables[i] < typed_userdata.zero_border ? typed_userdata.zero_border : integration_variables[i]; \
              complex_template<cubareal> evaluated_integrand = (*typed_userdata.integrand_container)(bordered_integration_variables); \
              result[0] = evaluated_integrand.real(); \
              result[1] = evaluated_integrand.imag(); \
              return 0; \
          } else { \
              /* pass array "integration_variables" directly */ \
              complex_template<cubareal> evaluated_integrand = (*typed_userdata.integrand_container)(integration_variables); \
              result[0] = evaluated_integrand.real(); \
              result[1] = evaluated_integrand.imag(); \
              return 0; \
          } \
        } \
        static const int ncomp = 2; \
        struct userdata_t \
        { \
            const secdecutil::IntegrandContainer<complex_template<cubareal>, cubareal const * const> * integrand_container; \
            const cubareal& zero_border; \
        } typed_userdata{nullptr,zero_border}; \
        CUBA_STRUCT_BODY \
 \
        std::function<secdecutil::UncorrelatedDeviation<complex_template<cubareal>> \
          (const secdecutil::IntegrandContainer<complex_template<cubareal>, cubareal const * const>&)> \
          get_together_integrate() \
          { \
            return [this] (const secdecutil::IntegrandContainer<complex_template<cubareal>, cubareal const * const>& integrand_container) { \
              CUBA_INTEGRATE_BODY \
              integrand_container.process_errors(); \
              return secdecutil::UncorrelatedDeviation<complex_template<cubareal>>({integral.at(0),integral.at(1)},{error.at(0),error.at(1)}); \
            }; \
        }; \
      public: \
        int flags; \
        cubareal zero_border; \
        static constexpr bool cuda_compliant_integrator = false; \
      };
      COMPLEX_CUBAREAL_INTEGRATOR(std::complex)
      #ifdef SECDEC_WITH_CUDA
        COMPLEX_CUBAREAL_INTEGRATOR(thrust::complex)
      #endif
      #undef COMPLEX_CUBA_INTEGRATOR

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
            std::string statefiledir; \
            constexpr static void* spin = nullptr; \
            constexpr static int integrator_type = 1; \
            std::shared_ptr<long long int> neval; \
            char togethermode = '0'; \
            \
            Vegas \
            ( \
                cubareal epsrel = 1e-2, \
                cubareal epsabs = 1e-7, \
                int flags = 0, \
                int seed = 0, \
                long long int mineval = 0, \
                long long int maxeval = 1e6, \
                cubareal zero_border = 0., \
                long long int nstart = 10000, \
                long long int nincrease = 5000, \
                long long int nbatch = 1000, \
                std::shared_ptr<long long int> neval = std::make_shared<long long int>(0), \
                std::string statefiledir = "", \
                char togethermode = '0' \
            ) : \
                epsrel(epsrel),epsabs(epsabs), \
                seed(seed),mineval(mineval),maxeval(maxeval), \
                nstart(nstart),nincrease(nincrease),nbatch(nbatch),neval(neval),statefiledir(statefiledir),togethermode(togethermode) \
            { \
                this->flags = flags; \
                this->zero_border = zero_border; \
            }; \
            Vegas \
            ( \
             const Vegas& original \
            ) : \
                Vegas(original.epsrel,original.epsabs,original.flags, \
                original.seed,original.mineval,original.maxeval, \
                original.zero_border,original.nstart, \
                original.nincrease,original.nbatch,original.neval,original.statefiledir,original.togethermode) \
            { \
                this->copy_together_flag(original); \
            };



        #define VEGAS_CALL(HAVE_ZERO_BORDER) \
        ::llVegas \
        ( \
            this->ndim, \
            this->ncomp, \
            this->template cuba_integrand_prototype<HAVE_ZERO_BORDER>, \
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
            statefiledir.empty() ? "" : (statefiledir + '/' + togethermode).c_str(), \
            spin, \
            neval.get(), \
            &fail, \
            this->integral.data(), \
            this->error.data(), \
            this->prob.data() \
        )

        #define VEGAS_INTEGRATE_BODY \
        /* Cuba output values */ \
        int fail; \
        \
        /* Cuba call */ \
        if (this->zero_border == 0) /* no zero_border */ \
            VEGAS_CALL(false); \
        else \
            VEGAS_CALL(true); \
        if(togethermode != '0' ) \
            togethermode++;
        

      template <typename T>
      struct Vegas : CubaIntegrator<T> {
        VEGAS_STRUCT_BODY
        void call_cuba(){
          VEGAS_INTEGRATE_BODY
        };
      };

      #define COMPLEX_VEGAS(complex_template) \
      template <typename T> \
      struct Vegas<complex_template<T>> : CubaIntegrator<complex_template<T>> { \
        std::unique_ptr<Integrator<T,T>> get_real_integrator(){ \
          return std::unique_ptr<Integrator<T,T>>( \
                                                     new Vegas<T>( \
                                                                     epsrel,epsabs,this->flags, \
                                                                     seed,mineval,maxeval,this->zero_border, \
                                                                     nstart,nincrease,nbatch,neval,statefiledir,'1' \
                                                                 ) \
                                                 ); \
        }; \
        VEGAS_STRUCT_BODY \
        void call_cuba(){ \
          VEGAS_INTEGRATE_BODY \
        }; \
      };

      COMPLEX_VEGAS(std::complex)
      #ifdef SECDEC_WITH_CUDA
        COMPLEX_VEGAS(thrust::complex)
      #endif
      #undef COMPLEX_VEGAS

      #undef VEGAS_CALL
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
            std::string statefiledir; \
            constexpr static void* spin = nullptr; \
            constexpr static int integrator_type = 2; \
            std::shared_ptr<long long int> neval; \
            char togethermode = '0'; \
            \
            Suave \
            ( \
                cubareal epsrel = 1e-2, \
                cubareal epsabs = 1e-7, \
                int flags = 0, \
                int seed = 0, \
                long long int mineval = 0, \
                long long int maxeval = 1e6, \
                cubareal zero_border = 0., \
                long long int nnew = 1000, \
                long long int nmin = 10, \
                cubareal flatness = 25., \
                std::shared_ptr<long long int> neval = std::make_shared<long long int>(0), \
                std::string statefiledir = "", \
                char togethermode = '0' \
            ) : \
                epsrel(epsrel),epsabs(epsabs), \
                seed(seed),mineval(mineval),maxeval(maxeval), \
                nnew(nnew),nmin(nmin),flatness(flatness),neval(neval), \
                statefiledir(statefiledir),togethermode(togethermode) \
            { \
                this->flags = flags; \
                this->zero_border = zero_border; \
            }; \
            \
            Suave \
            ( \
                const Suave& original \
            ) : \
                Suave(original.epsrel,original.epsabs, \
                original.flags,original.seed,original.mineval, \
                original.maxeval,original.zero_border,original.nnew, \
                original.nmin,original.flatness,original.neval, \
                original.statefiledir,original.togethermode) \
            { \
                this->copy_together_flag(original); \
            };

        #define SUAVE_CALL(HAVE_ZERO_BORDER) \
        ::llSuave \
        ( \
            this->ndim, \
            this->ncomp, \
            this->template cuba_integrand_prototype<HAVE_ZERO_BORDER>, \
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
            statefiledir.empty() ? "" : (statefiledir + '/' + togethermode).c_str(), \
            spin, \
            &nregions, \
            neval.get(), \
            &fail, \
            this->integral.data(), \
            this->error.data(), \
            this->prob.data() \
        )

        #define SUAVE_INTEGRATE_BODY \
        /* Cuba output values */ \
        int fail; \
        int nregions; \
        \
        /* Cuba call */ \
        if (this->zero_border == 0) /* no zero_border */ \
            SUAVE_CALL(false); \
        else \
            SUAVE_CALL(true); \
        if(togethermode != '0' ) \
            togethermode++;

      template <typename T>
      struct Suave : CubaIntegrator<T> {
        SUAVE_STRUCT_BODY
        void call_cuba(){
          SUAVE_INTEGRATE_BODY
        }
      };

      #define COMPLEX_SUAVE(complex_template) \
      template <typename T> \
      struct Suave<complex_template<T>> : CubaIntegrator<complex_template<T>> { \
        std::unique_ptr<Integrator<T,T>> get_real_integrator(){ \
          return std::unique_ptr<Integrator<T,T>>( \
                                                     new Suave<T>( \
                                                                     epsrel,epsabs,this->flags, \
                                                                     seed,mineval,maxeval,this->zero_border, \
                                                                     nnew,nmin,flatness,neval,statefiledir,'1' \
                                                                 ) \
                                                 ); \
        }; \
        SUAVE_STRUCT_BODY \
        void call_cuba(){ \
          SUAVE_INTEGRATE_BODY \
        } \
      };

      COMPLEX_SUAVE(std::complex)
      #ifdef SECDEC_WITH_CUDA
        COMPLEX_SUAVE(thrust::complex)
      #endif
      #undef COMPLEX_SUAVE

      #undef SUAVE_CALL
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
            std::string statefiledir; \
            constexpr static void * spin = nullptr; \
            constexpr static int integrator_type = 3; \
            std::shared_ptr<long long int> neval; \
            char togethermode = '0'; \
            \
            Divonne \
            ( \
                cubareal epsrel = 1e-2, \
                cubareal epsabs = 1e-7, \
                int flags = 0, \
                int seed = 0, \
                long long int mineval = 0, \
                long long int maxeval = 1e6, \
                cubareal zero_border = 0., \
                int key1 = 2000, \
                int key2 = 1, \
                int key3 = 1, \
                int maxpass = 4, \
                cubareal border = 0., \
                cubareal maxchisq = 1., \
                cubareal mindeviation = .15, \
                std::shared_ptr<long long int> neval = std::make_shared<long long int>(0), \
                std::string statefiledir = "", \
                char togethermode = '0' \
            ) : \
                epsrel(epsrel),epsabs(epsabs), \
                seed(seed),mineval(mineval),maxeval(maxeval), \
                key1(key1), key2(key2), key3(key3), maxpass(maxpass), \
                border(border), maxchisq(maxchisq), mindeviation(mindeviation),neval(neval), \
                statefiledir(statefiledir),togethermode(togethermode) \
            { \
                this->flags = flags; \
                this->zero_border = zero_border; \
            }; \
            \
            Divonne \
            ( \
                const Divonne& original \
            ) : \
                Divonne(original.epsrel,original.epsabs,original.flags,\
                original.seed,original.mineval,original.maxeval, \
                original.zero_border,original.key1,original.key2, \
                original.key3,original.maxpass,original.border, \
                original.maxchisq,original.mindeviation,original.neval, \
                original.statefiledir,original.togethermode) \
            { \
                this->copy_together_flag(original); \
            };

        #define DIVNONNE_CALL(HAVE_ZERO_BORDER) \
        ::llDivonne \
        ( \
            this->ndim, \
            this->ncomp, \
            this->template cuba_integrand_prototype<HAVE_ZERO_BORDER>, \
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
            statefiledir.empty() ? "" : (statefiledir + '/' + togethermode).c_str(), \
            spin, \
            &nregions, \
            neval.get(), \
            &fail, \
            this->integral.data(), \
            this->error.data(), \
            this->prob.data() \
        )

        #define DIVONNE_INTEGRATE_BODY \
        /* Cuba output values */ \
        int nregions; \
        int fail; \
        \
        /* Cuba call */ \
        if (this->zero_border == 0) /* no zero_border */ \
            DIVNONNE_CALL(false); \
        else \
            DIVNONNE_CALL(true); \
        if(togethermode != '0' ) \
            togethermode++;

      template <typename T>
      struct Divonne : CubaIntegrator<T> {
        DIVONNE_STRUCT_BODY
        void call_cuba(){
          DIVONNE_INTEGRATE_BODY
        }
      };

      #define COMPLEX_DIVONNE(complex_template) \
      template <typename T> \
      struct Divonne<complex_template<T>> : CubaIntegrator<complex_template<T>> { \
        std::unique_ptr<Integrator<T,T>> get_real_integrator(){ \
          return std::unique_ptr<Integrator<T,T>>( \
                                                     new Divonne<T>( \
                                                                       epsrel,epsabs,this->flags, \
                                                                       seed,mineval,maxeval,this->zero_border, \
                                                                       key1, key2, key3, maxpass, \
                                                                       border, maxchisq, mindeviation,neval,statefiledir,'1' \
                                                                   ) \
                                                 ); \
        }; \
        DIVONNE_STRUCT_BODY \
        void call_cuba(){ \
          DIVONNE_INTEGRATE_BODY \
        } \
      };

      COMPLEX_DIVONNE(std::complex)
      #ifdef SECDEC_WITH_CUDA
        COMPLEX_DIVONNE(thrust::complex)
      #endif
      #undef COMPLEX_DIVONNE

      #undef DIVONNE_CALL
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
            std::string statefiledir; \
            constexpr static void* spin = nullptr; \
            constexpr static int integrator_type = 4; \
            std::shared_ptr<long long int> neval; \
            char togethermode = '0'; \
            \
            Cuhre \
            ( \
                cubareal epsrel = 1e-2, \
                cubareal epsabs = 1e-7, \
                int flags = 0, \
                long long int mineval = 0, \
                long long int maxeval = 1e6, \
                cubareal zero_border = 0., \
                int key = 0, \
                std::shared_ptr<long long int> neval = std::make_shared<long long int>(0), \
                std::string statefiledir = "", \
                char togethermode = '0' \
            ) : \
                epsrel(epsrel),epsabs(epsabs), \
                mineval(mineval),maxeval(maxeval),key(key),neval(neval), \
                statefiledir(statefiledir),togethermode(togethermode) \
            { \
                this->flags = flags; \
                this->zero_border = zero_border; \
            }; \
            \
            Cuhre \
            ( \
                const Cuhre& original \
            ) : \
                Cuhre(original.epsrel,original.epsabs,original.flags, \
                original.mineval,original.maxeval,original.zero_border,original.key,original.neval, \
                original.statefiledir,original.togethermode) \
            { \
                this->copy_together_flag(original); \
            };

        #define CUHRE_CALL(HAVE_ZERO_BORDER) \
        ::llCuhre \
        ( \
            this->ndim, \
            this->ncomp, \
            this->template cuba_integrand_prototype<HAVE_ZERO_BORDER>, \
            this->userdata, \
            nvec, \
            epsrel, \
            epsabs, \
            this->flags, \
            mineval, \
            maxeval, \
            key, \
            statefiledir.empty() ? "" : (statefiledir + '/' + togethermode).c_str(), \
            spin, \
            &nregions, \
            neval.get(), \
            &fail, \
            this->integral.data(), \
            this->error.data(), \
            this->prob.data() \
        )

        #define CUHRE_INTEGRATE_BODY \
        /* Cuba output values */ \
        int nregions; \
        int fail; \
        \
        /* Cuba call */ \
        if (this->zero_border == 0) /* no zero_border */ \
            CUHRE_CALL(false); \
        else \
            CUHRE_CALL(true); \
        if(togethermode != '0' ) \
            togethermode++;

      template <typename T>
      struct Cuhre : CubaIntegrator<T> {
        CUHRE_STRUCT_BODY
        void call_cuba(){
          CUHRE_INTEGRATE_BODY
        }
      };

      #define COMPLEX_CUHRE(complex_template) \
      template <typename T> \
      struct Cuhre<complex_template<T>> : CubaIntegrator<complex_template<T>> { \
        std::unique_ptr<Integrator<T,T>> get_real_integrator(){ \
          return std::unique_ptr<Integrator<T,T>>( \
                                                     new Cuhre<T>( \
                                                                     epsrel,epsabs,this->flags, \
                                                                     mineval,maxeval,this->zero_border, \
                                                                     key,neval,statefiledir,'1' \
                                                                 ) \
                                                 ); \
        }; \
        CUHRE_STRUCT_BODY \
        void call_cuba(){ \
          CUHRE_INTEGRATE_BODY \
        } \
      };

      COMPLEX_CUHRE(std::complex)
      #ifdef SECDEC_WITH_CUDA
        COMPLEX_CUHRE(thrust::complex)
      #endif
      #undef COMPLEX_CUHRE

      #undef CUHRE_CALL
      #undef CUHRE_STRUCT_BODY
      #undef CUHRE_INTEGRATE_BODY
    }
}

#endif

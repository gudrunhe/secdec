#ifndef SecDecUtil_cquad_hpp_included
#define SecDecUtil_cquad_hpp_included

/*
 * This file implements a convenience wrapper around the
 * cquad integrator form the gnu scientific library (gsl).
 */

#include <complex>
#include <gsl/gsl_integration.h> // gsl_integration_cquad
#include <gsl/gsl_errno.h> // gsl_strerror, GSL_SUCCESS
#include <iostream>
#include <stdexcept>
#include <string>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/integrators/integrator.hpp>
#include <secdecutil/uncertainties.hpp>


namespace secdecutil
{

   namespace gsl
   {

      // this error is thrown if an error with the gsl occurs
      struct gsl_error : public std::runtime_error { using std::runtime_error::runtime_error; };

      void custom_gsl_error_handler(const char * reason, const char * file, int line, int gsl_errno)
      {
          throw gsl_error( std::string(gsl_strerror(gsl_errno)) + std::string(": ") + std::string(reason));
      };


     /*
      * C++ wrapper class for the gsl's cquad integrator.
      *
      * cquad always uses the float type "double", which can
      * be different from the type "T" used by the IntegrandContainer.
      */

      // real version for general type
      template<typename T>
      struct CQuad : Integrator<T,T>
      {
      protected:
        static double integrand_prototype_for_gsl(double x, void * params)
        {
          auto& integrand_container = *( reinterpret_cast<const secdecutil::IntegrandContainer<T, T const * const> *>(params) );

          // "x" has type double, but the integrand expects type T --> cast to type T
          T integration_variable = x;

          // initialize result with NaN --> result will be NaN if integrand throws an error
          double evaluated_integrand = std::nan("");

          evaluated_integrand = integrand_container.integrand(&integration_variable);

          return evaluated_integrand;
        };
        gsl_function integrand_for_gsl{/* function */ integrand_prototype_for_gsl, /* params --> set in "get_integrate" */ nullptr};
        int ndim;
        const double a = 0.0;
        const double b = 1.0;
        double value;
        double abserr;
        size_t nevals;
        std::shared_ptr<gsl_integration_cquad_workspace> workspace;

      public:
        double epsrel;
        double epsabs;
        const size_t n;
        bool verbose;

        const std::shared_ptr<gsl_integration_cquad_workspace>& get_workspace() const
        {
            return workspace;
        };

        // Constructors
        CQuad
        (
            double epsrel = 1e-2,
            double epsabs = 1e-7,
            size_t n = 100, // number of intervals to be kept simultaneously
            bool verbose = false
        ) :
            epsrel(epsrel),epsabs(epsabs),n(n),verbose(verbose)
        {
            gsl_error_handler_t * original_error_handler = gsl_set_error_handler_off();
            gsl_set_error_handler(custom_gsl_error_handler);
            workspace.reset( gsl_integration_cquad_workspace_alloc(n) , gsl_integration_cquad_workspace_free );
            gsl_set_error_handler(original_error_handler);
        };

        // construct with allocated workspace
        CQuad
        (
            double epsrel,
            double epsabs,
            size_t n, // number of intervals to be kept simultaneously
            bool verbose,
            const std::shared_ptr<gsl_integration_cquad_workspace>& workspace
        ) :
            epsrel(epsrel),epsabs(epsabs),n(n),verbose(verbose),workspace(workspace)
        {};

        // copy Constructor
        CQuad
        (
            const CQuad& original
        ) :
            CQuad(original.epsrel,original.epsabs,original.n,original.verbose)
        {};


      protected:
        std::function<secdecutil::UncorrelatedDeviation<T>
          (const secdecutil::IntegrandContainer<T, T const * const>&)> get_integrate()
        {
          return [this] (const secdecutil::IntegrandContainer<T, T const * const>& integrand_container)
            {
              ndim = integrand_container.number_of_integration_variables;

              // can only have one integration variable
              if (ndim > 1)
                  throw std::invalid_argument("\"CQuad\" can only be used for one dimensional integrands (got ndim=" + std::to_string(ndim) + ").");

              if (verbose)
              {
                  std::cout << "CQuad input parameters:" << std::endl;
                  std::cout << "  epsrel " << epsrel << std::endl;
                  std::cout << "  epsabs " << epsabs << std::endl;
                  std::cout << "  n " << n << std::endl;
                  std::cout << std::endl;
              }

              // nasty removal of constness --> restored in gsl_integrand_prototype
              integrand_for_gsl.params = const_cast<void*>( reinterpret_cast<const void*>(&integrand_container) );

              // call the cquad routine from the gsl
              gsl_error_handler_t * original_error_handler = gsl_set_error_handler_off();
              gsl_set_error_handler(custom_gsl_error_handler);
              gsl_integration_cquad(&integrand_for_gsl, a, b, epsabs, epsrel, workspace.get(), &value, &abserr, &nevals);
              gsl_set_error_handler(original_error_handler);

              // pack result
              secdecutil::UncorrelatedDeviation<T> result(value,abserr);

              if (verbose)
              {
                  std::cout << "CQuad result:" << std::endl;
                  std::cout << "  " << result << std::endl << std::endl;
              }

              return result;
            };
        };

      };

      // complex version
      template<typename T>
      struct CQuad<std::complex<T>> : Integrator<std::complex<T>,T>, CQuad<T>
      {
        public:
          std::unique_ptr<Integrator<T,T>> get_real_integrator()
          {
              return std::unique_ptr<Integrator<T,T>>( new CQuad<T>(this->epsrel,this->epsabs,this->n,this->verbose,this->workspace) );
          };

          // Constructor
          CQuad
          (
              double epsrel = 1e-2,
              double epsabs = 1e-7,
              size_t n = 100, // number of intervals to be kept simultaneously
              bool verbose = false
          ) :
              CQuad<T>(epsrel,epsabs,n,verbose)
          {};

      };

   }

}

#endif

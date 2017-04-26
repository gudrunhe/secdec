#ifndef SecDecUtil_qmc_hpp_included
#define SecDecUtil_qmc_hpp_included

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

        template<typename return_t>
        struct Qmc : Integrator<return_t,return_t>
        {
        protected:

            using input_t = return_t;
            std::function<secdecutil::UncorrelatedDeviation<return_t>(const secdecutil::IntegrandContainer<return_t, input_t const * const>&)>
            get_integrate()
            {
                std::function<secdecutil::UncorrelatedDeviation<return_t>(const secdecutil::IntegrandContainer<return_t, input_t const * const>&)> integrate_function = [this] (const secdecutil::IntegrandContainer<return_t, input_t const * const>& integrand_container)
                {
                    ::integrators::result<return_t> result;
                    result = integrator.integrate(integrand_container.integrand,integrand_container.number_of_integration_variables);
                    return secdecutil::UncorrelatedDeviation<return_t> { result.integral, result.error };
                };
                return integrate_function;
            };
        public:

            ::integrators::Qmc<return_t,input_t> integrator;

        };

        template<typename return_t>
        struct Qmc<std::complex<return_t>> : Integrator<std::complex<return_t>, return_t>
        {
        protected:

            using input_t = return_t;
            std::unique_ptr<Integrator<return_t, input_t>> get_real_integrator()
            {
                throw std::runtime_error("Separate integration of real and imaginary part has no advantage for this non-adaptive integrator. Please use \"together = true\".");
            };

            std::function<secdecutil::UncorrelatedDeviation<std::complex<return_t>>
            (const secdecutil::IntegrandContainer<std::complex<return_t>, input_t const * const>&)>
            get_together_integrate()
            {
                std::function<secdecutil::UncorrelatedDeviation<std::complex<return_t>>(const secdecutil::IntegrandContainer<std::complex<return_t>, input_t const * const>&)> integrate_function = [this] (const secdecutil::IntegrandContainer<std::complex<return_t>, input_t const * const>& integrand_container)
                {
                    ::integrators::result<std::complex<return_t>> result;
                    result = integrator.integrate(integrand_container.integrand,integrand_container.number_of_integration_variables);
                    return secdecutil::UncorrelatedDeviation<std::complex<return_t>> { result.integral, result.error };
                };
                return integrate_function;
            };
            
        public:

            ::integrators::Qmc<std::complex<return_t>,input_t> integrator;

            Qmc() { this->together = true;  };
            
        };

    }

}

#endif

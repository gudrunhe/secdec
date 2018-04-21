#ifndef SecDecUtil_integrator_hpp_included
#define SecDecUtil_integrator_hpp_included

/*
 * This file defines a general integrator interface for real and complex-valued functions.
 *
 * Real integrators only have to implement the function "get_integrate()",
 * which returns a function taking an integrand container and returning an UncorrelatedDeviation.
 *
 * Complex integrators have to implement the function "get_together_integrate()"
 * to integrate real and imaginary part at the same time
 * and/or "get_real_integrator()", which should return a unique pointer to a real-valued version of
 * the integrator. The latter can then be used to integrate real and imaginaty part separately
 * if the boolean member variable "together" is set to false.
 */

#include <complex>
#include <memory>
#include <stdexcept>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/uncertainties.hpp>

namespace secdecutil
{

  template<typename return_t, typename input_t>
  struct Integrator
  {
  protected:
    virtual std::function<secdecutil::UncorrelatedDeviation<return_t>
      (const secdecutil::IntegrandContainer<return_t, input_t const * const>&)>
      get_integrate() = 0;

  public:
    const std::function<secdecutil::UncorrelatedDeviation<return_t>
      (const secdecutil::IntegrandContainer<return_t, input_t const * const>&)>
      integrate;

    Integrator() :
    integrate
    (
      [ this ] (const secdecutil::IntegrandContainer<return_t, input_t const * const>& integrand_container)
      {
          return get_integrate()(integrand_container);
      }
    ) {};

    virtual ~Integrator() = default;

  };

  template<typename return_t, typename input_t>
  struct Integrator<std::complex<return_t>, input_t>
  {
  protected:
    virtual std::unique_ptr<Integrator<return_t, input_t>> get_real_integrator()
    {
      throw std::runtime_error("Separate integration of real and imaginary part is not available because pointer to real-valued integrator is not implemented for this integrator. Try \"together = true\".");
    }

    virtual std::function<secdecutil::UncorrelatedDeviation<std::complex<return_t>>
      (const secdecutil::IntegrandContainer<std::complex<return_t>, input_t const * const>&)>
      get_together_integrate()
    {
      throw std::runtime_error("Simultaneous integration of real and imaginary part is not implemented for this integrator. Try \"together = false\".");
    }

  public:

    bool together;
    const std::function<secdecutil::UncorrelatedDeviation<std::complex<return_t>> (const secdecutil::IntegrandContainer<std::complex<return_t>, input_t const * const>&)> integrate;

    Integrator() :
    together(false),
    integrate
    (
      [ this ] (const secdecutil::IntegrandContainer<std::complex<return_t>, input_t const * const>& integrand_container)
      {
        if (together) {

          return get_together_integrate()(integrand_container);

        } else {

          auto real_integrator = this->get_real_integrator();
          auto real_part = real_integrator->integrate(complex_to_real(integrand_container,std::real));
          auto imag_part = real_integrator->integrate(complex_to_real(integrand_container,std::imag));
          return secdecutil::UncorrelatedDeviation<std::complex<return_t>>({real_part.value,imag_part.value},{real_part.uncertainty,imag_part.uncertainty});

        }
    }
    ) {};

    virtual ~Integrator() = default;

  };

/*
 * The class "MultiIntegrator" defines an integrator
 * that switches between two integrators depending on
 * the dimension of the integrand: If the integrand
 * is lower dimensional than "critical_dim" then
 * "low_dim_integrator", otherwise "high_dim_integrator"
 * is used.
 */

  template<typename return_t, typename input_t>
  struct MultiIntegrator : Integrator<return_t,input_t>
  {
    Integrator<return_t,input_t>& low_dim_integrator;
    Integrator<return_t,input_t>& high_dim_integrator;
    int critical_dim;

    // Constructor
    MultiIntegrator
    (
      Integrator<return_t,input_t>& low_dim_integrator,
      Integrator<return_t,input_t>& high_dim_integrator,
      int critical_dim
    ) :
      low_dim_integrator(low_dim_integrator),
      high_dim_integrator(high_dim_integrator),
      critical_dim(critical_dim)
    {};

    std::function<secdecutil::UncorrelatedDeviation<return_t>
      (const secdecutil::IntegrandContainer<return_t, input_t const * const>&)>
      get_integrate()
      {
        return [this](const secdecutil::IntegrandContainer<return_t, input_t const * const>& ic)
          {
            if (ic.number_of_integration_variables < critical_dim)
              return low_dim_integrator.integrate(ic);
            else
              return high_dim_integrator.integrate(ic);
          };
      }

  };

  template<typename return_t, typename input_t>
  struct MultiIntegrator<std::complex<return_t>,input_t> : Integrator<std::complex<return_t>,input_t>
  {
    Integrator<std::complex<return_t>,input_t>& low_dim_integrator;
    Integrator<std::complex<return_t>,input_t>& high_dim_integrator;
    int critical_dim;

    std::function<secdecutil::UncorrelatedDeviation<std::complex<return_t>>
      (const secdecutil::IntegrandContainer<std::complex<return_t>, input_t const * const>&)>
      get_together_integrate()
      {
        return [this] (const secdecutil::IntegrandContainer<std::complex<return_t>, input_t const * const>& ic)
        {
          if (ic.number_of_integration_variables < critical_dim)
            return low_dim_integrator.integrate(ic);
          else
            return high_dim_integrator.integrate(ic);
        };
    };

    // Constructor
    MultiIntegrator
    (
      Integrator<std::complex<return_t>,input_t>& low_dim_integrator,
      Integrator<std::complex<return_t>,input_t>& high_dim_integrator,
      int critical_dim
    ) :
      low_dim_integrator(low_dim_integrator),
      high_dim_integrator(high_dim_integrator),
      critical_dim(critical_dim)
    { this->together = true; };

  };

}

#endif

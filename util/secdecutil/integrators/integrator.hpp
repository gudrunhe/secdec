#ifndef SecDecUtil_integrator_hpp_included
#define SecDecUtil_integrator_hpp_included

/*
 * This file defines a general integrator interface for real and complex-valued functions.
 *
 * Real integrators only have to implement the function integrate(),
 * which returns a function taking an integrand container and returning an UncorrelatedDeviation.
 * 
 * Complex integrators have to implement the function integrate_together()
 * to integrate real and imaginary part at the same time
 * and/or get_real_integrator(), which should return a pointer to a real-valued version of the integrator.
 * The latter can then be used to integrate real and imaginaty part separately
 * if the boolean together is set to false.
 */

#include <complex>
#include <memory>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/uncertainties.hpp>

namespace secdecutil
{

  template<typename return_t, typename input_t>
  struct Integrator
  {
    virtual std::function<secdecutil::UncorrelatedDeviation<return_t>
      (const secdecutil::IntegrandContainer<return_t, input_t const * const>&)>
      integrate() = 0;
  };
	
  template<typename return_t, typename input_t>
  struct Integrator<std::complex<return_t>, input_t>
  {
  protected:
    virtual std::shared_ptr<Integrator<return_t, input_t>> get_real_integrator(){
      throw std::runtime_error("Separate integration of real and imaginary parts not available because pointer to real-valued integrator is not implemented.");
    }

    virtual std::function<secdecutil::UncorrelatedDeviation<std::complex<return_t>>
      (const secdecutil::IntegrandContainer<std::complex<return_t>, input_t const * const>&)>
      integrate_together()
    {
      throw std::runtime_error("Simultaneous integration of real and imaginary parts not implemented.");
    }
	
  public:	
    bool together = true;
	
    std::function<secdecutil::UncorrelatedDeviation<std::complex<return_t>>
      (const secdecutil::IntegrandContainer<std::complex<return_t>, input_t const * const>&)>
      integrate() {
      return [ this ]
	(const secdecutil::IntegrandContainer<std::complex<return_t>, input_t const * const>& integrand_container)
	{
	  if (together)
	    {
	      return integrate_together()(integrand_container);
	    }
	  else
	    {
	      auto real_integrator = this->get_real_integrator();
	      auto real_part = real_integrator->integrate()(complex_to_real(integrand_container,std::real));
	      auto imag_part = real_integrator->integrate()(complex_to_real(integrand_container,std::imag));
	      return secdecutil::UncorrelatedDeviation<std::complex<return_t>>({real_part.value,imag_part.value},{real_part.uncertainty,imag_part.uncertainty});
	    }
	};
    };
  };

}

#endif

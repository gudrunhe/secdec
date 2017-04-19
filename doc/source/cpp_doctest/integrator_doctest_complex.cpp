#include <iostream>
#include <complex>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/uncertainties.hpp>
#include <secdecutil/integrators/cuba.hpp>

int main()
{
    using input_t = const double * const;
    using return_t = std::complex<double>;

    secdecutil::cuba::Vegas<return_t> integrator;
    std::function<return_t(input_t)> f = [] (input_t x) { return return_t{x[0],x[1]}; };
    secdecutil::IntegrandContainer<return_t,input_t> c(2,f);

    integrator.together = false; // integrate real and imaginary part separately (default)
    secdecutil::UncorrelatedDeviation<return_t> result_separate = integrator.integrate(c);

    integrator.together = true; // integrate real and imaginary part simultaneously
    secdecutil::UncorrelatedDeviation<return_t> result_together = integrator.integrate(c);

    std::cout << "result_separate: " << result_separate << std::endl;
    std::cout << "result_together: " << result_together << std::endl;
    
}

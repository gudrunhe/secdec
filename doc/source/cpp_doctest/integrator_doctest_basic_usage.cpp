#include <iostream>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/uncertainties.hpp>
#include <secdecutil/integrators/cuba.hpp>

int main()
{
    using input_t = const double * const;
    using return_t = double;

    secdecutil::cuba::Vegas<return_t> integrator;
    integrator.epsrel = 1e-4;
    integrator.maxeval = 1e7;

    secdecutil::IntegrandContainer<return_t,input_t> c(2, [] (input_t x) { return x[0]*x[1]; });
    secdecutil::UncorrelatedDeviation<return_t> result = integrator.integrate(c);

    std::cout << "result: " << result << std::endl;
}

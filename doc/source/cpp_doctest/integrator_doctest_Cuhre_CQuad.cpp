#include <iostream>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/uncertainties.hpp>
#include <secdecutil/integrators/integrator.hpp>
#include <secdecutil/integrators/cuba.hpp>
#include <secdecutil/integrators/cquad.hpp>

int main()
{
    using input_base_t = double;
    using input_t = const input_base_t * const;
    using return_t = double;

    secdecutil::cuba::Cuhre<return_t> cuhre;
    cuhre.epsrel = 1e-4;
    cuhre.maxeval = 1e7;

    secdecutil::CQuad<return_t> cquad;
    cquad.epsrel = 1e-10;
    cquad.epsabs = 1e-13;

    secdecutil::MultiIntegrator<return_t,input_base_t> integrator(cquad,cuhre,2);

    secdecutil::IntegrandContainer<return_t,input_t> one_dimensional(1, [] (input_t x) { return x[0]; });
    secdecutil::IntegrandContainer<return_t,input_t> two_dimensional(2, [] (input_t x) { return x[0]*x[1]; });

    secdecutil::UncorrelatedDeviation<return_t> result_1d = integrator.integrate(one_dimensional); // uses cquad
    secdecutil::UncorrelatedDeviation<return_t> result_2d = integrator.integrate(two_dimensional); // uses cuhre

    std::cout << "result_1d: " << result_1d << std::endl;
    std::cout << "result_2d: " << result_2d << std::endl;
}

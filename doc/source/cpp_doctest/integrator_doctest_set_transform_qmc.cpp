#include <iostream>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/uncertainties.hpp>
#include <secdecutil/integrators/qmc.hpp>

using input_base_t = double;
using input_t = input_base_t const * const;
using return_t = double;
using container_t = secdecutil::IntegrandContainer<return_t,input_t>;
using result_t = secdecutil::UncorrelatedDeviation<return_t>;

const int seed = 12345, maxdim = 4;

int main()
{
    /*
     * minimal instantiation
     */
    secdecutil::integrators::Qmc
    <
        return_t, // the return type of the integrand
        maxdim, // the highest dimension this integrator will be used for
        ::integrators::transforms::Baker::type // the integral transform
    > integrator_baker;
    integrator_baker.randomgenerator.seed(seed);

    /*
     * disable adaptation
     */
    secdecutil::integrators::Qmc
    <
        return_t, // the return type of the integrand
        maxdim, // the highest dimension this integrator will be used for
        ::integrators::transforms::Korobov<4,1>::type, // the integral transform
        container_t, // the functor type to be passed to this integrator
        ::integrators::fitfunctions::None::type // the fit function
    > integrator_korobov4x1;
    integrator_korobov4x1.randomgenerator.seed(seed);

    /*
     * enable adaptation
     */
    secdecutil::integrators::Qmc
    <
        return_t, // the return type of the integrand
        maxdim, // the highest dimension this integrator will be used for
        ::integrators::transforms::Sidi<3>::type, // the integral transform
        container_t, // the functor type to be passed to this integrator
        ::integrators::fitfunctions::PolySingular::type // the fit function
    > integrator_sidi3_adaptive;
    integrator_sidi3_adaptive.randomgenerator.seed(seed);

    // define the integrand as a functor
    container_t integrand(
                             4, // dimension
                             [] (input_t x, secdecutil::ResultInfo* result_info) { return x[0]*x[1]*x[2]*x[3]; } // integrand function
                         );

    // compute the integral with different settings
    result_t result_baker = integrator_baker.integrate(integrand);
    result_t result_korobov4x1 = integrator_korobov4x1.integrate(integrand);
    result_t result_sidi3_adaptive = integrator_sidi3_adaptive.integrate(integrand);

    // print the results
    std::cout << "baker: " << result_baker << std::endl;
    std::cout << "Korobov (weights 4, 1): " << result_korobov4x1 << std::endl;
    std::cout << "Sidi (weight 3, adaptive): " << result_sidi3_adaptive << std::endl;
}

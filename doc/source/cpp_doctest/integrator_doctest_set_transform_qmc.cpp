#include <iostream>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/uncertainties.hpp>
#include <secdecutil/integrators/qmc.hpp>

int main()
{
    using input_base_t = double;
    using input_t = input_base_t const * const;
    using return_t = double;
    using container_t = secdecutil::IntegrandContainer<return_t,input_t>;
    using result_t = secdecutil::UncorrelatedDeviation<return_t>;

    int seed = 12345;

    secdecutil::integrators::Qmc
    <
        return_t,
        container_t,
        ::integrators::transforms::Tent<input_base_t> // custom integral transform
    > integrator_tent;
    integrator_tent.randomgenerator.seed(seed);

    secdecutil::integrators::Qmc <return_t> integrator_default; // uses the default integral transform
    integrator_default.randomgenerator.seed(seed);

    container_t integrand(4, [] (input_t x) { return x[0]*x[1]*x[2]*x[3]; });

    result_t result_default = integrator_default.integrate(integrand);
    result_t result_tent = integrator_tent.integrate(integrand);

    std::cout << "result with default transform: " << result_default << std::endl;
    std::cout << "result with tent transform: " << result_tent << std::endl;
}

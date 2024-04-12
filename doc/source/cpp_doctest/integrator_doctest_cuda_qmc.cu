#include <iostream>
#include <secdecutil/integrand_container.hpp>
#include <secdecutil/uncertainties.hpp>
#include <secdecutil/integrators/qmc.hpp>

using input_base_t = double;
using input_t = input_base_t const * const;
using return_t = double;
using container_t = secdecutil::IntegrandContainer<return_t,input_t>;
using result_t = secdecutil::UncorrelatedDeviation<return_t>;

/*
 * `container_t` cannot be used on the GPU --> define a different container type
 */
struct cuda_integrand_t
{
    const static unsigned number_of_integration_variables = 4;

    // integrand function
    #ifdef __CUDACC__
        __host__ __device__
    #endif
    return_t operator()(input_t x)
    {
        return x[0]*x[1]*x[2]*x[3];
    };

    void process_errors() const{ /* error handling */}
} cuda_integrand;

const int seed = 12345, maxdim = 4;

int main()
{
    /*
     * Qmc capable of sampling on the GPU
     */
    secdecutil::integrators::Qmc
    <
        return_t, // the return type of the integrand
        maxdim, // the highest dimension this integrator will be used for
        ::integrators::transforms::Sidi<3>::type, // the integral transform
        cuda_integrand_t, // the functor type to be passed to this integrator
        ::integrators::fitfunctions::PolySingular::type // the fit function (optional)
    > integrator_sidi3_adaptive_gpu;
    integrator_sidi3_adaptive_gpu.randomgenerator.seed(seed);

    // compute the integral with different settings
    result_t result_sidi3_adaptive_gpu = integrator_sidi3_adaptive_gpu.integrate(cuda_integrand);

    // print the results
    std::cout << "Sidi (weight 3, adaptive): " << result_sidi3_adaptive_gpu << std::endl;
}

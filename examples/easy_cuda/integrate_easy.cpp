#include <cmath> // std::log
#include <iostream> // std::cout
#include <numeric> // std::accumulate
#include <vector> // std::vector

#include <secdecutil/integrators/qmc.hpp> // secdecutil::integrators::Qmc
#include <secdecutil/series.hpp> // secdecutil::Series
#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation
#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply

#include <easy.hpp> // the namespace "easy"

int main()
{
    // Generate the integrands
    const std::vector<easy::nested_series_t<easy::cuda_integrand_t>> sector_integrands = easy::make_cuda_integrands({}, {});

    // Add integrands of sectors (together flag)
    const easy::nested_series_t<easy::cuda_together_integrand_t> all_sectors =
        std::accumulate(++sector_integrands.begin(), sector_integrands.end(), easy::cuda_together_integrand_t()+*sector_integrands.begin());

    // Integrate
    secdecutil::integrators::Qmc<easy::integrand_return_t,easy::maximal_number_of_integration_variables,integrators::transforms::Korobov<3>::type,easy::cuda_together_integrand_t> integrator;
    const easy::nested_series_t<secdecutil::UncorrelatedDeviation<easy::integrand_return_t>> result_without_prefactor =
        secdecutil::deep_apply( all_sectors, integrator.integrate );

    // Multiply by the prefactor
    const easy::nested_series_t<secdecutil::UncorrelatedDeviation<easy::integrand_return_t>> full_result = result_without_prefactor * easy::prefactor({}, {});

    // Print result
    std::cout.precision(15);
    std::cout << "Numerical Result: " << full_result << std::endl;
    std::cout << "Analytic Result: " << "(" << 1 << ")*eps^-1 + (" << 1.0-std::log(2.0) << + ") + O(eps)" << std::endl;
}

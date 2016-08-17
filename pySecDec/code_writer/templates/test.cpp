#include <iostream>
#include <stdexcept>
#include <vector>
#include <cstdint>
#include <cuba.h>
#include <numeric> // std::accumulate
#include <functional> // std::bind
#include <type_traits> // std::remove_const
#include <typeinfo>

#include <secdecutil/series.hpp> // Series
#include <secdecutil/uncertainties.hpp> // GaussianUncertainty
#include <secdecutil/sector_container.hpp> // SectorContainer to IntegrandContainer
#include <secdecutil/integrand_container.hpp> // IntegrandContainer
#include <secdecutil/deep_apply.hpp> // deep_apply

#include "%(name)s.hpp"

void print_integral_info()
{
    std::cout << "-- print_integral_info --" << std::endl;
    std::cout << "%(name)s::number_of_sectors " << %(name)s::number_of_sectors << std::endl;
    std::cout << "%(name)s::number_of_regulators " << %(name)s::number_of_regulators << std::endl;
    std::cout << "%(name)s::number_of_real_parameters " << %(name)s::number_of_real_parameters << std::endl;
    std::cout << "%(name)s::number_of_complex_parameters " << %(name)s::number_of_complex_parameters << std::endl;

    std::cout << "%(name)s::lowest_orders";
    for ( const auto& lowest_order : %(name)s::lowest_orders )
        std::cout << " " << lowest_order;
    std::cout << std::endl;

    std::cout << "%(name)s::highest_orders";
    for ( const auto& highest_order : %(name)s::highest_orders )
        std::cout << " " << highest_order;
    std::cout << std::endl;

    std::cout << "%(name)s::requested_orders";
    for ( const auto& requested_order : %(name)s::requested_orders )
        std::cout << " " << requested_order;
    std::cout << std::endl;
}

// TODO - userdata should be just nest.integrand not nest
int cuba_integrand_prototype(const int *ndim, const cubareal integration_variables[], const int *ncomp, cubareal result[], void *userdata)
{
    auto& nest = *( reinterpret_cast<secdecutil::IntegrandContainer<%(name)s::integrand_return_t, %(name)s::real_t const * const> *>(userdata) );
    result[0] = nest.integrand(integration_variables).real(); // TODO: omit ".real()" if return type is real
    return 0;
};

template<typename integrand_return_t, typename real_t, typename complex_t>
std::function<secdecutil::GaussianUncertainty<cubareal>(secdecutil::IntegrandContainer<integrand_return_t, real_t const * const>)>
cuba_integrate()
{
    return [ ] (secdecutil::IntegrandContainer<integrand_return_t, real_t const * const> nest)
    {
        std::cout << "-- Integrating --" << std::endl;

        struct cuba_parameters_t {
            int ncomp = 1;
            int nvec = 1;
            cubareal epsrel = 1e-2;
            cubareal epsabs = 1e-16;
            int verbose = 0;
            int seed = 0;
            int mineval = 0;
            int maxeval = 50000;
            int nstart = 1000;
            int nincrease = 500;
            int nbatch = 1000;
            int gridno = 0;
            const char * statefile = nullptr;
            void* spin = nullptr;
        } cuba_parameters;

        // Cuba output values
        std::array<cubareal, 1 > integral; // ncomp
        std::array<cubareal, 1 > error; // ncomp
        std::array<cubareal, 1 > prob; // ncomp
        int comp, nregions, neval, fail;

        // Cuba call
        // TODO - userdata should be just nest.integrand not nest
        Vegas(
              nest.number_of_integration_variables,
              cuba_parameters.ncomp,
              cuba_integrand_prototype,
              reinterpret_cast<void*>(&nest), // userdata
              cuba_parameters.nvec,
              cuba_parameters.epsrel,
              cuba_parameters.epsabs,
              cuba_parameters.verbose,
              cuba_parameters.seed,
              cuba_parameters.mineval,
              cuba_parameters.maxeval,
              cuba_parameters.nstart,
              cuba_parameters.nincrease,
              cuba_parameters.nbatch,
              cuba_parameters.gridno,
              cuba_parameters.statefile,
              cuba_parameters.spin,
              &neval,
              &fail,
              integral.data(),
              error.data(),
              prob.data()
              );

        std::cout << "VEGAS RESULT: " << " neval: " << neval << " fail: " << fail << std::endl;

        for( unsigned comp = 0; comp < cuba_parameters.ncomp; comp++)
            std::cout << integral.at(comp) << " +- " << error.at(comp) << " p = " << prob.at(comp) << std::endl;

        return secdecutil::GaussianUncertainty<cubareal>(integral.at(0),error.at(0)); // TODO - allow more components
    };
};

int main()
{
    // TODO - write method to parse arguments and check validity
    // User Specified Phase-space point
    const std::vector<%(name)s::real_t> real_parameters = { 0.9, 0.1 };
    const std::vector<%(name)s::complex_t> complex_parameters = {  };
    if ( real_parameters.size() != %(name)s::number_of_real_parameters )
        throw std::logic_error("Did not set the correct number of real parameters");
    if ( complex_parameters.size() != + %(name)s::number_of_complex_parameters )
        throw std::logic_error("Did not set the correct number of complex parameters");

    print_integral_info();

    const auto sector_integrands = %(name)s::make_integrands(real_parameters, complex_parameters);

    // const auto integrand_container_sum = sector_integrands.at(0) + sector_integrands.at(1); // Example how to add integrand containers

    // Add integrands of sectors (together flag)
    const auto all_sectors = std::accumulate(++sector_integrands.begin(), sector_integrands.end(), *sector_integrands.begin() );

    // Integrate
    auto result_all = secdecutil::deep_apply( all_sectors,  cuba_integrate<%(name)s::integrand_return_t,%(name)s::real_t,%(name)s::complex_t>() );

    std::cout << "-- All -- " << std::endl;
    std::cout << result_all << std::endl;

    return 0;
}

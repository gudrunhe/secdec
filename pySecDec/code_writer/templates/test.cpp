#include <iostream>
#include <stdexcept>
#include <vector>
#include <cstdint>
#include <numeric> // std::accumulate
#include <functional> // std::bind
#include <type_traits> // std::remove_const
#include <typeinfo>

#include <secdecutil/integrators/cuba.hpp> // Vegas
#include <secdecutil/series.hpp> // Series
#include <secdecutil/uncertainties.hpp> // UncorrelatedDeviation
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

    const auto sector_integrands = %(name)s::make_integrands(real_parameters, complex_parameters);

    // const auto integrand_container_sum = sector_integrands.at(0) + sector_integrands.at(1); // Example how to add integrand containers

    // Add integrands of sectors (together flag)
    const auto all_sectors = std::accumulate(++sector_integrands.begin(), sector_integrands.end(), *sector_integrands.begin() );

    // Integrate
    auto integrator = secdecutil::cuba::Vegas<%(name)s::integrand_return_t>();
    integrator.flags = 2; // verbose output --> see cuba manual
    auto result_all = secdecutil::deep_apply( all_sectors,  integrator.integrate );

    std::cout << "------------" << std::endl << std::endl;

    std::cout << "-- integral info -- " << std::endl;
    print_integral_info();
    std::cout << std::endl;

    std::cout << "-- integral without prefactor -- " << std::endl;
    std::cout << result_all << std::endl << std::endl;

    std::cout << "-- prefactor -- " << std::endl;
    auto prefactor = %(name)s::prefactor(real_parameters, complex_parameters);
    std::cout << prefactor << std::endl << std::endl;

    std::cout << "-- full result (prefactor*integral) -- " << std::endl;
    std::cout << prefactor*result_all << std::endl;

    return 0;
}

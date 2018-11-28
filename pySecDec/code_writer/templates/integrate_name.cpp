#include <iostream> // std::cout
#include <numeric> // std::accumulate
#include <vector> // std::vector

#include <secdecutil/integrators/cuba.hpp> // secdecutil::cuba::Vegas, secdecutil::cuba::Suave, secdecutil::cuba::Cuhre, secdecutil::cuba::Divonne
#include <secdecutil/integrators/qmc.hpp> // secdecutil::integrators::Qmc
#include <secdecutil/series.hpp> // secdecutil::Series
#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation
#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply

#include "%(name)s.hpp"

void print_integral_info()
{
    std::cout << "-- print_integral_info --" << std::endl;
    std::cout << "%(name)s::number_of_sectors " << %(name)s::number_of_sectors << std::endl;

    std::cout << "%(name)s::number_of_regulators " << %(name)s::number_of_regulators << std::endl;
    std::cout << "%(name)s::names_of_regulators ";
    for ( const auto& name : %(name)s::names_of_regulators )
        std::cout << " " << name;
    std::cout << std::endl;

    std::cout << "%(name)s::number_of_real_parameters " << %(name)s::number_of_real_parameters << std::endl;
    std::cout << "%(name)s::names_of_real_parameters ";
    for ( const auto& name : %(name)s::names_of_real_parameters )
        std::cout << " " << name;
    std::cout << std::endl;

    std::cout << "%(name)s::number_of_complex_parameters " << %(name)s::number_of_complex_parameters << std::endl;
    std::cout << "%(name)s::names_of_complex_parameters ";
    for ( const auto& name : %(name)s::names_of_complex_parameters )
        std::cout << " " << name;
    std::cout << std::endl;

    std::cout << "%(name)s::lowest_orders";
    for ( const auto& lowest_order : %(name)s::lowest_orders )
        std::cout << " " << lowest_order;
    std::cout << std::endl;

    std::cout << "%(name)s::highest_orders";
    for ( const auto& highest_order : %(name)s::highest_orders )
        std::cout << " " << highest_order;
    std::cout << std::endl;

    std::cout << "%(name)s::lowest_prefactor_orders";
    for ( const auto& highest_order : %(name)s::lowest_prefactor_orders )
        std::cout << " " << highest_order;
    std::cout << std::endl;

    std::cout << "%(name)s::highest_prefactor_orders";
    for ( const auto& highest_order : %(name)s::highest_prefactor_orders )
        std::cout << " " << highest_order;
    std::cout << std::endl;

    std::cout << "%(name)s::requested_orders";
    for ( const auto& requested_order : %(name)s::requested_orders )
        std::cout << " " << requested_order;
    std::cout << std::endl;
}

int main()
{
    // User Specified Phase-space point
    const std::vector<%(name)s::real_t> real_parameters = { /* EDIT: insert real parameter values here */ };
    const std::vector<%(name)s::complex_t> complex_parameters = { /* EDIT: insert complex parameter values here */ };

    // Generate the integrands (optimization of the contour if applicable)
    const std::vector<%(name)s::nested_series_t<%(name)s::integrand_t>> sector_integrands = %(name)s::make_integrands(real_parameters, complex_parameters);

    // Add integrands of sectors (together flag)
    const %(name)s::nested_series_t<%(name)s::integrand_t> all_sectors = std::accumulate(++sector_integrands.begin(), sector_integrands.end(), *sector_integrands.begin());

    // Integrate
    secdecutil::cuba::Vegas<%(name)s::integrand_return_t> integrator;
    integrator.flags = 2; // verbose output --> see cuba manual
    const %(name)s::nested_series_t<secdecutil::UncorrelatedDeviation<%(name)s::integrand_return_t>> result_all = secdecutil::deep_apply( all_sectors, integrator.integrate );

    std::cout << "------------" << std::endl << std::endl;

    std::cout << "-- integral info -- " << std::endl;
    print_integral_info();
    std::cout << std::endl;

    std::cout << "-- integral without prefactor -- " << std::endl;
    std::cout << result_all << std::endl << std::endl;

    std::cout << "-- prefactor -- " << std::endl;
    const %(name)s::nested_series_t<%(name)s::integrand_return_t> prefactor = %(name)s::prefactor(real_parameters, complex_parameters);
    std::cout << prefactor << std::endl << std::endl;

    std::cout << "-- full result (prefactor*integral) -- " << std::endl;
    std::cout << prefactor*result_all << std::endl;

    return 0;
}

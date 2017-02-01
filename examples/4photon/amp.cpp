#include <iostream>
#include <stdexcept>
#include <vector>
#include <cstdint>
#include <numeric> // std::accumulate
#include <functional> // std::bind
#include <type_traits> // std::remove_const
#include <typeinfo>
#include <cmath>

#include <secdecutil/integrators/cuba.hpp> // Vegas
#include <secdecutil/series.hpp> // Series
#include <secdecutil/uncertainties.hpp> // UncorrelatedDeviation
#include <secdecutil/sector_container.hpp> // SectorContainer to IntegrandContainer
#include <secdecutil/integrand_container.hpp> // IntegrandContainer
#include <secdecutil/deep_apply.hpp> // deep_apply

#include "yyyy_box6Dim/yyyy_box6Dim.hpp"
#include "yyyy_bubble/yyyy_bubble.hpp"

void print_integral_info()
{
    std::cout << "-- print_integral_info --" << std::endl;
    std::cout << "yyyy_box6Dim::number_of_sectors " << yyyy_box6Dim::number_of_sectors << std::endl;

    std::cout << "yyyy_box6Dim::number_of_regulators " << yyyy_box6Dim::number_of_regulators << std::endl;
    std::cout << "yyyy_box6Dim::names_of_regulators ";
    for ( const auto& name : yyyy_box6Dim::names_of_regulators )
        std::cout << " " << name;
    std::cout << std::endl;

    std::cout << "yyyy_box6Dim::number_of_real_parameters " << yyyy_box6Dim::number_of_real_parameters << std::endl;
    std::cout << "yyyy_box6Dim::names_of_real_parameters ";
    for ( const auto& name : yyyy_box6Dim::names_of_real_parameters )
        std::cout << " " << name;
    std::cout << std::endl;

    std::cout << "yyyy_box6Dim::number_of_complex_parameters " << yyyy_box6Dim::number_of_complex_parameters << std::endl;
    std::cout << "yyyy_box6Dim::names_of_complex_parameters ";
    for ( const auto& name : yyyy_box6Dim::names_of_complex_parameters )
        std::cout << " " << name;
    std::cout << std::endl;

    std::cout << "yyyy_box6Dim::lowest_orders";
    for ( const auto& lowest_order : yyyy_box6Dim::lowest_orders )
        std::cout << " " << lowest_order;
    std::cout << std::endl;

    std::cout << "yyyy_box6Dim::highest_orders";
    for ( const auto& highest_order : yyyy_box6Dim::highest_orders )
        std::cout << " " << highest_order;
    std::cout << std::endl;

    std::cout << "yyyy_box6Dim::lowest_prefactor_orders";
    for ( const auto& highest_order : yyyy_box6Dim::lowest_prefactor_orders )
        std::cout << " " << highest_order;
    std::cout << std::endl;

    std::cout << "yyyy_box6Dim::highest_prefactor_orders";
    for ( const auto& highest_order : yyyy_box6Dim::highest_prefactor_orders )
        std::cout << " " << highest_order;
    std::cout << std::endl;

    std::cout << "yyyy_box6Dim::requested_orders";
    for ( const auto& requested_order : yyyy_box6Dim::requested_orders )
        std::cout << " " << requested_order;
    std::cout << std::endl;
}

int main()
{
     double s = .9, t = -.1;
     double u = -s -t;
    // bubble
    // User Specified Phase-space point
    const std::vector<yyyy_bubble::real_t> bubble_parameters_t = { t };
    const std::vector<yyyy_bubble::real_t> bubble_parameters_u = { u };
    const std::vector<yyyy_bubble::complex_t> bubble_complex_parameters = {  };
    if ( (bubble_parameters_t.size() != yyyy_bubble::number_of_real_parameters) || (bubble_parameters_u.size() != yyyy_bubble::number_of_real_parameters) )
        throw std::logic_error("Did not set the correct number of real parameters");
    if ( bubble_complex_parameters.size() != + yyyy_bubble::number_of_complex_parameters )
        throw std::logic_error("Did not set the correct number of complex parameters");
    const auto bub_t_sector_integrands = yyyy_bubble::make_integrands(bubble_parameters_t, bubble_complex_parameters);
    const auto bub_u_sector_integrands = yyyy_bubble::make_integrands(bubble_parameters_u, bubble_complex_parameters);
// box
    const std::vector<yyyy_box6Dim::real_t> box_parameters = { t,u };
    const std::vector<yyyy_box6Dim::complex_t> box_complex_parameters = {  };
    if ( box_parameters.size() != yyyy_box6Dim::number_of_real_parameters )
        throw std::logic_error("Did not set the correct number of real parameters");
    if ( box_complex_parameters.size() != + yyyy_box6Dim::number_of_complex_parameters )
        throw std::logic_error("Did not set the correct number of complex parameters");

    const auto box_sector_integrands = yyyy_box6Dim::make_integrands(box_parameters, box_complex_parameters);

    // const auto integrand_container_sum = sector_integrands.at(0) + sector_integrands.at(1); // Example how to add integrand containers

    // Add integrands of sectors (together flag)
    const auto all_sectors_bub_t = std::accumulate(++bub_t_sector_integrands.begin(), bub_t_sector_integrands.end(), *bub_t_sector_integrands.begin() );
    const auto all_sectors_bub_u = std::accumulate(++bub_u_sector_integrands.begin(), bub_u_sector_integrands.end(), *bub_u_sector_integrands.begin() );
    const auto all_sectors_box = std::accumulate(++box_sector_integrands.begin(), box_sector_integrands.end(), *box_sector_integrands.begin() );

    // define the integrator
    auto complex_integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    complex_integrator.flags = 2; // verbose output
    complex_integrator.epsrel = 1e-4;
    complex_integrator.epsabs = 1e-4;

    // Integrate bubble_t
    auto result_bub_t = secdecutil::deep_apply( all_sectors_bub_t,  complex_integrator.integrate )
                      * yyyy_bubble::prefactor(bubble_parameters_t, bubble_complex_parameters)
                      ;
    // Integrate bubble_u
    auto result_bub_u = secdecutil::deep_apply( all_sectors_bub_u,  complex_integrator.integrate )
                      * yyyy_bubble::prefactor(bubble_parameters_u, bubble_complex_parameters)
                      ;
    // Integrate box

    auto result_box = secdecutil::deep_apply( all_sectors_box,  complex_integrator.integrate )
                    * yyyy_box6Dim::prefactor(box_parameters, box_complex_parameters)
                    ;

    auto amp = -8.*( 1. + (t*t + u*u)/s * result_box +  (t-u)/s*( result_bub_u-result_bub_t ) );

    std::function<secdecutil::UncorrelatedDeviation<std::complex<double>>(secdecutil::UncorrelatedDeviation<std::complex<double>>)> conjugate =
        [] (secdecutil::UncorrelatedDeviation<std::complex<double>> x) {
            return secdecutil::UncorrelatedDeviation<std::complex<double>>{std::conj(x.value), x.uncertainty};
        };

//
    std::cout << "------------" << std::endl << std::endl;

//    std::cout << "-- integral info -- " << std::endl;
//    print_integral_info();
//    std::cout << std::endl;
    constexpr auto M_pi = 3.14159265358979323846;

    std::cout << " amplitude M++-- " << std::endl;
    std::cout << amp  /* *secdecutil::deep_apply(amp, conjugate) */ << std::endl << std::endl;

    std::cout << " analytic result " << std::endl;
    auto analytic = -8.*( 1. + (t-u)/s*std::log(t/u)  + (t*t + u*u)/(2.*s*s)*( M_pi*M_pi + std::log(t/u)*std::log(t/u) ) );
    std::cout << analytic << std::endl << std::endl;


    return 0;
}

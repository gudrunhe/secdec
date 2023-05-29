#include "catch_amalgamated.hpp"
using Catch::Approx;

#include <complex> // std::complex
#include <numeric> // std::accumulate
#include <secdecutil/integrators/cuba.hpp> // Vegas
#include <secdecutil/deep_apply.hpp> // deep_apply
#include <secdecutil/series.hpp> // Series

#define QUOTE(ARG) #ARG
#define QUOTE_EXPAND(ARG) QUOTE(ARG)
#define INTEGRAL_NAME difference

#include QUOTE_EXPAND(INTEGRAL_NAME.hpp)

TEST_CASE( "check result", "[difference]" ) {
    
    // User Specified Phase-space point
    const std::vector<INTEGRAL_NAME::real_t> real_parameters = { };
    const std::vector<INTEGRAL_NAME::complex_t> complex_parameters = { };

    // define and comnfigure integrator
    auto integrator = secdecutil::cuba::Vegas<difference::integrand_return_t>();
    integrator.flags = 2; // verbose output --> see cuba manual
    const double epsrel = 1e-2; integrator.epsrel = epsrel;

    // Construct the amplitudes
    std::vector<INTEGRAL_NAME::nested_series_t<INTEGRAL_NAME::sum_t>> unwrapped_amplitudes =
        INTEGRAL_NAME::make_amplitudes(real_parameters, complex_parameters, "../difference/difference_data", integrator);

    // Pack amplitudes into handler
    INTEGRAL_NAME::handler_t<INTEGRAL_NAME::amplitudes_t> amplitudes
    (
        unwrapped_amplitudes, integrator.epsrel, integrator.epsabs
        // further optional arguments: epsrel, epsabs, maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
    );
    
    // integrate
    const std::vector<INTEGRAL_NAME::nested_series_t<secdecutil::UncorrelatedDeviation<INTEGRAL_NAME::integrand_return_t>>> result = amplitudes.evaluate();
    auto result_with_prefactor = result.at(0);

    // analytic result obtained with Mathematica
    constexpr std::complex<double> I(0,1);
    secdecutil::Series<std::complex<double>> target_result
    (

        0, // lowest order in epsilon
        2, // highest computed order in epsilon

        {
             1.64493406684822643647241516664602518923 + 3.1415926535897932384626433832795028842 * I, // eps ** 0
             2.08781123053685858754509217178101012328 - 6.2831853071795864769252867665590057684 * I, // eps ** 1
            -5.94029019737039970544633397517750766917 + 4.2570651807194096861418776386549427857 * I, // eps ** 2
//           5.77945251635494087034720012662916969501 - 2.2309450542592328953584685107508798034 * I, // eps ** 3 // this is slow
        },

        true, // series is truncated above; i.e. "+ O(eps**4)"
        "eps" // the expansion parameter

    );

    // basic checks
    REQUIRE(       result_with_prefactor.get_order_min() == target_result.get_order_min()       );
    REQUIRE(       result_with_prefactor.get_order_max() == target_result.get_order_max()       );
    REQUIRE( result_with_prefactor.get_truncated_above() == target_result.get_truncated_above() );
    REQUIRE(   result_with_prefactor.expansion_parameter == target_result.expansion_parameter   );

    std::cout << "----------------" << std::endl << std::endl;

    for (int order = target_result.get_order_min() ; order <= target_result.get_order_max() ; ++ order)
    {
        std::cout << "checking order \"eps^" << order << "\" ..." << std::endl;

        // check that the uncertainties are reasonable
        REQUIRE( result_with_prefactor.at(order).uncertainty.real() <= abs(2*epsrel * target_result.at(order).real()) );
        REQUIRE( result_with_prefactor.at(order).uncertainty.imag() <= abs(2*epsrel * target_result.at(order).imag()) );

        // check that the desired uncertainties are reached
        REQUIRE( result_with_prefactor.at(order).uncertainty.real() <= abs(epsrel * result_with_prefactor.at(order).value) );
        REQUIRE( result_with_prefactor.at(order).uncertainty.imag() <= abs(epsrel * result_with_prefactor.at(order).value) );

        // check integral value
        REQUIRE(  result_with_prefactor.at(order).value.real() == Approx( target_result.at(order).real() ).epsilon( epsrel )  );
        REQUIRE(  result_with_prefactor.at(order).value.imag() == Approx( target_result.at(order).imag() ).epsilon( epsrel )  );

        std::cout << "----------------" << std::endl << std::endl;
    }

};

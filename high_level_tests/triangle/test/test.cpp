#include "catch_amalgamated.hpp"
using Catch::Approx;

#include "triangle.hpp"

#include <complex> // std::complex
#include <numeric> // std::accumulate
#include <secdecutil/integrators/cuba.hpp> // Vegas
#include <secdecutil/deep_apply.hpp> // deep_apply
#include <secdecutil/series.hpp> // Series

#define QUOTE(ARG) #ARG
#define QUOTE_EXPAND(ARG) QUOTE(ARG)
#define INTEGRAL_NAME triangle

#include QUOTE_EXPAND(INTEGRAL_NAME.hpp)

constexpr std::complex<double> I(0,1);

TEST_CASE( "check number of parameters", "[triangle]" ) {
    
    // define and configure dummy integrator
    auto integrator = secdecutil::cuba::Vegas<triangle::integrand_return_t>();
    integrator.flags = 0; // verbose output --> see cuba manual
    integrator.maxeval = 2e6;
    
    REQUIRE_THROWS_AS( INTEGRAL_NAME::make_amplitudes(/* real_parameters */ {}, /* complex_parameters */ {}, "../triangle/triangle_data", integrator), std::logic_error );

    try {

      triangle::make_amplitudes(/* real_parameters */ {}, /* complex_parameters */ {}, "../triangle/triangle_data", integrator);

    } catch (std::logic_error error) {

      REQUIRE( error.what() == std::string("Called \"triangle_integral::make_integrands\" with 0 \"real_parameters\" (2 expected).") );

    }

};

const secdecutil::Series<double> target_prefactor
                                 (
                                        0, // lowest order in epsilon
                                        2, // highest computed order in epsilon
                                        {
                                              1.0,                    // eps ** 0
                                              0.84556867019693427879, // eps ** 1
                                              1.6473613217057587792,  // eps ** 2
                                         },
                                         true, // series is truncated above; i.e. "+ O(eps)"
                                         "eps" // the expansion parameter
                                 );

//TEST_CASE( "check prefactor", "[triangle]" ) {
//
//    SECTION( "s = 9, msq = 1" ) {
//
//        double s = 9, msq = 1;
//        const auto prefactor = triangle::prefactor(/* real_parameters */ {s,msq}, /* complex_parameters */ {});
//        REQUIRE( prefactor == target_prefactor );
//
//    }
//
//    SECTION( "s = 0.9, msq = 0.1" ) {
//
//        double s = 0.9, msq = 0.1;
//        const auto prefactor = triangle::prefactor(/* real_parameters */ {s,msq}, /* complex_parameters */ {});
//        REQUIRE( prefactor == target_prefactor );
//
//    }
//
//}

void check_pySecDec_triangle(double s, double msq, secdecutil::Series<std::complex<double>> target_result_with_prefactor, const secdecutil::Series<double> epsrels)
    {

        // define and configure integrator
        auto integrator = secdecutil::cuba::Vegas<triangle::integrand_return_t>();
        integrator.flags = 0; // verbose output --> see cuba manual
        integrator.maxeval = 2e6;
        
        // Construct the amplitudes
        std::vector<INTEGRAL_NAME::nested_series_t<INTEGRAL_NAME::sum_t>> unwrapped_amplitudes =
            INTEGRAL_NAME::make_amplitudes({s,msq}, {}, "../triangle/triangle_data", integrator);

        // Pack amplitudes into handler
        INTEGRAL_NAME::handler_t<INTEGRAL_NAME::amplitudes_t> amplitudes
        (
            unwrapped_amplitudes, integrator.epsrel, integrator.epsabs
            // further optional arguments: epsrel, epsabs, maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
        );
        for (int order = target_result_with_prefactor.get_order_min() ; order <= target_result_with_prefactor.get_order_max() ; ++ order)
        {
            amplitudes.expression.at(0).at(order).epsrel = epsrels.at(order);
        }

        // integrate
        const std::vector<INTEGRAL_NAME::nested_series_t<secdecutil::UncorrelatedDeviation<INTEGRAL_NAME::integrand_return_t>>> result = amplitudes.evaluate();
        auto result_with_prefactor = result.at(0);

        std::cout << "----------------" << std::endl << std::endl;

        for (int order = target_result_with_prefactor.get_order_min() ; order <= target_result_with_prefactor.get_order_max() ; ++ order)
        {
            std::cout << "checking order \"eps^" << order << "\" ..." << std::endl;

            // check that the uncertainties are reasonable
            REQUIRE( result_with_prefactor.at(order).uncertainty.real() <= abs(2*epsrels.at(order) * target_result_with_prefactor.at(order).real()) );
            REQUIRE( result_with_prefactor.at(order).uncertainty.imag() <= abs(2*epsrels.at(order) * target_result_with_prefactor.at(order).imag()) );

            // check that the desired uncertainties are reached
            REQUIRE( result_with_prefactor.at(order).uncertainty.real() <= abs(epsrels.at(order) * result_with_prefactor.at(order).value) );
            REQUIRE( result_with_prefactor.at(order).uncertainty.imag() <= abs(epsrels.at(order) * result_with_prefactor.at(order).value) );

            // check integral value
            REQUIRE(  result_with_prefactor.at(order).value.real() == Approx( target_result_with_prefactor.at(order).real() ).epsilon( epsrels.at(order) )  );
            REQUIRE(  result_with_prefactor.at(order).value.imag() == Approx( target_result_with_prefactor.at(order).imag() ).epsilon( epsrels.at(order) )  );

            std::cout << "----------------" << std::endl << std::endl;
        }

}

TEST_CASE( "check result with prefactor", "[triangle]" ) {

    SECTION( "s = 9, msq = 1" ) {

        // set the relative target uncertainties
        const secdecutil::Series<double> epsrels(-2,0,
                                                    {
                                                        5e-3, // order eps ** -2
                                                        1e-2, // order eps ** -1
                                                        1e-1, // order eps **  0
                                                    }
                                                 );

        // result obtained with secdec-3
        secdecutil::Series<std::complex<double>> target_result
        (

            -2, // lowest order in epsilon
             0, // highest computed order in epsilon

            {
                -0.03805 - 0.07465 * I, // eps ** -2
                 0.3117  + 0.2378  * I, // eps ** -1
                -1.24    + 0.165   * I, // eps **  0
            },

            true, // series is truncated above; i.e. "+ O(eps)"
            "eps" // the expansion parameter

        );

        check_pySecDec_triangle(/* s */ 9.0, /* msq */ 1.0, target_prefactor*target_result, epsrels);

    }

    SECTION( "s = 0.9, msq = 0.1" ) { // numerically difficult

        // set the relative target uncertainties
        const secdecutil::Series<double> epsrels(-2,0,
                                                    {
                                                        1e-2, // order eps ** -2
                                                        5e-2, // order eps ** -1
                                                        5e-1, // order eps **  0
                                                    }
                                                 );

        // result obtained with secdec-3
        secdecutil::Series<std::complex<double>> target_result
        (

            -2, // lowest order in epsilon
             0, // highest computed order in epsilon

            {
                - 3.80 -  7.465 * I, // eps ** -2
                 13.6  - 10.6   * I, // eps ** -1
                -20.   + 46.    * I, // eps **  0
            },

            true, // series is truncated above; i.e. "+ O(eps)"
            "eps" // the expansion parameter

        );

        check_pySecDec_triangle(/* s */ 0.9, /* msq */ 0.1, target_prefactor*target_result, epsrels);

    }

};

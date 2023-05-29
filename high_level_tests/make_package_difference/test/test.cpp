#include "catch_amalgamated.hpp"
using Catch::Approx;

#include "difference.hpp"

#include <complex> // std::complex
#include <numeric> // std::accumulate
#include <secdecutil/integrators/cuba.hpp> // Vegas
#include <secdecutil/deep_apply.hpp> // deep_apply
#include <secdecutil/series.hpp> // Series

TEST_CASE( "check result", "[difference]" ) {

    // get the integrands
    const auto sector_integrands = difference::make_integrands(/* real_parameters */ {}, /* complex_parameters */ {});

    // add integrands of sectors (together flag)
    const auto all_sectors = std::accumulate(++sector_integrands.begin(), sector_integrands.end(), *sector_integrands.begin() );

    // define and comnfigure integrator
    auto integrator = secdecutil::cuba::Vegas<difference::integrand_return_t>();
    integrator.flags = 2; // verbose output --> see cuba manual
    const double epsrel = 1e-2; integrator.epsrel = epsrel;

    // integrate
    auto result = secdecutil::deep_apply( all_sectors,  integrator.integrate );
    //          * difference::prefactor(/* real_parameters */ {}, /* complex_parameters */ {});
    // the prefactor is just one --> don't need to multiply it


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
    REQUIRE(       result.get_order_min() == target_result.get_order_min()       );
    REQUIRE(       result.get_order_max() == target_result.get_order_max()       );
    REQUIRE( result.get_truncated_above() == target_result.get_truncated_above() );
    REQUIRE(   result.expansion_parameter == target_result.expansion_parameter   );

    std::cout << "----------------" << std::endl << std::endl;

    for (int order = target_result.get_order_min() ; order <= target_result.get_order_max() ; ++ order)
    {
        std::cout << "checking order \"eps^" << order << "\" ..." << std::endl;

        // check that the uncertainties are reasonable
        REQUIRE( result.at(order).uncertainty.real() <= std::abs(2*epsrel * target_result.at(order).real()) );
        REQUIRE( result.at(order).uncertainty.imag() <= std::abs(2*epsrel * target_result.at(order).imag()) );

        // check that the desired uncertainties are reached
        REQUIRE( result.at(order).uncertainty.real() <= std::abs(epsrel * result.at(order).value.real()) );
        REQUIRE( result.at(order).uncertainty.imag() <= std::abs(epsrel * result.at(order).value.imag()) );

        // check integral value
        REQUIRE(  result.at(order).value.real() == Approx( target_result.at(order).real() ).epsilon( epsrel )  );
        REQUIRE(  result.at(order).value.imag() == Approx( target_result.at(order).imag() ).epsilon( epsrel )  );

        std::cout << "----------------" << std::endl << std::endl;
    }

};

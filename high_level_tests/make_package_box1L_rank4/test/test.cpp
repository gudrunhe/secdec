#include "catch_amalgamated.hpp"
using Catch::Approx;

#include <complex> // std::complex
#include <numeric> // std::accumulate
#include <secdecutil/integrators/cuba.hpp> // Vegas
#include <secdecutil/deep_apply.hpp> // deep_apply
#include <secdecutil/series.hpp> // Series

#define QUOTE(ARG) #ARG
#define QUOTE_EXPAND(ARG) QUOTE(ARG)
#define INTEGRAL_NAME box1L_rank4

#include QUOTE_EXPAND(INTEGRAL_NAME.hpp)

TEST_CASE( "check result", "[INTEGRAL_NAME]" ) {

    // User Specified Phase-space point
    const std::vector<INTEGRAL_NAME::real_t> real_parameters = { 16.0, -75.0, 1.0 };
    const std::vector<INTEGRAL_NAME::complex_t> complex_parameters = {  };

    // get the integrands
    const auto sector_integrands = INTEGRAL_NAME::make_integrands(real_parameters, complex_parameters);

    // add integrands of sectors (together flag)
    const auto all_sectors = std::accumulate(++sector_integrands.begin(), sector_integrands.end(), *sector_integrands.begin() );

    // define and comnfigure integrator
    auto integrator = secdecutil::cuba::Cuhre<INTEGRAL_NAME::integrand_return_t>();
    // integrator.flags = 2; // verbose output --> see cuba manual
    integrator.maxeval = 1e7;
    integrator.together = true;
    const double epsrel = 1e-10; integrator.epsrel = epsrel;
    const double epsabs = 1e-10; integrator.epsabs = epsabs;

    // integrate
    auto result_without_prefactor = secdecutil::deep_apply( all_sectors,  integrator.integrate );
    auto prefactor = INTEGRAL_NAME::prefactor(real_parameters, complex_parameters);
    auto result_with_prefactor = result_without_prefactor * prefactor;


    // target result, obtained in a long run of pySecDec
    constexpr std::complex<double> I(0,1);
    secdecutil::Series<std::complex<double>> target_result_without_prefactor
    (

        0, // lowest order in epsilon
        1, // highest computed order in epsilon

        {
              97.52083333333 +  0.0          * I, // eps ** 0
            -124.937368867   - 11.9030583278 * I, // eps ** 1
        },

        true, // series is truncated above; i.e. "+ O(eps**2)"
        "eps" // the expansion parameter
    );
    secdecutil::Series<std::complex<double>> target_prefactor
    (

        -1, // lowest order in epsilon
         0, // highest computed order in epsilon

        {
             1.0,              // eps ** -1
            -0.57721566490153, // eps **  0
        },

        true, // series is truncated above; i.e. "+ O(eps**2)"
        "eps" // the expansion parameter
    );
    secdecutil::Series<std::complex<double>> target_result_with_prefactor
    (

        -1, // lowest order in epsilon
         0, // highest computed order in epsilon

        {
              97.52083333333 +  0.0            * I, // eps ** -1
            -181.22792152123 - 11.903058327787 * I, // eps **  0
        },

        true, // series is truncated above; i.e. "+ O(eps**2)"
        "eps" // the expansion parameter
    );

    // basic checks
    REQUIRE(       result_without_prefactor.get_order_min() == target_result_without_prefactor.get_order_min()       );
    REQUIRE(       result_without_prefactor.get_order_max() == target_result_without_prefactor.get_order_max()       );
    REQUIRE( result_without_prefactor.get_truncated_above() == target_result_without_prefactor.get_truncated_above() );
    REQUIRE(   result_without_prefactor.expansion_parameter == target_result_without_prefactor.expansion_parameter   );
    REQUIRE(          result_with_prefactor.get_order_min() == target_result_with_prefactor.get_order_min()          );
    REQUIRE(          result_with_prefactor.get_order_max() == target_result_with_prefactor.get_order_max()          );
    REQUIRE(    result_with_prefactor.get_truncated_above() == target_result_with_prefactor.get_truncated_above()    );
    REQUIRE(      result_with_prefactor.expansion_parameter == target_result_with_prefactor.expansion_parameter      );

    std::cout << "----------------" << std::endl << std::endl;

    for (int order = target_result_with_prefactor.get_order_min() ; order <= target_result_with_prefactor.get_order_max() ; ++ order)
    {
        std::cout << "checking order \"eps^" << order << "\" ..." << std::endl;

        // check that the uncertainties are reasonable
        REQUIRE( result_without_prefactor.at(order+1).uncertainty.real() <= std::abs(2*epsrel * target_result_without_prefactor.at(order+1).real()) );
        if (  target_result_without_prefactor.at(order+1).imag() != 0.0  )
            REQUIRE( result_without_prefactor.at(order+1).uncertainty.imag() <= std::abs(2*epsrel * target_result_without_prefactor.at(order+1).imag()) );
        REQUIRE(      result_with_prefactor.at(order).uncertainty.real() <= std::abs(2*epsrel * target_result_with_prefactor.at(order).real())      );
        if (  target_result_with_prefactor.at(order).imag() != 0.0  )
            REQUIRE(      result_with_prefactor.at(order).uncertainty.imag() <= std::abs(2*epsrel * target_result_with_prefactor.at(order).imag())      );

        // check that the desired uncertainty is reached before multiplying the prefactor
        REQUIRE( result_without_prefactor.at(order+1).uncertainty.real() <= std::abs(epsrel * result_without_prefactor.at(order+1).value.real()) );
        if (  target_result_without_prefactor.at(order+1).imag() == 0.0  )
            REQUIRE(      result_without_prefactor.at(order+1).uncertainty.imag() <= epsabs      );
        else
            REQUIRE( result_without_prefactor.at(order+1).uncertainty.imag() <= std::abs(epsrel * result_without_prefactor.at(order+1).value.imag()) );

        // check values
        REQUIRE(  result_without_prefactor.at(order+1).value.real() == Approx( target_result_without_prefactor.at(order+1).real() ).epsilon( epsrel )  );
        if (  target_result_without_prefactor.at(order+1).imag() == 0.0  )
            REQUIRE(  std::abs( result_without_prefactor.at(order+1).value.imag() ) <= epsabs  );
        else
            REQUIRE(  result_without_prefactor.at(order+1).value.imag() == Approx( target_result_without_prefactor.at(order+1).imag() ).epsilon( epsrel )  );
        REQUIRE(  prefactor.at(order).real() == Approx( target_prefactor.at(order).real() ).epsilon( 1e-13 )  );
        REQUIRE(  prefactor.at(order).imag() == Approx( target_prefactor.at(order).imag() ).epsilon( 1e-13 )  );
        REQUIRE(  result_with_prefactor.at(order).value.real() == Approx( target_result_with_prefactor.at(order).real() ).epsilon( 10.0*epsrel )  );
        if (  target_result_with_prefactor.at(order).imag() == 0.0  )
            REQUIRE(  result_with_prefactor.at(order).value.imag() <= epsabs  );
        else
            REQUIRE(  result_with_prefactor.at(order).value.imag() == Approx( target_result_with_prefactor.at(order).imag() ).epsilon( 10.0*epsrel )  );

        std::cout << "----------------" << std::endl << std::endl;
    }

};

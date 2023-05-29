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

    // define and configure integrator
    auto integrator = secdecutil::cuba::Cuhre<INTEGRAL_NAME::integrand_return_t>();
    // integrator.flags = 2; // verbose output --> see cuba manual
    integrator.maxeval = 1e7;
    integrator.together = true;
    const double epsrel = 1e-10; integrator.epsrel = epsrel;
    const double epsabs = 1e-10; integrator.epsabs = epsabs;
    
    // Construct the amplitudes
    std::vector<INTEGRAL_NAME::nested_series_t<INTEGRAL_NAME::sum_t>> unwrapped_amplitudes =
        INTEGRAL_NAME::make_amplitudes(real_parameters, complex_parameters, "../box1L_rank4/box1L_rank4_data", integrator);

    // Pack amplitudes into handler
    INTEGRAL_NAME::handler_t<INTEGRAL_NAME::amplitudes_t> amplitudes
    (
        unwrapped_amplitudes, integrator.epsrel, integrator.epsabs
        // further optional arguments: epsrel, epsabs, maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
    );
    
    // integrate
    const std::vector<INTEGRAL_NAME::nested_series_t<secdecutil::UncorrelatedDeviation<INTEGRAL_NAME::integrand_return_t>>> result = amplitudes.evaluate();
    auto result_with_prefactor = result.at(0);


    // target result, obtained in a long run of pySecDec
    constexpr std::complex<double> I(0,1);
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
    REQUIRE(          result_with_prefactor.get_order_min() == target_result_with_prefactor.get_order_min()          );
    REQUIRE(          result_with_prefactor.get_order_max() == target_result_with_prefactor.get_order_max()          );
    REQUIRE(    result_with_prefactor.get_truncated_above() == target_result_with_prefactor.get_truncated_above()    );
    REQUIRE(      result_with_prefactor.expansion_parameter == target_result_with_prefactor.expansion_parameter      );

    std::cout << "----------------" << std::endl << std::endl;

    for (int order = target_result_with_prefactor.get_order_min() ; order <= target_result_with_prefactor.get_order_max() ; ++ order)
    {
        std::cout << "checking order \"eps^" << order << "\" ..." << std::endl;

        // check that the uncertainties are reasonable
        REQUIRE(      result_with_prefactor.at(order).uncertainty.real() <= std::abs(3*epsrel * target_result_with_prefactor.at(order).real())      );
        if (  target_result_with_prefactor.at(order).imag() != 0.0  )
            REQUIRE(      result_with_prefactor.at(order).uncertainty.imag() <= std::abs(3*epsrel * target_result_with_prefactor.at(order).imag())      );

        // check values
        REQUIRE(  result_with_prefactor.at(order).value.real() == Approx( target_result_with_prefactor.at(order).real() ).epsilon( 10.0*epsrel )  );
        if (  target_result_with_prefactor.at(order).imag() == 0.0  )
            REQUIRE(  result_with_prefactor.at(order).value.imag() <= epsabs  );
        else
            REQUIRE(  result_with_prefactor.at(order).value.imag() == Approx( target_result_with_prefactor.at(order).imag() ).epsilon( 10.0*epsrel )  );

        std::cout << "----------------" << std::endl << std::endl;
    }

};

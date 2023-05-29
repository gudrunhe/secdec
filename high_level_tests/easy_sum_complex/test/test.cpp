#include "catch_amalgamated.hpp"
using Catch::Approx;

#include <complex> // std::complex
#include <numeric> // std::accumulate
#include <secdecutil/integrators/integrator.hpp> // MultiIntegrator
#include <secdecutil/integrators/cuba.hpp> // Vegas
#include <secdecutil/integrators/cquad.hpp> // CQuad
#include <secdecutil/integrators/qmc.hpp> // Qmc
#include <secdecutil/deep_apply.hpp> // deep_apply
#include <secdecutil/series.hpp> // Series

#define QUOTE(ARG) #ARG
#define QUOTE_EXPAND(ARG) QUOTE(ARG)
#define INTEGRAL_NAME easy_sum_complex

#include QUOTE_EXPAND(INTEGRAL_NAME.hpp)

TEST_CASE( "check result qmc", "[INTEGRAL_NAME]" ) {

    SECTION("default integral transform") {

        // User Specified Phase-space point
        const std::vector<INTEGRAL_NAME::real_t> real_parameters = { 0.1 };
        const std::vector<INTEGRAL_NAME::complex_t> complex_parameters = { };

        // define and configure integrator
        secdecutil::integrators::Qmc<INTEGRAL_NAME::integrand_return_t,
        INTEGRAL_NAME::maximal_number_of_integration_variables,
        ::integrators::transforms::Korobov<3>::type
        #ifdef SECDEC_WITH_CUDA
            ,INTEGRAL_NAME::cuda_together_integrand_t
        #endif
        > integrator;
        integrator.maxeval = 1e6;
        const double epsrel = 1e-4; integrator.epsrel = epsrel;
        const double epsabs = 1e-7; integrator.epsabs = epsabs;

        // Construct the amplitudes
        std::vector<INTEGRAL_NAME::nested_series_t<INTEGRAL_NAME::sum_t>> unwrapped_amplitudes =
            INTEGRAL_NAME::make_amplitudes(real_parameters, complex_parameters, "../easy_sum_complex/easy_sum_complex_data", integrator);

        // Pack amplitudes into handler
        INTEGRAL_NAME::handler_t<INTEGRAL_NAME::amplitudes_t> amplitudes
        (
            unwrapped_amplitudes, integrator.epsrel, integrator.epsabs
            // further optional arguments: epsrel, epsabs, maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
        );
        
        // integrate
        const std::vector<INTEGRAL_NAME::nested_series_t<secdecutil::UncorrelatedDeviation<INTEGRAL_NAME::integrand_return_t>>> result = amplitudes.evaluate();
        auto result_with_prefactor_1 = result.at(0);
        auto result_with_prefactor_2 = result.at(1);


        // target result, obtained in a long run of pySecDec
        constexpr std::complex<double> I(0,1);
        secdecutil::Series<std::complex<double>> target_result_with_prefactor_1
        (

            -1, // lowest order in epsilon
             0, // highest computed order in epsilon

            {
                 0.2                     + 0.0 * I, // eps ** -1
                 0.22962348064032504712 + 0.0 * I, // eps **  0
            },

            true, // series is truncated above; i.e. "+ O(eps**2)"
            "eps" // the expansion parameter
        );
        secdecutil::Series<std::complex<double>> target_result_with_prefactor_2
        (

            -2, // lowest order in epsilon
             0, // highest computed order in epsilon

            {
                 0.05                     + 0.0 * I, // eps ** -2
                 0.015342640972002734529  + 0.0 * I, // eps ** -1
                 0.0033313156240476989125 + 0.0 * I  // eps **  0
            },

            true, // series is truncated above; i.e. "+ O(eps**2)"
            "eps" // the expansion parameter
        );

        // basic checks
        REQUIRE(          result_with_prefactor_1.get_order_min() == target_result_with_prefactor_1.get_order_min()          );
        REQUIRE(          result_with_prefactor_1.get_order_max() == target_result_with_prefactor_1.get_order_max()          );
        REQUIRE(    result_with_prefactor_1.get_truncated_above() == target_result_with_prefactor_1.get_truncated_above()    );
        REQUIRE(      result_with_prefactor_1.expansion_parameter == target_result_with_prefactor_1.expansion_parameter      );
        
        REQUIRE(          result_with_prefactor_2.get_order_min() == target_result_with_prefactor_2.get_order_min()          );
        REQUIRE(          result_with_prefactor_2.get_order_max() == target_result_with_prefactor_2.get_order_max()          );
        REQUIRE(    result_with_prefactor_2.get_truncated_above() == target_result_with_prefactor_2.get_truncated_above()    );
        REQUIRE(      result_with_prefactor_2.expansion_parameter == target_result_with_prefactor_2.expansion_parameter      );

        std::cout << "----------------" << std::endl << std::endl;

        for (int order = target_result_with_prefactor_1.get_order_min() ; order <= target_result_with_prefactor_1.get_order_max() ; ++ order)
        {
            std::cout << "checking order \"eps^" << order << "\" ..." << std::endl;

            // check that the uncertainties are reasonable
            if (  target_result_with_prefactor_1.at(order).real() != 0.0  )
                REQUIRE(      result_with_prefactor_1.at(order).uncertainty.real() <= std::abs(2*epsrel * target_result_with_prefactor_1.at(order).real())      );
            if (  target_result_with_prefactor_1.at(order).imag() != 0.0  )
                REQUIRE(      result_with_prefactor_1.at(order).uncertainty.imag() <= std::abs(2*epsrel * target_result_with_prefactor_1.at(order).imag())      );

            // check values
            REQUIRE(  result_with_prefactor_1.at(order).value.real() == Approx( target_result_with_prefactor_1.at(order).real() ).epsilon( 10.0*epsrel )  );
            if (  target_result_with_prefactor_1.at(order).imag() == 0.0  )
                REQUIRE(  result_with_prefactor_1.at(order).value.imag() <= epsabs  );
            else
                REQUIRE(  result_with_prefactor_1.at(order).value.imag() == Approx( target_result_with_prefactor_1.at(order).imag() ).epsilon( 10.0*epsrel )  );

            std::cout << "----------------" << std::endl << std::endl;
        }
        
        for (int order = target_result_with_prefactor_2.get_order_min() ; order <= target_result_with_prefactor_2.get_order_max() ; ++ order)
        {
            std::cout << "checking order \"eps^" << order << "\" ..." << std::endl;

            // check that the uncertainties are reasonable
            if (  target_result_with_prefactor_2.at(order).real() != 0.0  )
                REQUIRE(      result_with_prefactor_2.at(order).uncertainty.real() <= std::abs(2*epsrel * target_result_with_prefactor_2.at(order).real())      );
            if (  target_result_with_prefactor_2.at(order).imag() != 0.0  )
                REQUIRE(      result_with_prefactor_2.at(order).uncertainty.imag() <= std::abs(2*epsrel * target_result_with_prefactor_2.at(order).imag())      );

            // check values
            REQUIRE(  result_with_prefactor_2.at(order).value.real() == Approx( target_result_with_prefactor_2.at(order).real() ).epsilon( 10.0*epsrel )  );
            if (  target_result_with_prefactor_2.at(order).imag() == 0.0  )
                REQUIRE(  result_with_prefactor_2.at(order).value.imag() <= epsabs  );
            else
                REQUIRE(  result_with_prefactor_2.at(order).value.imag() == Approx( target_result_with_prefactor_2.at(order).imag() ).epsilon( 10.0*epsrel )  );

            std::cout << "----------------" << std::endl << std::endl;
        }
    }

};

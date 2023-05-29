#include "catch_amalgamated.hpp"
using Catch::Approx;

#include <complex> // std::complex
#include <numeric> // std::accumulate
#include <secdecutil/integrators/qmc.hpp> // Qmc
#include <secdecutil/deep_apply.hpp> // deep_apply
#include <secdecutil/series.hpp> // Series

#define QUOTE(ARG) #ARG
#define QUOTE_EXPAND(ARG) QUOTE(ARG)
#define INTEGRAL_NAME one_integration_variable

#include QUOTE_EXPAND(INTEGRAL_NAME.hpp)

TEST_CASE( "check result with qmc", "[INTEGRAL_NAME]" ) {

    SECTION("default integral transform") {

        // User Specified Phase-space point
        const std::vector<INTEGRAL_NAME::real_t> real_parameters = {};
        const std::vector<INTEGRAL_NAME::complex_t> complex_parameters = {};

        // define and configure integrator
        secdecutil::integrators::Qmc<INTEGRAL_NAME::integrand_return_t,
        INTEGRAL_NAME::maximal_number_of_integration_variables,
        ::integrators::transforms::Korobov<3>::type
        #ifdef SECDEC_WITH_CUDA
            ,INTEGRAL_NAME::cuda_together_integrand_t
        #endif
        > integrator;
        integrator.minn = 1e5;
        integrator.maxeval = 1e6;

        const double epsrel = 1e-11;
        const double epsabs = 1e-10;
        integrator.epsrel = 1e-11;
        integrator.epsabs = 1e-10;
        integrator.randomgenerator.seed(798431);

        // Construct the amplitudes
        std::vector<INTEGRAL_NAME::nested_series_t<INTEGRAL_NAME::sum_t>> unwrapped_amplitudes =
            INTEGRAL_NAME::make_amplitudes(real_parameters, complex_parameters, "../one_integration_variable/one_integration_variable_data", integrator);

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
        secdecutil::Series<double> target_result_with_prefactor
        (

            -1, // lowest order in epsilon
             1, // highest computed order in epsilon

            {
                 0.0 , // eps ** -1
                -0.5 , // eps **  0
                -0.25, // eps **  1
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
            if (  target_result_with_prefactor.at(order) != 0.0  )
                REQUIRE(  result_with_prefactor.at(order).uncertainty.real() <= std::abs(2*epsrel * target_result_with_prefactor.at(order))  );

            // check values
            if (  target_result_with_prefactor.at(order) == 0.0  )
                REQUIRE(  result_with_prefactor.at(order).value.real() <= epsabs  );
            else
                REQUIRE(  result_with_prefactor.at(order).value.real() == Approx( target_result_with_prefactor.at(order) ).epsilon( 10.0*epsrel )  );

            std::cout << "----------------" << std::endl << std::endl;
        }

    }

    SECTION("Sidi's integral transform") {

        // User Specified Phase-space point
        const std::vector<INTEGRAL_NAME::real_t> real_parameters = {};
        const std::vector<INTEGRAL_NAME::complex_t> complex_parameters = {};

        // define and configure integrator
        secdecutil::integrators::Qmc<INTEGRAL_NAME::integrand_return_t,
        INTEGRAL_NAME::maximal_number_of_integration_variables,
        ::integrators::transforms::Sidi<3>::type
        #ifdef SECDEC_WITH_CUDA
            ,INTEGRAL_NAME::cuda_together_integrand_t
        #endif
        > integrator;
        integrator.minn = 1e5;
        integrator.maxeval = 1e6;

        const double epsrel = 1e-11;
        const double epsabs = 1e-10;
        integrator.epsrel = 1e-11;
        integrator.epsabs = 1e-10;
        integrator.randomgenerator.seed(798431);

        // Construct the amplitudes
        std::vector<INTEGRAL_NAME::nested_series_t<INTEGRAL_NAME::sum_t>> unwrapped_amplitudes =
            INTEGRAL_NAME::make_amplitudes(real_parameters, complex_parameters, "../one_integration_variable/one_integration_variable_data", integrator);

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
        secdecutil::Series<double> target_result_with_prefactor
        (

            -1, // lowest order in epsilon
             1, // highest computed order in epsilon

            {
                 0.0 , // eps ** -1
                -0.5 , // eps **  0
                -0.25, // eps **  1
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
            if (  target_result_with_prefactor.at(order) != 0.0  )
                REQUIRE(  result_with_prefactor.at(order).uncertainty.real() <= std::abs(2*epsrel * target_result_with_prefactor.at(order))  );

            // check values
            if (  target_result_with_prefactor.at(order) == 0.0  )
                REQUIRE(  result_with_prefactor.at(order).value.real() <= epsabs  );
            else
                REQUIRE(  result_with_prefactor.at(order).value.real() == Approx( target_result_with_prefactor.at(order) ).epsilon( 10.0*epsrel )  );

            std::cout << "----------------" << std::endl << std::endl;
        }

    }

};

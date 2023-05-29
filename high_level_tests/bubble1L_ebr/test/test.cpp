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
#define INTEGRAL_NAME bubble1L_ebr

#include QUOTE_EXPAND(INTEGRAL_NAME.hpp)

TEST_CASE( "check result qmc", "[INTEGRAL_NAME]" ) {

    SECTION("default integral transform") {

        // User Specified Phase-space point
        const std::vector<INTEGRAL_NAME::real_t> real_parameters = { 4., 0.002, 1. };
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
            INTEGRAL_NAME::make_amplitudes(real_parameters, complex_parameters, "../bubble1L_ebr/bubble1L_ebr_data", integrator);

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

            -2, // lowest order in epsilon
             0, // highest computed order in epsilon

            {
                 0.0    + 0.0    * I, // eps ** -2
                 1.0    + 0.0    * I, // eps ** -1
                 4.509087643808781e-02 + 3.138451060936178e+00 * I, // eps **  0 // obtained from high precision run of pysecdec
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
            REQUIRE(      result_with_prefactor.at(order).uncertainty.real() <= std::abs(2*epsrel * target_result_with_prefactor.at(order).real())      );
            if (  target_result_with_prefactor.at(order).imag() != 0.0  )
                REQUIRE(      result_with_prefactor.at(order).uncertainty.imag() <= std::abs(2*epsrel * target_result_with_prefactor.at(order).imag())      );

            // check values
            REQUIRE(  result_with_prefactor.at(order).value.real() == Approx( target_result_with_prefactor.at(order).real() ).epsilon( 10.0*epsrel )  );
            if (  target_result_with_prefactor.at(order).imag() == 0.0  )
                REQUIRE(  result_with_prefactor.at(order).value.imag() <= epsabs  );
            else
                REQUIRE(  result_with_prefactor.at(order).value.imag() == Approx( target_result_with_prefactor.at(order).imag() ).epsilon( 10.0*epsrel )  );

            std::cout << "----------------" << std::endl << std::endl;
        }
    }

};

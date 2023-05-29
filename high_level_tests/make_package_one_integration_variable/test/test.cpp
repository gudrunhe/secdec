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

        // get the integrands
        #ifdef SECDEC_WITH_CUDA
            const auto sector_integrands = INTEGRAL_NAME::make_cuda_integrands(real_parameters, complex_parameters);
        #else
            const auto sector_integrands = INTEGRAL_NAME::make_integrands(real_parameters, complex_parameters);
        #endif

        // add integrands of sectors (together flag)
        const auto all_sectors = std::accumulate(++sector_integrands.begin(), sector_integrands.end(),
        #ifdef SECDEC_WITH_CUDA
            INTEGRAL_NAME::cuda_together_integrand_t()+
        #endif
        *sector_integrands.begin() );

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

        // integrate
        auto result_without_prefactor = secdecutil::deep_apply( all_sectors,  integrator.integrate );
        auto prefactor = INTEGRAL_NAME::prefactor(real_parameters, complex_parameters);
        auto result_with_prefactor = result_without_prefactor * prefactor;


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
                REQUIRE(  result_with_prefactor.at(order).uncertainty <= std::abs(2*epsrel * target_result_with_prefactor.at(order))  );

            // check values
            if (  target_result_with_prefactor.at(order) == 0.0  )
                REQUIRE(  result_with_prefactor.at(order).value <= epsabs  );
            else
                REQUIRE(  result_with_prefactor.at(order).value == Approx( target_result_with_prefactor.at(order) ).epsilon( 10.0*epsrel )  );

            std::cout << "----------------" << std::endl << std::endl;
        }

    }

    SECTION("Sidi's integral transform") {

        // User Specified Phase-space point
        const std::vector<INTEGRAL_NAME::real_t> real_parameters = {};
        const std::vector<INTEGRAL_NAME::complex_t> complex_parameters = {};

        // get the integrands
        #ifdef SECDEC_WITH_CUDA
            const auto sector_integrands = INTEGRAL_NAME::make_cuda_integrands(real_parameters, complex_parameters);
        #else
            const auto sector_integrands = INTEGRAL_NAME::make_integrands(real_parameters, complex_parameters);
        #endif

        // add integrands of sectors (together flag)
        const auto all_sectors = std::accumulate(++sector_integrands.begin(), sector_integrands.end(),
        #ifdef SECDEC_WITH_CUDA
            INTEGRAL_NAME::cuda_together_integrand_t()+
        #endif
        *sector_integrands.begin() );

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

        // integrate
        auto result_without_prefactor = secdecutil::deep_apply( all_sectors,  integrator.integrate );
        auto prefactor = INTEGRAL_NAME::prefactor(real_parameters, complex_parameters);
        auto result_with_prefactor = result_without_prefactor * prefactor;


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
                REQUIRE(  result_with_prefactor.at(order).uncertainty <= std::abs(2*epsrel * target_result_with_prefactor.at(order))  );

            // check values
            if (  target_result_with_prefactor.at(order) == 0.0  )
                REQUIRE(  result_with_prefactor.at(order).value <= epsabs  );
            else
                REQUIRE(  result_with_prefactor.at(order).value == Approx( target_result_with_prefactor.at(order) ).epsilon( 10.0*epsrel )  );

            std::cout << "----------------" << std::endl << std::endl;
        }

    }

};

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
#define INTEGRAL_NAME bubble1L_dot_rank4

#include QUOTE_EXPAND(INTEGRAL_NAME.hpp)


TEST_CASE( "check result qmc", "[INTEGRAL_NAME]" ) {

    SECTION("default integral transform") {

        // User Specified Phase-space point
        const std::vector<INTEGRAL_NAME::real_t> real_parameters = { 1.275, 1.275 };
        const std::vector<INTEGRAL_NAME::complex_t> complex_parameters = { 30.886875, 30.886875, 123.5475 };

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
        integrator.maxeval = 1e7;
        const double epsrel = 1e-4; integrator.epsrel = epsrel;
        const double epsabs = 1e-7; integrator.epsabs = epsabs;

        // integrate
        auto result_without_prefactor = secdecutil::deep_apply( all_sectors,  integrator.integrate );
        auto prefactor = INTEGRAL_NAME::prefactor(real_parameters, complex_parameters);
        auto result_with_prefactor = result_without_prefactor * prefactor;


        // target result, obtained in a long run of pySecDec
        constexpr std::complex<double> I(0,1);
        secdecutil::Series<std::complex<double>> target_result_with_prefactor
        (

            -1, // lowest order in epsilon
             0, // highest computed order in epsilon

            {
                 1.0    + 0.0    * I, // eps ** -1
                -1.2708 + 2.4179 * I, // eps **  0
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

    SECTION("Korobov12 integral transform") {

        // User Specified Phase-space point
        const std::vector<INTEGRAL_NAME::real_t> real_parameters = { 1.275, 1.275 };
        const std::vector<INTEGRAL_NAME::complex_t> complex_parameters = { 30.886875, 30.886875, 123.5475 };

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
        ::integrators::transforms::Korobov<12>::type
        #ifdef SECDEC_WITH_CUDA
            ,INTEGRAL_NAME::cuda_together_integrand_t
        #endif
        > integrator;
        integrator.maxeval = 1e7;
        const double epsrel = 1e-4; integrator.epsrel = epsrel;
        const double epsabs = 1e-7; integrator.epsabs = epsabs;

        // integrate
        auto result_without_prefactor = secdecutil::deep_apply( all_sectors,  integrator.integrate );
        auto prefactor = INTEGRAL_NAME::prefactor(real_parameters, complex_parameters);
        auto result_with_prefactor = result_without_prefactor * prefactor;


        // target result, obtained in a long run of pySecDec
        constexpr std::complex<double> I(0,1);
        secdecutil::Series<std::complex<double>> target_result_with_prefactor
        (

            -1, // lowest order in epsilon
             0, // highest computed order in epsilon

            {
                 1.0    + 0.0    * I, // eps ** -1
                -1.2708 + 2.4179 * I, // eps **  0
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

TEST_CASE( "check result cuhre cquad", "[INTEGRAL_NAME]" ) {

    // User Specified Phase-space point
    const std::vector<INTEGRAL_NAME::real_t> real_parameters = { 1.275, 1.275 };
    const std::vector<INTEGRAL_NAME::complex_t> complex_parameters = { 30.886875, 30.886875, 123.5475 };

    // get the integrands
    const auto sector_integrands = INTEGRAL_NAME::make_integrands(real_parameters, complex_parameters);

    // add integrands of sectors (together flag)
    const auto all_sectors = std::accumulate(++sector_integrands.begin(), sector_integrands.end(), *sector_integrands.begin() );

    // define and configure integrator
    secdecutil::gsl::CQuad<INTEGRAL_NAME::integrand_return_t> cquad;
    secdecutil::cuba::Cuhre<INTEGRAL_NAME::integrand_return_t> cuhre;
    cuhre.maxeval = 1e7;
    cuhre.together = true;
    const double epsrel = 1e-4; cquad.epsrel = cuhre.epsrel = epsrel;
    const double epsabs = 1e-7; cquad.epsabs = cuhre.epsabs = epsabs;
    secdecutil::MultiIntegrator<INTEGRAL_NAME::integrand_return_t,INTEGRAL_NAME::real_t> integrator(cquad,cuhre,2);

    // integrate
    auto result_without_prefactor = secdecutil::deep_apply( all_sectors,  integrator.integrate );
    auto prefactor = INTEGRAL_NAME::prefactor(real_parameters, complex_parameters);
    auto result_with_prefactor = result_without_prefactor * prefactor;


    // target result, obtained in a long run of pySecDec
    constexpr std::complex<double> I(0,1);
    secdecutil::Series<std::complex<double>> target_result_with_prefactor
    (

        -1, // lowest order in epsilon
         0, // highest computed order in epsilon

        {
             1.0    + 0.0    * I, // eps ** -1
            -1.2708 + 2.4179 * I, // eps **  0
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

};

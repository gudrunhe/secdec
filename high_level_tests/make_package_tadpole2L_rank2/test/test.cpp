#include "catch_amalgamated.hpp"
using Catch::Approx;

#include <complex> // std::complex
#include <numeric> // std::accumulate
#include <secdecutil/integrators/integrator.hpp> // MultiIntegrator
#include <secdecutil/integrators/cuba.hpp> // Vegas
#include <secdecutil/integrators/cquad.hpp> // CQuad
#include <secdecutil/deep_apply.hpp> // deep_apply
#include <secdecutil/series.hpp> // Series

#define QUOTE(ARG) #ARG
#define QUOTE_EXPAND(ARG) QUOTE(ARG)
#define INTEGRAL_NAME tadpole2L_rank2

#include QUOTE_EXPAND(INTEGRAL_NAME.hpp)

TEST_CASE( "check result", "[INTEGRAL_NAME]" ) {

    // User Specified Phase-space point
    const std::vector<INTEGRAL_NAME::real_t> real_parameters = { 1., 1., 1. };
    const std::vector<INTEGRAL_NAME::complex_t> complex_parameters = {  };

    // get the integrands
    const auto sector_integrands = INTEGRAL_NAME::make_integrands(real_parameters, complex_parameters);

    // add integrands of sectors (together flag)
    const auto all_sectors = std::accumulate(++sector_integrands.begin(), sector_integrands.end(), *sector_integrands.begin() );

    // define and configure integrator
    secdecutil::gsl::CQuad<INTEGRAL_NAME::integrand_return_t> cquad;
    secdecutil::cuba::Cuhre<INTEGRAL_NAME::integrand_return_t> cuhre;
   /*
    // verbosity
    cuhre.flags = 2;
    cquad.verbose = true;
    */
    cuhre.maxeval = 1e6;
    cuhre.together = true;
    const double epsrel = 1e-8;  cquad.epsrel = cuhre.epsrel = epsrel;
    const double epsabs = 1e-10; cquad.epsabs = cuhre.epsabs = epsabs;
    secdecutil::MultiIntegrator<INTEGRAL_NAME::integrand_return_t,INTEGRAL_NAME::real_t> integrator(cquad,cuhre,2);

    // integrate
    auto result_without_prefactor = secdecutil::deep_apply( all_sectors,  integrator.integrate );
    auto prefactor = INTEGRAL_NAME::prefactor(real_parameters, complex_parameters);


    // target result, obtained in a long run of pySecDec
    constexpr std::complex<double> I(0,1);
    secdecutil::Series<std::complex<double>> target_result_without_prefactor
    (

        -1, // lowest order in epsilon
         1, // highest computed order in epsilon

        {
              10.0,       // eps ** -1
             - 4.0,       // eps **  0
             -34.5127841  // eps **  1
        },

        true, // series is truncated above; i.e. "+ O(eps**2)"
        "eps" // the expansion parameter
    );

    secdecutil::Series<std::complex<double>> target_prefactor
    (

        -1, // lowest order in epsilon
         1, // highest computed order in epsilon

        {
             -0.25,              // eps ** -1
             -0.461392167549234, // eps **  0
             -1.87323249797567   // eps **  1
        },

        true, // series is truncated above; i.e. "+ O(eps**2)"
        "eps" // the expansion parameter
    );

    // basic checks
    REQUIRE(          result_without_prefactor.get_order_min() == target_result_without_prefactor.get_order_min()          );
    REQUIRE(          result_without_prefactor.get_order_max() == target_result_without_prefactor.get_order_max()          );
    REQUIRE(    result_without_prefactor.get_truncated_above() == target_result_without_prefactor.get_truncated_above()    );
    REQUIRE(      result_without_prefactor.expansion_parameter == target_result_without_prefactor.expansion_parameter      );

    std::cout << "----------------" << std::endl << std::endl;

    for (int order = target_result_without_prefactor.get_order_min() ; order <= target_result_without_prefactor.get_order_max() ; ++ order)
    {
        std::cout << "checking order \"eps^" << order << "\" ..." << std::endl;

        // check that the uncertainties are reasonable
        REQUIRE(      result_without_prefactor.at(order).uncertainty.real() <= std::abs(2*epsrel * target_result_without_prefactor.at(order).real())      );
        if (  target_result_without_prefactor.at(order).imag() != 0.0  )
            REQUIRE(      result_without_prefactor.at(order).uncertainty.imag() <= std::abs(2*epsrel * target_result_without_prefactor.at(order).imag())      );

        // check values
        SECTION(" INTEGRAL "){
            REQUIRE(  result_without_prefactor.at(order).value.real() == Approx( target_result_without_prefactor.at(order).real() ).epsilon( 10.0*epsrel )  );
            if (  target_result_without_prefactor.at(order).imag() == 0.0  )
                REQUIRE(  result_without_prefactor.at(order).value.imag() <= epsabs  );
            else
                REQUIRE(  result_without_prefactor.at(order).value.imag() == Approx( target_result_without_prefactor.at(order).imag() ).epsilon( 10.0*epsrel )  );
        }
        SECTION(" PREFACTOR "){
            REQUIRE(  prefactor.at(order).real() == Approx( target_prefactor.at(order).real() ).epsilon( epsrel )  );
            if (  target_prefactor.at(order).imag() == 0.0  )
                REQUIRE(  prefactor.at(order).imag() <= epsabs  );
            else
                REQUIRE(  prefactor.at(order).imag() == Approx( target_prefactor.at(order).imag() ).epsilon( epsrel )  );
        }

        std::cout << "----------------" << std::endl << std::endl;
    }

};

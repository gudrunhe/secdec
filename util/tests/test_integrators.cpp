#include "catch.hpp"
#include "../secdecutil/integrand_container.hpp"
#include "../secdecutil/integrators/integrator.hpp"
#include "../secdecutil/integrators/cuba.hpp"
#include "../secdecutil/uncertainties.hpp"

#include <complex>
#include <cuba.h>
#include <string>

void test_integrator_real(secdecutil::Integrator<cubareal,cubareal>& integrator, cubareal epsrel, int dimensionality = 10){
  
    const std::function<cubareal(cubareal const * const)> integrand =
    [dimensionality] (cubareal const * const variables)
    {
        cubareal out = 1.;
        for (int i=0; i<dimensionality; ++i)
            out *= 6. * variables[i] * (1. - variables[i]);
        return out;
    };

    const auto integrand_container = secdecutil::IntegrandContainer<cubareal, cubareal const * const>(dimensionality,integrand);
    cubareal expected_result = 1.;

    auto computed_result = integrator.integrate(integrand_container);

    REQUIRE( computed_result.value > 0.9 );
    REQUIRE( computed_result.value < 1.1 );
    REQUIRE( computed_result.value == Approx( expected_result ).epsilon( epsrel ) );
    REQUIRE( computed_result.uncertainty <= epsrel );
};

void test_integrator_complex(secdecutil::Integrator<std::complex<cubareal>,cubareal>& integrator, cubareal epsrel){

    constexpr int dimensionality = 4;

    const std::function<std::complex<cubareal>(cubareal const * const)> integrand =
    [] (cubareal const * const variables)
    {
        cubareal out_real = 1.;
        for (int i=0; i<dimensionality; ++i)
            out_real *= 6. * variables[i] * (1. - variables[i]);

        cubareal out_imag = variables[0] * (1. - variables[0]);

        return std::complex<cubareal>{out_real,out_imag};
    };

    const auto integrand_container = secdecutil::IntegrandContainer<std::complex<cubareal>, cubareal const * const>(dimensionality,integrand);
    std::complex<cubareal> expected_result = {1.,1./6.};

    auto computed_result = integrator.integrate(integrand_container);
    
    REQUIRE( computed_result.value.real() > expected_result.real() - 0.1 );
    REQUIRE( computed_result.value.real() < expected_result.real() + 0.1 );

    REQUIRE( computed_result.value.imag() > expected_result.imag() - 0.1 );
    REQUIRE( computed_result.value.imag() < expected_result.imag() + 0.1 );

    REQUIRE( computed_result.value.real() == Approx( expected_result.real() ).epsilon( epsrel ) );
    REQUIRE( computed_result.value.imag() == Approx( expected_result.imag() ).epsilon( epsrel ) );

    REQUIRE( computed_result.uncertainty.real() <= epsrel );
    REQUIRE( computed_result.uncertainty.imag() <= epsrel );
};


TEST_CASE( "Test Vegas integrator with real", "[Integrator][Cuba][Vegas]" ) {
  cubareal epsrel = 1e-3;
  auto integrator = secdecutil::cuba::Vegas<cubareal>(epsrel);
  test_integrator_real(integrator, epsrel);
};

TEST_CASE( "Test Vegas integrator with complex", "[Integrator][Cuba][Vegas]" ) {
  cubareal epsrel = 1e-3;
  auto integrator = secdecutil::cuba::Vegas<std::complex<cubareal>>(epsrel);
  SECTION( "Integrate real and imag together" ) {
    integrator.together = true;
    test_integrator_complex(integrator, epsrel);
  }
  SECTION( "Integrate real and imag separately" ) {
    integrator.together = false;
    test_integrator_complex(integrator, epsrel);
  }
};

TEST_CASE( "Test Suave integrator with real", "[Integrator][Cuba][Suave]" ) {
  cubareal epsrel = 1e-3;
  auto integrator = secdecutil::cuba::Suave<cubareal>(epsrel);
  test_integrator_real(integrator, epsrel);
};

TEST_CASE( "Test Suave integrator with complex", "[Integrator][Cuba][Suave]" ) {
  cubareal epsrel = 1e-3;
  auto integrator = secdecutil::cuba::Suave<std::complex<cubareal>>(epsrel);
  SECTION( "Integrate real and imag together" ) {
    integrator.together = true;
    test_integrator_complex(integrator, epsrel);
  }
  SECTION( "Integrate real and imag separately" ) {
    integrator.together = false;
    test_integrator_complex(integrator, epsrel);
  }
};

TEST_CASE( "Test Divonne integrator with real", "[Integrator][Cuba][Divonne]" ) {
  cubareal epsrel = 1e-3;
  auto integrator = secdecutil::cuba::Divonne<cubareal>(epsrel);
  cubacores(1,1000);
  integrator.key1 = -2000; // Tuned to pass this test
  integrator.seed = 1;
  test_integrator_real(integrator, epsrel);
};

TEST_CASE( "Test Divonne integrator with complex", "[Integrator][Cuba][Divonne]" ) {
  cubareal epsrel = 1e-3;
  auto integrator = secdecutil::cuba::Divonne<std::complex<cubareal>>(epsrel);
  SECTION( "Integrate real and imag together" ) {
    integrator.together = true;
    test_integrator_complex(integrator, epsrel);
  }
  SECTION( "Integrate real and imag separately" ) {
    integrator.together = false;
    test_integrator_complex(integrator, epsrel);
  }
};

 TEST_CASE( "Test Cuhre integrator with real", "[Integrator][Cuba][Cuhre]" ) {
   int dimensionality = 5; // Cuhre can't handle the 10-dim example.
   cubareal epsrel = 1e-3;
   auto integrator = secdecutil::cuba::Cuhre<cubareal>(epsrel);
   integrator.key = 7;
   integrator.maxeval = 1e8;
   test_integrator_real(integrator, epsrel, dimensionality);
 };

TEST_CASE( "Test Cuhre integrator with complex", "[Integrator][Cuba][Cuhre]" ) {
  cubareal epsrel = 1e-3;
  auto integrator = secdecutil::cuba::Cuhre<std::complex<cubareal>>(epsrel);
  SECTION( "Integrate real and imag together" ) {
    integrator.together = true;
    test_integrator_complex(integrator, epsrel);
  }
  SECTION( "Integrate real and imag separately" ) {
    integrator.together = false;
    test_integrator_complex(integrator, epsrel);
  }
 };

 TEST_CASE( "Test exceptions in complex base class", "[Integrator]" ) {

  secdecutil::Integrator<std::complex<cubareal>,cubareal> empty_integrator;
  cubareal epsrel = 1e-4;

  SECTION( "together = true" ) {

    empty_integrator.together = true;
    REQUIRE_THROWS_AS( test_integrator_complex(empty_integrator, epsrel), std::runtime_error );

    try {

      test_integrator_complex(empty_integrator, epsrel);

    } catch (std::runtime_error error) {

      REQUIRE( error.what() == std::string("Simultaneous integration of real and imaginary part is not implemented for this integrator. Try \"together = false\".") );

    }

  }

  SECTION( "together = false" ) {

    empty_integrator.together = false;
    REQUIRE_THROWS_AS( test_integrator_complex(empty_integrator, epsrel), std::runtime_error );

    try {

      test_integrator_complex(empty_integrator, epsrel);

    } catch (std::runtime_error error) {

      REQUIRE( error.what() == std::string("Separate integration of real and imaginary part is not available because pointer to real-valued integrator is not implemented for this integrator. Try \"together = true\".") );

    }

  }

};

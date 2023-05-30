#include "catch_amalgamated.hpp"
using Catch::Approx;

#include "../secdecutil/integrand_container.hpp"
#include "../secdecutil/integrators/integrator.hpp"
#include "../secdecutil/integrators/cuba.hpp"
#include "../secdecutil/uncertainties.hpp"

#include <complex>
#ifdef SECDEC_WITH_CUDA
    #include <thrust/complex.h>
    template <typename ...T> using complex_template = thrust::complex<T...>;
#else
    template <typename ...T> using complex_template = std::complex<T...>;
#endif
#include <cmath>
#include <cuba.h>
#include <string>

template<typename real_t>
void test_integrator_real(secdecutil::Integrator<real_t,real_t>& integrator, cubareal epsrel, int dimensionality = 10, bool test_zero_border = false){

    const std::function<real_t(real_t const * const, secdecutil::ResultInfo*)> integrand =
    [dimensionality] (real_t const * const variables, secdecutil::ResultInfo* result_info)
    {
        real_t out = 1.;
        for (int i=0; i<dimensionality; ++i)
            out *= 6. * variables[i] * (1. - variables[i]);
        return out;
    };

    const auto integrand_container = secdecutil::IntegrandContainer<real_t, real_t const * const>(dimensionality,integrand);
    cubareal expected_result;
    if (test_zero_border)
        expected_result = std::pow(0.5 + 3./4., dimensionality);
    else
        expected_result = 1.;

    auto computed_result = integrator.integrate(integrand_container);

    if (test_zero_border) {
        REQUIRE( computed_result.value > expected_result - 0.1 );
        REQUIRE( computed_result.value < expected_result + 0.1 );
        REQUIRE( computed_result.value == Approx( expected_result ).epsilon( epsrel ) );
        REQUIRE( computed_result.uncertainty <= epsrel * computed_result.value );
    } else {
        REQUIRE( computed_result.value > 0.9 );
        REQUIRE( computed_result.value < 1.1 );
        REQUIRE( computed_result.value == Approx( expected_result ).epsilon( epsrel ) );
        REQUIRE( computed_result.uncertainty <= epsrel * computed_result.value );
    }
};

template<typename real_t>
void test_integrator_complex(secdecutil::Integrator<complex_template<real_t>,real_t>& integrator, cubareal epsrel){

    constexpr int dimensionality = 4;

    const std::function<complex_template<real_t>(real_t const * const, secdecutil::ResultInfo*)> integrand =
    [] (real_t const * const variables, secdecutil::ResultInfo* result_info)
    {
        real_t out_real = 1.;
        for (int i=0; i<dimensionality; ++i)
            out_real *= 6. * variables[i] * (1. - variables[i]);

        real_t out_imag = variables[0] * (1. - variables[0]);

        return complex_template<real_t>{out_real,out_imag};
    };

    const auto integrand_container = secdecutil::IntegrandContainer<complex_template<real_t>, real_t const * const>(dimensionality,integrand);
    complex_template<cubareal> expected_result = {1.,1./6.};

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
  auto integrator = secdecutil::cuba::Vegas<double>(epsrel);
  test_integrator_real(integrator, epsrel);
};

TEST_CASE( "Test Vegas integrator with complex", "[Integrator][Cuba][Vegas]" ) {
  cubareal epsrel = 1e-3;
  auto integrator = secdecutil::cuba::Vegas<complex_template<cubareal>>(epsrel);
  SECTION( "Integrate real and imag together" ) {
    integrator.together = true;
    test_integrator_complex(integrator, epsrel);
  }
  SECTION( "Integrate real and imag separately" ) {
    integrator.together = false;
    test_integrator_complex(integrator, epsrel);
  }
};

TEST_CASE( "Test Vegas integrator with long double", "[Integrator][Cuba][Vegas]" ) {
  cubareal epsrel = 1e-3;
  auto integrator = secdecutil::cuba::Vegas<long double>(epsrel);
  test_integrator_real(integrator, epsrel);
};

TEST_CASE( "Test Vegas integrator with complex long double", "[Integrator][Cuba][Vegas]" ) {
  cubareal epsrel = 1e-3;
  auto integrator = secdecutil::cuba::Vegas<complex_template<long double>>(epsrel);
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
  auto integrator = secdecutil::cuba::Suave<complex_template<cubareal>>(epsrel);
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
  const int num_idle_cores = 1;
  const int num_max_cores = 1000;
  cubacores(&num_idle_cores,&num_max_cores);
  integrator.key1 = -2000; // Tuned to pass this test
  integrator.seed = 1;
  test_integrator_real(integrator, epsrel);
};

TEST_CASE( "Test Divonne integrator with complex", "[Integrator][Cuba][Divonne]" ) {
  cubareal epsrel = 1e-3;
  auto integrator = secdecutil::cuba::Divonne<complex_template<cubareal>>(epsrel);
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
  auto integrator = secdecutil::cuba::Cuhre<complex_template<cubareal>>(epsrel);
  SECTION( "Integrate real and imag together" ) {
    integrator.together = true;
    test_integrator_complex(integrator, epsrel);
  }
  SECTION( "Integrate real and imag separately" ) {
    integrator.together = false;
    test_integrator_complex(integrator, epsrel);
  }
 };

 TEST_CASE( "Test Cuhre result_info contour deformation sign check error", "[Integrator][Cuba][Cuhre]" ) {
    int dimensionality = 1;
    const std::function<double(double const * const, secdecutil::ResultInfo*)> integrand = [] (double const * const variables, secdecutil::ResultInfo* result_info)
        { result_info->return_value = secdecutil::ResultInfo::ReturnValue::sign_check_error_contour_deformation; return 0; };
    const auto integrand_container = secdecutil::IntegrandContainer<double, double const * const>(dimensionality,integrand);

    auto integrator = secdecutil::cuba::Cuhre<cubareal>();

    REQUIRE_THROWS_AS( integrator.integrate(integrand_container) , secdecutil::sign_check_error );
    REQUIRE_THROWS_WITH( integrator.integrate(integrand_container) , Catch::Matchers::ContainsSubstring( "contour deformation" ) );
  };

  TEST_CASE( "Test Vegas result_info positive polynomial sign check error", "[Integrator][Cuba][Cuhre]" ) {
    int dimensionality = 1;
    const std::function<double(double const * const, secdecutil::ResultInfo*)> integrand = [] (double const * const variables, secdecutil::ResultInfo* result_info)
        { result_info->return_value = secdecutil::ResultInfo::ReturnValue::sign_check_error_positive_polynomial; return 0; };
    const auto integrand_container = secdecutil::IntegrandContainer<double, double const * const>(dimensionality,integrand);

    auto integrator = secdecutil::cuba::Vegas<cubareal>();

    REQUIRE_THROWS_AS( integrator.integrate(integrand_container) , secdecutil::sign_check_error );
    REQUIRE_THROWS_WITH( integrator.integrate(integrand_container) , Catch::Matchers::ContainsSubstring( "positive polynomial" ) );
  };

 TEST_CASE( "Test exceptions in complex base class", "[Integrator]" ) {

  secdecutil::Integrator<complex_template<cubareal>,cubareal> empty_integrator;
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

TEST_CASE( "Test integrators one dimensional", "[Integrator][Cuba][Vegas][Suave][Divonne][Cuhre]" ) {

    // some CUBA integrators fail in 1D without an internal workaround
    int dimensionality = 1;
    cubareal epsrel = 1e-3;

    SECTION( "Vegas" ) {

        auto integrator = secdecutil::cuba::Vegas<cubareal>(epsrel);
        test_integrator_real(integrator, epsrel, dimensionality);

    }

    SECTION( "Suave" ) {

        auto integrator = secdecutil::cuba::Suave<cubareal>(epsrel);
        test_integrator_real(integrator, epsrel, dimensionality);

    }

    SECTION( "Divonne" ) {

        auto integrator = secdecutil::cuba::Divonne<cubareal>(epsrel);
        test_integrator_real(integrator, epsrel, dimensionality);

    }

    SECTION( "Cuhre" ) {

        auto integrator = secdecutil::cuba::Cuhre<cubareal>(epsrel);
        test_integrator_real(integrator, epsrel, dimensionality);

    }

};

TEST_CASE( "Test integrators zero_border", "[Integrator][Cuba][Vegas][Suave][Divonne][Cuhre]" ) {

    // some CUBA integrators fail in 1D without an internal workaround
    int dimensionality = 3;
    cubareal epsrel = 1e-4;

    SECTION( "Vegas" ) {

        auto integrator = secdecutil::cuba::Vegas<cubareal>(epsrel);
        integrator.zero_border = 0.5;
        test_integrator_real(integrator, epsrel, dimensionality, /* test_zero_border */ true);

    }

    SECTION( "Suave" ) {

        auto integrator = secdecutil::cuba::Suave<cubareal>(epsrel);
        integrator.zero_border = 0.5; integrator.nnew = 1e4;
        test_integrator_real(integrator, epsrel, dimensionality, /* test_zero_border */ true);

    }

    SECTION( "Divonne" ) {

        auto integrator = secdecutil::cuba::Divonne<cubareal>(epsrel);
        integrator.zero_border = 0.5;
        test_integrator_real(integrator, epsrel, dimensionality, /* test_zero_border */ true);

    }

    SECTION( "Cuhre" ) {

        auto integrator = secdecutil::cuba::Cuhre<cubareal>(epsrel);
        integrator.zero_border = 0.5;
        test_integrator_real(integrator, epsrel, dimensionality, /* test_zero_border */ true);

    }

};

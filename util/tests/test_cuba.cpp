#include "catch.hpp"
#include "../secdecutil/cuba.hpp"
#include "../secdecutil/integrand_container.hpp"
#include "../secdecutil/uncertainties.hpp"

#include <complex>
#include <cuba.h>

TEST_CASE( "Test Vegas integrator with real", "[Cuba][Vegas]" ) {

    constexpr int dimensionality = 10;

    std::function<cubareal(cubareal const * const)> integrand =
    [] (cubareal const * const variables)
    {
        cubareal out = 1.;
        for (int i=0; i<dimensionality; ++i)
            out *= 6. * variables[i] * (1. - variables[i]);
        return out;
    };

    auto integrand_container = secdecutil::IntegrandContainer<cubareal, cubareal const * const>(dimensionality,integrand);
    cubareal expected_result = 1.;

    cubareal epsrel = 1e-3;
    auto integrator = secdecutil::cuba::Vegas<cubareal>(epsrel);
    auto computed_result = integrator.integrate(integrand_container);

    REQUIRE( computed_result.value > 0.9 );
    REQUIRE( computed_result.value < 1.1 );
    REQUIRE( computed_result.value == Approx( expected_result ).epsilon( epsrel ) );
    REQUIRE( computed_result.uncertainty <= epsrel );
};

TEST_CASE( "Test Vegas integrator with complex", "[Cuba][Vegas]" ) {

    constexpr int dimensionality = 4;

    std::function<std::complex<cubareal>(cubareal const * const)> integrand =
    [] (cubareal const * const variables)
    {
        cubareal out_real = 1.;
        for (int i=0; i<dimensionality; ++i)
            out_real *= 6. * variables[i] * (1. - variables[i]);

        cubareal out_imag = variables[0] * (1. - variables[0]);

        return std::complex<cubareal>{out_real,out_imag};
    };

    auto integrand_container = secdecutil::IntegrandContainer<std::complex<cubareal>, cubareal const * const>(dimensionality,integrand);
    std::complex<cubareal> expected_result = {1.,1./6.};

    cubareal epsrel = 1e-3;
    auto integrator = secdecutil::cuba::Vegas<std::complex<cubareal>>(epsrel);
    integrator.maxeval = 1000000;
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

#include "catch_amalgamated.hpp"
using Catch::Approx;

#include <cmath> // std::log
#include <complex> // std::complex
#include <numeric> // std::accumulate
#include <secdecutil/integrators/integrator.hpp> // MultiIntegrator
#include <secdecutil/integrators/cuba.hpp> // Vegas
#include <secdecutil/integrators/cquad.hpp> // CQuad
#include <secdecutil/deep_apply.hpp> // deep_apply
#include <secdecutil/series.hpp> // Series

#define QUOTE(ARG) #ARG
#define QUOTE_EXPAND(ARG) QUOTE(ARG)
#define INTEGRAL_NAME regulator_in_powerlist

#include QUOTE_EXPAND(INTEGRAL_NAME.hpp)

/*
 * constants
 */

const INTEGRAL_NAME::complex_t I(0,1);
const INTEGRAL_NAME::real_t EulerGamma = 0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495L;
const INTEGRAL_NAME::real_t pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068L;

/*
 * target result
 */
// obtained in a long run of pySecDec
const secdecutil::Series<INTEGRAL_NAME::integrand_return_t> target_result_without_prefactor_without_kinematics
{
    -4, // lowest order in epsilon
    -2, // highest computed order in epsilon

    {
          7.0/3.0  , // eps ** -4
          0.0      , // eps ** -3
        -30.7054359, // eps ** -2
    },

    true, // series is truncated above; i.e. "+ O(eps**-1)"
    "eps" // the expansion parameter
};

INTEGRAL_NAME::nested_series_t<INTEGRAL_NAME::integrand_return_t> target_prefactor(INTEGRAL_NAME::real_t s)
{
    /*
     * expansion of "(-s)**(-4*eps - 1)*exp(-I*pi*(2*eps + 3))*gamma(4*eps + 1)/gamma(eps)**2":
     * + eps**2/s
     * + eps**3*(-4*log(-s)/s - 2*EulerGamma/s - 2*I*pi/s)
     * + eps**4*(8*log(-s)**2/s + 8*EulerGamma*log(-s)/s + 8*I*pi*log(-s)/s - 5*pi**2/(6*s) + 2*EulerGamma**2/s + 4*EulerGamma*I*pi/s)
     * + O(eps**5)
    */

    // need the log with nonstandard "-I*delta" prescription
    const INTEGRAL_NAME::integrand_return_t log_minus_s = s<0 ? std::log(-s) : std::log(s) - I*pi;

    return
    {
        2, // lowest order in epsilon
        4, // highest computed order in epsilon

        {
            1./s, // eps ** 2
            -4.*log_minus_s/s - 2.*EulerGamma/s - 2.*I*pi/s, // eps ** 3
            8.*log_minus_s*log_minus_s/s + 8.*EulerGamma*log_minus_s/s + 8.*I*pi*log_minus_s/s - 5.*pi*pi/(6.*s) + 2.*EulerGamma*EulerGamma/s + 4.*EulerGamma*I*pi/s, // eps ** 4
        },

        true, // series is truncated above; i.e. "+ O(eps**5)"
        "eps" // the expansion parameter
    };
}

/*
 * global main test routine
 */
void test_integral(INTEGRAL_NAME::real_t s)
{
    // User Specified Phase-space point
    const std::vector<INTEGRAL_NAME::real_t> real_parameters = { s };
    const std::vector<INTEGRAL_NAME::complex_t> complex_parameters = {  };

    // define and configure integrator
    secdecutil::gsl::CQuad<INTEGRAL_NAME::integrand_return_t> cquad;
    secdecutil::cuba::Cuhre<INTEGRAL_NAME::integrand_return_t> cuhre;
    cuhre.maxeval = 1e7;
    cuhre.together = true;
    const double epsrel = 1e-3; cquad.epsrel = cuhre.epsrel = epsrel;
    const double epsabs = 1e-4; cquad.epsabs = cuhre.epsabs = epsabs;
    cuhre.flags = 2; cquad.verbose = true; // be verbose
    secdecutil::MultiIntegrator<INTEGRAL_NAME::integrand_return_t,INTEGRAL_NAME::real_t> integrator(cquad,cuhre,2);

    // Construct the amplitudes
    std::vector<INTEGRAL_NAME::nested_series_t<INTEGRAL_NAME::sum_t>> unwrapped_amplitudes =
        INTEGRAL_NAME::make_amplitudes(real_parameters, complex_parameters, "../regulator_in_powerlist/regulator_in_powerlist_data", integrator);

    // Pack amplitudes into handler
    INTEGRAL_NAME::handler_t<INTEGRAL_NAME::amplitudes_t> amplitudes
    (
        unwrapped_amplitudes, epsrel, epsabs
        // further optional arguments: epsrel, epsabs, maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
    );
    
    // integrate
    const std::vector<INTEGRAL_NAME::nested_series_t<secdecutil::UncorrelatedDeviation<INTEGRAL_NAME::integrand_return_t>>> result = amplitudes.evaluate();
    auto result_with_prefactor = result.at(0);
    
    // compute target result
    const INTEGRAL_NAME::nested_series_t<INTEGRAL_NAME::integrand_return_t> target_result_with_prefactor =
        target_result_without_prefactor_without_kinematics * target_prefactor(s);

    // basic checks
    CHECK(          result_with_prefactor.get_order_min() == target_result_with_prefactor.get_order_min()          );
    CHECK(          result_with_prefactor.get_order_max() == target_result_with_prefactor.get_order_max()          );
    CHECK(    result_with_prefactor.get_truncated_above() == target_result_with_prefactor.get_truncated_above()    );
    CHECK(      result_with_prefactor.expansion_parameter == target_result_with_prefactor.expansion_parameter      );

    std::cout << "obtained result:" << std::endl << result_with_prefactor << std::endl << std::endl;
    std::cout << "target result:" << std::endl << target_result_with_prefactor << std::endl << std::endl;

    std::cout << "----------------" << std::endl << std::endl;

    for (int order = target_result_with_prefactor.get_order_min() ; order <= target_result_with_prefactor.get_order_max() ; ++order)
    {
        std::cout << "checking order \"eps^" << order << "\" ... with prefactor" << std::endl;

        // check values (`2*epsrel` because the orders are mixed by the prefactor.)
        if (  target_result_with_prefactor.at(order).real() == 0.0  ) {
            CHECK(  result_with_prefactor.at(order).uncertainty.real() <= epsabs  );
            CHECK(  result_with_prefactor.at(order).value.real() <= 2. * result_with_prefactor.at(order).uncertainty.real()  );
        } else
            CHECK(  result_with_prefactor.at(order).value.real() == Approx( target_result_with_prefactor.at(order).real() ).epsilon( 2.0*epsrel )  );
        if (  target_result_with_prefactor.at(order).imag() == 0.0  ) {
            CHECK(  result_with_prefactor.at(order).uncertainty.imag() <= epsabs  );
            CHECK(  result_with_prefactor.at(order).value.imag() <= 2. * result_with_prefactor.at(order).uncertainty.imag()  );
        } else
            CHECK(  result_with_prefactor.at(order).value.imag() == Approx( target_result_with_prefactor.at(order).imag() ).epsilon( 2.0*epsrel )  );
    }

    std::cout << "----------------" << std::endl << std::endl;

};

/*
 * tests
 */

TEST_CASE( "check Euclidean point", "[INTEGRAL_NAME]" ) {
    test_integral( -2.0 );
};

TEST_CASE( "check physical point", "[INTEGRAL_NAME]" ) {
    test_integral( +2.0 );
};

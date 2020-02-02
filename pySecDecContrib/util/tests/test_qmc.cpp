#include "catch.hpp"
#include "../secdecutil/integrand_container.hpp"
#include "../secdecutil/integrators/integrator.hpp"
#include "../secdecutil/integrators/qmc.hpp"
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

  TEST_CASE( "Test result_info contour deformation sign check error with qmc", "[IntegrandContainer]" ) {
    const int dimensionality = 1;
    const std::function<double(double const * const, secdecutil::ResultInfo*)> integrand = [] (double const * const variables, secdecutil::ResultInfo* result_info)
        { result_info->return_value = secdecutil::ResultInfo::ReturnValue::sign_check_error_contour_deformation; return 0; };
    const auto integrand_container = secdecutil::IntegrandContainer<double, double const * const>(dimensionality,integrand);

    auto integrator = secdecutil::integrators::Qmc<double,dimensionality,integrators::transforms::Korobov<3>::type, decltype(integrand_container)>();

    REQUIRE_THROWS_AS( integrator.integrate(integrand_container) , secdecutil::sign_check_error );
    REQUIRE_THROWS_WITH( integrator.integrate(integrand_container) , Catch::Matchers::Contains( "contour deformation" ) );
  };

  TEST_CASE( "Test result_info positive polynomial sign check error with qmc", "[IntegrandContainer]" ) {
    const int dimensionality = 1;
    const std::function<double(double const * const, secdecutil::ResultInfo*)> integrand = [] (double const * const variables, secdecutil::ResultInfo* result_info)
        { result_info->return_value = secdecutil::ResultInfo::ReturnValue::sign_check_error_positive_polynomial; return 0; };
    const auto integrand_container = secdecutil::IntegrandContainer<double, double const * const>(dimensionality,integrand);

    auto integrator = secdecutil::integrators::Qmc<double,dimensionality,integrators::transforms::Korobov<3>::type, decltype(integrand_container)>();

    REQUIRE_THROWS_AS( integrator.integrate(integrand_container) , secdecutil::sign_check_error );
    REQUIRE_THROWS_WITH( integrator.integrate(integrand_container) , Catch::Matchers::Contains( "positive polynomial" ) );
  };

#include "catch_amalgamated.hpp"
using Catch::Approx;

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

#ifdef SECDEC_WITH_CUDA
    #define HOSTDEVICE __host__ __device__
#else
    #define HOSTDEVICE
#endif

struct sign_check_error_integrand_t
{

    const static unsigned number_of_integration_variables = 4;

    HOSTDEVICE double operator()(double const * const variables, secdecutil::ResultInfo* result_info)
    {
        result_info->return_value = secdecutil::ResultInfo::ReturnValue::sign_check_error_contour_deformation;
        return 0;
    };


} sign_check_error_integrand;

struct sign_check_error_contour_deformation_integrand_t
{

    const static unsigned number_of_integration_variables = 4;

    HOSTDEVICE double operator()(double const * const variables, secdecutil::ResultInfo* result_info)
    {
        result_info->return_value = secdecutil::ResultInfo::ReturnValue::sign_check_error_contour_deformation;
        return 0;
    };


} sign_check_error_contour_deformation_integrand;

struct sign_check_error_positive_polynomial_integrand_t
{

    const static unsigned number_of_integration_variables = 4;

    HOSTDEVICE double operator()(double const * const variables, secdecutil::ResultInfo* result_info)
    {
        result_info->return_value = secdecutil::ResultInfo::ReturnValue::sign_check_error_positive_polynomial;
        return 0;
    };


} sign_check_error_positive_polynomial_integrand;


TEST_CASE( "Test result_info contour deformation sign check error with qmc", "[IntegrandContainer]" ) {
  
    using integrand_t = secdecutil::IntegrandContainer</*integrand_return_t*/ double,/*x*/ double const * const,/*parameters*/ double>;
  
    const int dimensionality = 1;
    const integrand_t integrand_container = secdecutil::IntegrandContainer<double, double const * const>(dimensionality,sign_check_error_contour_deformation_integrand);
    auto integrator = secdecutil::integrators::Qmc<double,dimensionality,integrators::transforms::Korobov<3>::type, integrand_t>();

    REQUIRE_THROWS_AS( integrator.integrate(integrand_container) , secdecutil::sign_check_error );
    REQUIRE_THROWS_WITH( integrator.integrate(integrand_container) , Catch::Matchers::ContainsSubstring( "contour deformation" ) );
};

TEST_CASE( "Test result_info positive polynomial sign check error with qmc", "[IntegrandContainer]" ) {
    
    using integrand_t = secdecutil::IntegrandContainer</*integrand_return_t*/ double,/*x*/ double const * const,/*parameters*/ double>;

    const int dimensionality = 1;
    const integrand_t integrand_container = secdecutil::IntegrandContainer<double, double const * const>(dimensionality,sign_check_error_positive_polynomial_integrand);
    auto integrator = secdecutil::integrators::Qmc<double,dimensionality,integrators::transforms::Korobov<3>::type, decltype(integrand_container)>();

    REQUIRE_THROWS_AS( integrator.integrate(integrand_container) , secdecutil::sign_check_error );
    REQUIRE_THROWS_WITH( integrator.integrate(integrand_container) , Catch::Matchers::ContainsSubstring( "positive polynomial" ) );
};


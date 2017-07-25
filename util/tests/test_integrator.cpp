#include "catch.hpp"
#include "../secdecutil/integrand_container.hpp"
#include "../secdecutil/integrators/integrator.hpp"

#include <complex>

/*
 * Integrator that returns a definite value and error.
 * This integrator class is useful to check which
 * instance has been used.
 */
template<typename return_t, typename input_t>
struct DummyIntegrator : secdecutil::Integrator<return_t,input_t>
{
    // dummy integrator: always return this value and error
    return_t return_value;
    return_t return_error;

    // constructor
    DummyIntegrator(return_t return_value, return_t return_error) :
    return_value(return_value), return_error(return_error)
    {};

    std::function<secdecutil::UncorrelatedDeviation<return_t>
      (const secdecutil::IntegrandContainer<return_t, input_t const * const>&)>
      get_integrate()
      {
          const return_t& return_value_copy = return_value, return_error_copy = return_error;
          return [return_value_copy,return_error_copy](const secdecutil::IntegrandContainer<return_t, input_t const * const>& ic)
            {
                return secdecutil::UncorrelatedDeviation<return_t>{return_value_copy,return_error_copy};
            };
      }
};
template<typename return_t, typename input_t>
struct DummyIntegrator<std::complex<return_t>,input_t> : secdecutil::Integrator<std::complex<return_t>,input_t>
{
    // dummy integrator: always return this value and error
    std::complex<return_t> return_value;
    std::complex<return_t> return_error;

    // constructor
    DummyIntegrator(std::complex<return_t> return_value, std::complex<return_t> return_error) :
    return_value(return_value), return_error(return_error)
    { this->together = true; };

    std::function<secdecutil::UncorrelatedDeviation<std::complex<return_t>>
      (const secdecutil::IntegrandContainer<std::complex<return_t>, input_t const * const>&)>
      get_together_integrate()
      {
        return [this] (const secdecutil::IntegrandContainer<std::complex<return_t>, input_t const * const>& integrand_container) {
          return secdecutil::UncorrelatedDeviation<std::complex<return_t>>(return_value,return_error);
        };
    };
};

double real_func(double const * const x) { return *x+2.0; }
std::complex<double> complex_func(float const * const x) { return *x+2.0; }

TEST_CASE( "Test real MultiIntegrator", "[MultiIntegrator]" ) {

    DummyIntegrator<double,double> integrator1(11,1);
    DummyIntegrator<double,double> integrator2(12,2);
    secdecutil::MultiIntegrator<double,double> switching_integrator(integrator1,integrator2,/*critical_dim*/2);

    secdecutil::IntegrandContainer<double,double const * const> dim0(/*dim*/0,real_func);
    secdecutil::IntegrandContainer<double,double const * const> dim1(/*dim*/1,real_func);
    secdecutil::IntegrandContainer<double,double const * const> dim2(/*dim*/2,real_func);
    secdecutil::IntegrandContainer<double,double const * const> dim3(/*dim*/3,real_func);

    // expect integrator1
    secdecutil::UncorrelatedDeviation<double> result0 = switching_integrator.integrate(dim0);
    secdecutil::UncorrelatedDeviation<double> result1 = switching_integrator.integrate(dim1);
    secdecutil::UncorrelatedDeviation<double> expected_result0{11,1};
    secdecutil::UncorrelatedDeviation<double> expected_result1{11,1};

    // expect integrator2
    secdecutil::UncorrelatedDeviation<double> result2 = switching_integrator.integrate(dim2);
    secdecutil::UncorrelatedDeviation<double> result3 = switching_integrator.integrate(dim3);
    secdecutil::UncorrelatedDeviation<double> expected_result2{12,2};
    secdecutil::UncorrelatedDeviation<double> expected_result3{12,2};

    REQUIRE( result0.value == expected_result0.value ); REQUIRE( result0.uncertainty == expected_result0.uncertainty );
    REQUIRE( result1.value == expected_result1.value ); REQUIRE( result1.uncertainty == expected_result1.uncertainty );
    REQUIRE( result2.value == expected_result2.value ); REQUIRE( result2.uncertainty == expected_result2.uncertainty );
    REQUIRE( result3.value == expected_result3.value ); REQUIRE( result3.uncertainty == expected_result3.uncertainty );

}
TEST_CASE( "Test complex MultiIntegrator", "[MultiIntegrator]" ) {

    DummyIntegrator<std::complex<double>,float> integrator1(std::complex<double>{10,11},std::complex<double>{0,1});
    DummyIntegrator<std::complex<double>,float> integrator2(std::complex<double>{12,13},std::complex<double>{2,3});
    secdecutil::MultiIntegrator<std::complex<double>,float> switching_integrator(integrator1,integrator2,/*critical_dim*/2);

    secdecutil::IntegrandContainer<std::complex<double>,float const * const> dim0(/*dim*/0,complex_func);
    secdecutil::IntegrandContainer<std::complex<double>,float const * const> dim1(/*dim*/1,complex_func);
    secdecutil::IntegrandContainer<std::complex<double>,float const * const> dim2(/*dim*/2,complex_func);
    secdecutil::IntegrandContainer<std::complex<double>,float const * const> dim3(/*dim*/3,complex_func);

    // expect integrator1
    secdecutil::UncorrelatedDeviation<std::complex<double>> result0 = switching_integrator.integrate(dim0);
    secdecutil::UncorrelatedDeviation<std::complex<double>> result1 = switching_integrator.integrate(dim1);
    secdecutil::UncorrelatedDeviation<std::complex<double>> expected_result0{{10,11},{0,1}};
    secdecutil::UncorrelatedDeviation<std::complex<double>> expected_result1{{10,11},{0,1}};

    // expect integrator2
    secdecutil::UncorrelatedDeviation<std::complex<double>> result2 = switching_integrator.integrate(dim2);
    secdecutil::UncorrelatedDeviation<std::complex<double>> result3 = switching_integrator.integrate(dim3);
    secdecutil::UncorrelatedDeviation<std::complex<double>> expected_result2{{12,13},{2,3}};
    secdecutil::UncorrelatedDeviation<std::complex<double>> expected_result3{{12,13},{2,3}};

    REQUIRE( result0.value == expected_result0.value ); REQUIRE( result0.uncertainty == expected_result0.uncertainty );
    REQUIRE( result1.value == expected_result1.value ); REQUIRE( result1.uncertainty == expected_result1.uncertainty );
    REQUIRE( result2.value == expected_result2.value ); REQUIRE( result2.uncertainty == expected_result2.uncertainty );
    REQUIRE( result3.value == expected_result3.value ); REQUIRE( result3.uncertainty == expected_result3.uncertainty );

}

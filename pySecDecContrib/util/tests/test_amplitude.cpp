#include "catch_amalgamated.hpp"
#include "../secdecutil/amplitude.hpp"

#include "../secdecutil/integrators/qmc.hpp" // secdecutil::integrators::Qmc
#include "../secdecutil/integrators/cuba.hpp" // namespace secdecutil::cuba
#include "../secdecutil/integrand_container.hpp" // secdecutil::IntegrandContainer
#include "../secdecutil/uncertainties.hpp" // secdecutil::UncorrelatedDeviation

#include <functional> // std::bind, std::placeholders
#include <memory> // std::shared_ptr, std::make_shared
#include <vector> // std::vector

#ifdef SECDEC_WITH_CUDA
    #include <thrust/complex.h>
    #define HOSTDEVICE __host__ __device__
    typedef thrust::complex<double> complex_t;
#else
    #include <complex>
    #define HOSTDEVICE
    typedef std::complex<double> complex_t;
#endif

using namespace Catch::Matchers;

struct simple_integrand_t
{

    const static unsigned number_of_integration_variables = 4;

    HOSTDEVICE double operator()(double const * const x)
    {
        return x[0] * x[1] * x[2] * x[3]*x[3]; // integrates to 1/24 = 0.041666666...
    };

} simple_integrand;

struct other_integrand_t
{

    const static unsigned number_of_integration_variables = 5;

    HOSTDEVICE double operator()(double const * const x)
    {
        return 10. * (1. + x[0] * x[1] * x[2]) * x[3]*exp(x[3]) * x[4]*exp(x[4]); // integrates to 45/4 = 11.25
    };

} other_integrand;

struct complex_integrand_t
{

    const static unsigned number_of_integration_variables = 2;

    HOSTDEVICE complex_t operator()(double const * const x)
    {
        return x[0] * (complex_t(0.,1.) + x[1]); // integrates to 0.25+0.5i
    };

} complex_integrand;

struct complex_const_integrand_t
{

    const static unsigned number_of_integration_variables = 1;

    HOSTDEVICE complex_t operator()(double const * const x)
    {
        return complex_t(1., 0.); 
    };

} complex_const_integrand;

TEST_CASE( "Integration with QmcIntegral", "[Integral][QmcIntegral]" ) {

    using integrand_t = secdecutil::IntegrandContainer</*integrand_return_t*/ double,/*x*/ double const * const,/*parameters*/ double>;
    using integrator_t = secdecutil::integrators::Qmc</*integrand_return_t*/ double,/*maxdim*/4,integrators::transforms::Korobov<3>::type,integrand_t>;
    using integral_t = secdecutil::amplitude::QmcIntegral</*integrand_return_t*/ double,/*real_t*/ double, integrator_t, integrand_t>;
    
    // Set up integrator
    integrator_t integrator;
    integrator.randomgenerator.seed(42546);
    integrator.generatingvectors.clear();
    integrator.generatingvectors[65521] = {1,18303,27193,16899,31463,13841};
    integrator.generatingvectors[131071] = {1,49763,21432,15971,52704,48065};
    integrator.generatingvectors[196597] = {1,72610,13914,40202,16516,29544};
    integrator.generatingvectors[262139] = {1,76811,28708,119567,126364,5581};
    integrator.generatingvectors[327673] = {1,125075,70759,81229,99364,145331};
    
    // Set up integrand
    const integrand_t simple_integrand_container = integrand_t(simple_integrand.number_of_integration_variables, [](double const * const x, secdecutil::ResultInfo * result_info){return simple_integrand(x);});
    
    // Set up integral
    const std::shared_ptr<integrator_t> integrator_ptr = std::make_shared<integrator_t>(integrator);
    const std::shared_ptr<integral_t> integral_ptr = std::make_shared<integral_t>(integrator_ptr, simple_integrand_container);

    SECTION("getter functions before compute()") {

        REQUIRE(integral_ptr->get_number_of_function_evaluations() == 0);
        REQUIRE(integral_ptr->get_next_number_of_function_evaluations() == integrator_ptr->minn);

        REQUIRE_THROWS_AS(integral_ptr->get_integral_result(), secdecutil::amplitude::integral_not_computed_error);
        REQUIRE_THROWS_WITH(integral_ptr->get_integral_result(), Equals("class Integral: get_integral_result called before compute."));

        REQUIRE_THROWS_AS(integral_ptr->get_integration_time(), secdecutil::amplitude::integral_not_computed_error);
        REQUIRE_THROWS_WITH(integral_ptr->get_integration_time(), Equals("class Integral: get_integration_time called before compute."));

        REQUIRE(integral_ptr->get_scaleexpo() > 0.5);

    };

    SECTION("setting next number of function evaluations, compute(), and getter function after compute()") {

        const bool verbose = false;

        integral_ptr->set_next_number_of_function_evaluations(9000); // Note: number must be above qmc->minn
        REQUIRE(integral_ptr->get_next_number_of_function_evaluations() == 9000);

        integral_ptr->compute(verbose);

        // should have gone to the smallest generating vector
        REQUIRE(integral_ptr->get_number_of_function_evaluations() == 65521);
        REQUIRE(integral_ptr->get_next_number_of_function_evaluations() == 65521);

        const double value_first_estimate = integral_ptr->get_integral_result().value;
        const double uncertainty_first_estimate = integral_ptr->get_integral_result().uncertainty;
        REQUIRE(uncertainty_first_estimate < 1e-8);
        REQUIRE(uncertainty_first_estimate > 1e-15);
        REQUIRE_THAT(value_first_estimate, Catch::Matchers::WithinAbs(1./24., 3.*uncertainty_first_estimate)); // expect 3 sigma agreement 99.7% of the time

        auto integration_time = integral_ptr->get_integration_time(); // should not throw
        REQUIRE(typeid(integration_time) == typeid(double));

        REQUIRE(integral_ptr->get_scaleexpo() > 0.5);

        // should not decrease next number of function evaluations
        integral_ptr->set_next_number_of_function_evaluations(10);
        REQUIRE(integral_ptr->get_number_of_function_evaluations() == 65521);
        REQUIRE(integral_ptr->get_next_number_of_function_evaluations() == 65521);

        // should increase number_of_function_evaluations
        integral_ptr->set_next_number_of_function_evaluations(200000);
        REQUIRE(integral_ptr->get_number_of_function_evaluations() == 65521);
        REQUIRE(integral_ptr->get_next_number_of_function_evaluations() == 200000);

        integral_ptr->compute(verbose);

        // should have a smaller error on the integral now
        const double value_second_estimate = integral_ptr->get_integral_result().value;
        const double uncertainty_second_estimate = integral_ptr->get_integral_result().uncertainty;
        REQUIRE(uncertainty_second_estimate < 5e-10);
        REQUIRE_THAT(value_second_estimate, Catch::Matchers::WithinAbs(1./24., 3.*uncertainty_second_estimate)); // expect 3 sigma agreement 99.7% of the time

        // exceeding largest available QMC lattice // (now just computes the largest lattice)
        //integral_ptr->set_next_number_of_function_evaluations(400000);
        //REQUIRE_THROWS_AS(integral_ptr->compute(verbose), std::domain_error);
        //REQUIRE_THROWS_WITH(integral_ptr->compute(verbose), Equals("class QmcIntegral: The requested number_of_function_evaluations (400000) exceeds the largest available lattice (327673)."));

    };

};

TEST_CASE( "Integration with CubaIntegral", "[Integral][CubaIntegral]" ) {
    
    using integrand_t = secdecutil::IntegrandContainer</*integrand_return_t*/ double,/*x*/ double const * const,/*parameters*/ double>;
    using integrator_t = secdecutil::cuba::Vegas<double>;
    using integral_t = secdecutil::amplitude::CubaIntegral</*integrand_return_t*/ double,/*real_t*/ double, integrator_t, integrand_t>;

    integrator_t integrator;
    integrator.mineval = 12345;
    integrator.seed = 12345;
    
    // Set up integrand
    const integrand_t simple_integrand_container = integrand_t(simple_integrand.number_of_integration_variables, [](double const * const x, secdecutil::ResultInfo * result_info){return simple_integrand(x);});
    
    // Set up integral
    const std::shared_ptr<integrator_t> integrator_ptr = std::make_shared<integrator_t>(integrator);
    const std::shared_ptr<integral_t> integral_ptr = std::make_shared<integral_t>(integrator_ptr, simple_integrand_container);


    SECTION("getter functions before compute()") {

        REQUIRE(integral_ptr->get_number_of_function_evaluations() == 0);
        REQUIRE(integral_ptr->get_next_number_of_function_evaluations() == integrator_ptr->mineval);

        REQUIRE_THROWS_AS(integral_ptr->get_integral_result(), secdecutil::amplitude::integral_not_computed_error);
        REQUIRE_THROWS_WITH(integral_ptr->get_integral_result(), Equals("class Integral: get_integral_result called before compute."));

        REQUIRE_THROWS_AS(integral_ptr->get_integration_time(), secdecutil::amplitude::integral_not_computed_error);
        REQUIRE_THROWS_WITH(integral_ptr->get_integration_time(), Equals("class Integral: get_integration_time called before compute."));

        REQUIRE(integral_ptr->get_scaleexpo() == 0.5); // Monte Carlo scaling

    };

    SECTION("setting next number of function evaluations, compute(), and getter function after compute()") {
        
        const bool verbose = false;

        integral_ptr->set_next_number_of_function_evaluations(1e5);
        REQUIRE(integral_ptr->get_next_number_of_function_evaluations() == 1e5);

        integral_ptr->compute(verbose);

        // should have gone to the smallest generating vector
        REQUIRE(integral_ptr->get_number_of_function_evaluations() == 1e5);
        REQUIRE(integral_ptr->get_next_number_of_function_evaluations() == 1e5);

        const double value_first_estimate = integral_ptr->get_integral_result().value;
        const double uncertainty_first_estimate = integral_ptr->get_integral_result().uncertainty;
        REQUIRE(uncertainty_first_estimate < 1e-3);
        REQUIRE(uncertainty_first_estimate > 1e-5);
        REQUIRE_THAT(value_first_estimate, Catch::Matchers::WithinAbs(1./24.,3.*uncertainty_first_estimate)); // expect 3 sigma agreement 99.7% of the time

        auto integration_time = integral_ptr->get_integration_time(); // should not throw
        REQUIRE(typeid(integration_time) == typeid(double));

        REQUIRE(integral_ptr->get_scaleexpo() == 0.5);

        // should not decrease next number of function evaluations
        integral_ptr->set_next_number_of_function_evaluations(10);
        REQUIRE(integral_ptr->get_number_of_function_evaluations() == 1e5);
        REQUIRE(integral_ptr->get_next_number_of_function_evaluations() == 1e5);

        // should increase number_of_function_evaluations
        integral_ptr->set_next_number_of_function_evaluations(5e7);
        REQUIRE(integral_ptr->get_number_of_function_evaluations() == 1e5);
        REQUIRE(integral_ptr->get_next_number_of_function_evaluations() == 5e7);

        integral_ptr->compute(verbose);

        // should have a smaller error on the integral now
        const double value_second_estimate = integral_ptr->get_integral_result().value;
        const double uncertainty_second_estimate = integral_ptr->get_integral_result().uncertainty;
        REQUIRE(uncertainty_second_estimate < 5e-6);
        REQUIRE_THAT(value_second_estimate, Catch::Matchers::WithinAbs(1./24.,3.*uncertainty_second_estimate)); // expect 3 sigma agreement 99.7% of the time

    };

};


TEST_CASE( "Operator overloads of WeightedIntegral", "[WeightedIntegral]" ) {
    
    using integrand_t = secdecutil::IntegrandContainer</*integrand_return_t*/ double,/*x*/ double const * const,/*parameters*/ double>;
    using qmc_integrator_t = secdecutil::integrators::Qmc</*integrand_return_t*/ double,/*maxdim*/4,integrators::transforms::Korobov<3>::type,integrand_t>;
    using cuba_integrator_t = secdecutil::cuba::Vegas<double>;
    using integral_t = secdecutil::amplitude::Integral</*integrand_return_t*/ double,/*real_t*/ double>;
    using qmc_integral_t = secdecutil::amplitude::QmcIntegral</*integrand_return_t*/ double,/*real_t*/ double, qmc_integrator_t, integrand_t>;
    using cuba_integral_t = secdecutil::amplitude::CubaIntegral</*integrand_return_t*/ double,/*real_t*/ double, cuba_integrator_t, integrand_t>;
    
    qmc_integrator_t qmc_integrator;
    qmc_integrator.randomgenerator.seed(42546);
    
    cuba_integrator_t cuba_integrator;
    cuba_integrator.mineval = 12345;
    cuba_integrator.seed = 123143;

    const std::shared_ptr<qmc_integrator_t> qmc_integrator_ptr = std::make_shared<qmc_integrator_t>(qmc_integrator);
    const std::shared_ptr<cuba_integrator_t> cuba_integrator_ptr = std::make_shared<cuba_integrator_t>(cuba_integrator);

    const integrand_t simple_integrand_container = integrand_t(simple_integrand.number_of_integration_variables, [](double const * const x, secdecutil::ResultInfo * result_info){return simple_integrand(x);});
    const integrand_t other_integrand_container = integrand_t(other_integrand.number_of_integration_variables, [](double const * const x, secdecutil::ResultInfo * result_info){return other_integrand(x);});

    std::shared_ptr<integral_t> simple_integral_ptr = std::make_shared<qmc_integral_t>(qmc_integrator_ptr, simple_integrand_container);
    std::shared_ptr<integral_t> other_integral_ptr = std::make_shared<cuba_integral_t>(cuba_integrator_ptr, other_integrand_container);

    std::vector<secdecutil::amplitude::WeightedIntegral<integral_t,/*coefficient_t*/int>> weighted_simple_integral{{simple_integral_ptr,100}};
    std::vector<secdecutil::amplitude::WeightedIntegral<integral_t,/*coefficient_t*/int>> weighted_other_integral{{other_integral_ptr}};

    SECTION(" += ") {

        weighted_simple_integral += weighted_other_integral;
        REQUIRE( weighted_simple_integral.size() == 2 );
        REQUIRE( weighted_simple_integral.at(0).coefficient == 100 );
        REQUIRE( weighted_simple_integral.at(0).integral == simple_integral_ptr );
        REQUIRE( weighted_simple_integral.at(1).coefficient == 1 );
        REQUIRE( weighted_simple_integral.at(1).integral == other_integral_ptr );

    };

    SECTION(" + ") {

        auto sum = weighted_simple_integral + weighted_other_integral;
        REQUIRE( typeid(sum) == typeid(weighted_simple_integral) );
        REQUIRE( sum.size() == 2 );
        REQUIRE( sum.at(0).coefficient == 100 );
        REQUIRE( sum.at(0).integral == simple_integral_ptr );
        REQUIRE( sum.at(1).coefficient == 1 );
        REQUIRE( sum.at(1).integral == other_integral_ptr );

    };

    SECTION(" *= ") {

        weighted_simple_integral *= 5;
        REQUIRE( weighted_simple_integral.size() == 1 );
        REQUIRE( weighted_simple_integral.at(0).coefficient == 500 );
        REQUIRE( weighted_simple_integral.at(0).integral == simple_integral_ptr );

    };

    SECTION(" * ") {

        // coefficient * integral
        auto c_times_i = 80 * weighted_simple_integral;
        REQUIRE( typeid(c_times_i) == typeid(weighted_simple_integral) );
        REQUIRE( c_times_i.size() == 1 );
        REQUIRE( c_times_i.at(0).coefficient == 8000 );
        REQUIRE( c_times_i.at(0).integral == simple_integral_ptr );

        // integral * coefficient
        auto i_times_c = weighted_other_integral * 80;
        REQUIRE( typeid(i_times_c) == typeid(weighted_other_integral) );
        REQUIRE( i_times_c.size() == 1 );
        REQUIRE( i_times_c.at(0).coefficient == 80 );
        REQUIRE( i_times_c.at(0).integral == other_integral_ptr );

    };

};

TEST_CASE( "Optimized integration with WeightedIntegralHandler", "[WeightedIntegralHandler]" ) {
    
    using integrand_t = secdecutil::IntegrandContainer</*integrand_return_t*/ double,/*x*/ double const * const,/*parameters*/ double>;
    using qmc_integrator_t = secdecutil::integrators::Qmc</*integrand_return_t*/ double,/*maxdim*/4,integrators::transforms::Korobov<6>::type,integrand_t>;
    using cuba_integrator_t = secdecutil::cuba::Cuhre<double>;
    using integral_t = secdecutil::amplitude::Integral</*integrand_return_t*/ double,/*real_t*/ double>;
    using qmc_integral_t = secdecutil::amplitude::QmcIntegral</*integrand_return_t*/ double,/*real_t*/ double, qmc_integrator_t, integrand_t>;
    using cuba_integral_t = secdecutil::amplitude::CubaIntegral</*integrand_return_t*/ double,/*real_t*/ double, cuba_integrator_t, integrand_t>;
    
    qmc_integrator_t qmc_integrator;
    qmc_integrator.randomgenerator.seed(42546);
    qmc_integrator.verbosity=3;
    const std::shared_ptr<qmc_integrator_t> qmc_integrator_ptr = std::make_shared<qmc_integrator_t>(qmc_integrator); // this also checks the copy constructor

    
    cuba_integrator_t cuba_integrator;
    cuba_integrator.flags = 2;

    const std::shared_ptr<cuba_integrator_t> cuba_integrator_ptr = std::make_shared<cuba_integrator_t>(cuba_integrator);

    const integrand_t simple_integrand_container = integrand_t(simple_integrand.number_of_integration_variables, [](double const * const x, secdecutil::ResultInfo * result_info){return simple_integrand(x);});
    const integrand_t other_integrand_container = integrand_t(other_integrand.number_of_integration_variables, [](double const * const x, secdecutil::ResultInfo * result_info){return other_integrand(x);});

    std::shared_ptr<integral_t> simple_integral_ptr = std::make_shared<qmc_integral_t>(qmc_integrator_ptr, simple_integrand_container);
    std::shared_ptr<integral_t> other_integral_ptr = std::make_shared<cuba_integral_t>(cuba_integrator_ptr, other_integrand_container);

    using weighted_integral_sum_t = std::vector<secdecutil::amplitude::WeightedIntegral<integral_t,/*coefficient_t*/double>>;
    std::vector<weighted_integral_sum_t> integral_sums
        {
            weighted_integral_sum_t{{simple_integral_ptr, 2.5}},
            weighted_integral_sum_t{{simple_integral_ptr,12.5}} + weighted_integral_sum_t{{other_integral_ptr,1.2}},
            weighted_integral_sum_t{{simple_integral_ptr,12.5}} + weighted_integral_sum_t{{other_integral_ptr,1.2}} + weighted_integral_sum_t{{other_integral_ptr,-1.2}}
        };
    std::vector<double> integral_sum_solutions =
        {
            2.5/24.0,                 // weighted_integral_sum_t{{simple_integral, 2.5}},
            12.5/24.0 + 1.2*45.0/4.0, // weighted_integral_sum_t{{simple_integral,12.5}} + weighted_integral_sum_t{{other_integral,1.2}},
            12.5/24.0                 // weighted_integral_sum_t{{simple_integral,12.5}} + weighted_integral_sum_t{{other_integral,1.2}} + weighted_integral_sum_t{{other_integral,-1.2}}
        };

    SECTION("instantiation and options with container_t=std::vector") {

        using sum_handler_t = secdecutil::amplitude::WeightedIntegralHandler</*integrand_return_t*/ double, /*real_t*/ double, /*coefficient_t*/ double, /*container_t*/ std::vector>;

        sum_handler_t sum_handler
        (
            integral_sums,
            1e-12, // epsrel
            1e-7, // epsabs
            1e6, // maxeval
            1e3, // mineval
            50., // maxincreasefac
            0.1, // min_epsrel
            1e-7, // min_epsabs
            1e-15, // max_epsrel
            1e-18 // max_epsabs
        );

        for(sum_handler_t::sum_t& sum : sum_handler.expression)
        {
            REQUIRE( 1e-12 == sum.epsrel );
            REQUIRE( 1e-7 == sum.epsabs );
            REQUIRE( 1e6 == sum.maxeval );
            REQUIRE( 1e3 == sum.mineval );
            REQUIRE( 50. == sum.maxincreasefac );
            REQUIRE( 0.1 == sum.min_epsrel );
            REQUIRE( 1e-7 == sum.min_epsabs );
            REQUIRE( 1e-15 == sum.max_epsrel );
            REQUIRE( 1e-18 == sum.max_epsabs );
        }

        // set individual field
        sum_handler.expression.at(2).maxincreasefac = 1.4;
        for(size_t i = 0; i < sum_handler.expression.size(); ++i)
            REQUIRE( sum_handler.expression.at(i).maxincreasefac == (i == 2 ? 1.4 : 50.) );

    };

    SECTION("instantiation and options with container_t=secdecutil::Series") {

        using sum_handler_t = secdecutil::amplitude::WeightedIntegralHandler</*integrand_return_t*/ double, /*real_t*/ double, /*coefficient_t*/ double, /*container_t*/ secdecutil::Series>;

        sum_handler_t sum_handler
        (
            secdecutil::Series<weighted_integral_sum_t>(-1,1,integral_sums),
            1e-12, // epsrel
            1e-7, // epsabs
            1e6, // maxeval
            1e2, // mineval
            50., // maxincreasefac
            0.1, // min_epsrel
            1e-7, // min_epsabs
            1e-15, // max_epsrel
            1e-18 // max_epsabs
        );

        for(sum_handler_t::sum_t& sum : sum_handler.expression)
        {
            REQUIRE( 1e-12 == sum.epsrel );
            REQUIRE( 1e-7 == sum.epsabs );
            REQUIRE( 1e6 == sum.maxeval );
            REQUIRE( 1e2 == sum.mineval );
            REQUIRE( 50. == sum.maxincreasefac );
            REQUIRE( 0.1 == sum.min_epsrel );
            REQUIRE( 1e-7 == sum.min_epsabs );
            REQUIRE( 1e-15 == sum.max_epsrel );
            REQUIRE( 1e-18 == sum.max_epsabs );
        }

        // set individual field
        sum_handler.expression.at(0).min_epsabs = 0.14;
        for(int i = sum_handler.expression.get_order_min(); i <= sum_handler.expression.get_order_max(); ++i)
            REQUIRE( sum_handler.expression.at(i).min_epsabs == (i == 0 ? 0.14 : 1e-7) );

    };
    SECTION("computing the amplitude") {

        using sum_handler_t = secdecutil::amplitude::WeightedIntegralHandler</*integrand_return_t*/ double, /*real_t*/ double, /*coefficient_t*/ double, /*container_t*/ std::vector>;

        double epsrel = 1e-12;

        sum_handler_t sum_handler
        (
            integral_sums,
            epsrel,
            1e-20, // epsabs
            1e6, // maxeval
            1e2, // mineval, to ensure Cuhre doesn't reach requested precision in first integrate call
            50., // maxincreasefac
            1e-10, // min_epsrel
            1e-7, // min_epsabs
            1e-15, // max_epsrel
            1e-18 // max_epsabs
        );

        sum_handler.verbose = true;
        auto sum_results = sum_handler.evaluate();

        REQUIRE( typeid(sum_results) == typeid(std::vector<secdecutil::UncorrelatedDeviation<double>>) );
        for(size_t i = 0; i < sum_handler.expression.size(); ++i)
        {
            std::cout << "sum_results[" << i << "] = " << sum_results.at(i) << std::endl;
            REQUIRE_THAT( sum_results.at(i).value, Catch::Matchers::WithinAbs(integral_sum_solutions.at(i), 3.*epsrel) ); // expect 3 sigma agreement 99.7% of the time
        }

    };

    SECTION("amplitude precision with Vegas, together=false") {
        using integrand_t = secdecutil::IntegrandContainer</*integrand_return_t*/ complex_t,/*x*/ double const * const,/*parameters*/ double>;
        using integral_t = secdecutil::amplitude::Integral</*integrand_return_t*/ complex_t,/*real_t*/ double>;
        
        using cuba_integrator_t = secdecutil::cuba::Vegas<complex_t>;
        using cuba_integral_t = secdecutil::amplitude::CubaIntegral</*integrand_return_t*/ complex_t,/*real_t*/ double, cuba_integrator_t, integrand_t>;
        const std::shared_ptr<cuba_integrator_t> cuba_integrator_ptr = std::make_shared<cuba_integrator_t>();
        cuba_integrator_ptr->flags = 2;
        cuba_integrator_ptr->seed = 42546;

        const integrand_t complex_integrand_container = integrand_t(complex_integrand.number_of_integration_variables, [](double const * const x, secdecutil::ResultInfo * result_info){return complex_integrand(x);});
        std::shared_ptr<integral_t> complex_integral_ptr = std::make_shared<cuba_integral_t>(cuba_integrator_ptr, complex_integrand_container);

        using weighted_integral_sum_t = std::vector<secdecutil::amplitude::WeightedIntegral<integral_t,/*coefficient_t*/complex_t>>;
        std::vector<weighted_integral_sum_t> integral_sums
        {
            weighted_integral_sum_t{{complex_integral_ptr, complex_t(1.0,0.)}},
            weighted_integral_sum_t{{complex_integral_ptr, complex_t(0.,1e-3)}}
        };
        using sum_handler_t = secdecutil::amplitude::WeightedIntegralHandler</*integrand_return_t*/ complex_t, /*real_t*/ double, /*coefficient_t*/ complex_t, /*container_t*/ std::vector>;

        double epsrel = 1e-3;

        sum_handler_t sum_handler
        (
            integral_sums,
            epsrel,
            1e-20, // epsabs
            1e6, // maxeval
            1e3, // mineval
            50., // maxincreasefac
            1e-2, // min_epsrel
            1e-2, // min_epsabs
            1e-14, // max_epsrel
            1e-16 // max_epsabs
        );

        sum_handler.verbose = true;
        auto sum_results = sum_handler.evaluate();

        REQUIRE( typeid(sum_results) == typeid(std::vector<secdecutil::UncorrelatedDeviation<complex_t>>) );

        REQUIRE_THAT( sum_results.at(0).value.real(), Catch::Matchers::WithinAbs(0.25   ,3*sum_results.at(0).uncertainty.real()) ); // expect 3 sigma agreement 99.7% of the time
        REQUIRE_THAT( sum_results.at(0).value.imag(), Catch::Matchers::WithinAbs(0.5    ,3*sum_results.at(0).uncertainty.imag()) ); // expect 3 sigma agreement 99.7% of the time
        REQUIRE_THAT( sum_results.at(1).value.real(), Catch::Matchers::WithinAbs(-5e-4  ,3*sum_results.at(1).uncertainty.real()) ); // expect 3 sigma agreement 99.7% of the time
        REQUIRE_THAT( sum_results.at(1).value.imag(), Catch::Matchers::WithinAbs(2.5e-4 ,3*sum_results.at(1).uncertainty.imag()) ); // expect 3 sigma agreement 99.7% of the time

        REQUIRE( std::abs(sum_results.at(0).uncertainty)/std::abs(sum_results.at(0).value) < epsrel );
        REQUIRE( std::abs(sum_results.at(1).uncertainty)/std::abs(sum_results.at(1).value) < epsrel );
    };

    SECTION("amplitude precision with Suave, together=true") {
        using integrand_t = secdecutil::IntegrandContainer</*integrand_return_t*/ complex_t,/*x*/ double const * const,/*parameters*/ double>;
        using integral_t = secdecutil::amplitude::Integral</*integrand_return_t*/ complex_t,/*real_t*/ double>;
        
        using cuba_integrator_t = secdecutil::cuba::Suave<complex_t>;
        using cuba_integral_t = secdecutil::amplitude::CubaIntegral</*integrand_return_t*/ complex_t,/*real_t*/ double, cuba_integrator_t, integrand_t>;
        const std::shared_ptr<cuba_integrator_t> cuba_integrator_ptr = std::make_shared<cuba_integrator_t>();
        cuba_integrator_ptr->flags = 2;
        cuba_integrator_ptr->nmin = 100;
        cuba_integrator_ptr->seed = 42546;
        cuba_integrator_ptr->together = true;

        const integrand_t complex_integrand_container = integrand_t(complex_integrand.number_of_integration_variables, [](double const * const x, secdecutil::ResultInfo * result_info){return complex_integrand(x);});
        std::shared_ptr<integral_t> complex_integral_ptr = std::make_shared<cuba_integral_t>(cuba_integrator_ptr, complex_integrand_container);

        using weighted_integral_sum_t = std::vector<secdecutil::amplitude::WeightedIntegral<integral_t,/*coefficient_t*/complex_t>>;
        std::vector<weighted_integral_sum_t> integral_sums
        {
            weighted_integral_sum_t{{complex_integral_ptr, complex_t(1e-3,0.)}},
            weighted_integral_sum_t{{complex_integral_ptr, complex_t(0.,1)}}
        };
        using sum_handler_t = secdecutil::amplitude::WeightedIntegralHandler</*integrand_return_t*/ complex_t, /*real_t*/ double, /*coefficient_t*/ complex_t, /*container_t*/ std::vector>;

        double epsrel = 3e-3;

        sum_handler_t sum_handler
        (
            integral_sums,
            epsrel,
            1e-20, // epsabs
            1e6, // maxeval
            1e3, // mineval
            50., // maxincreasefac
            1e-2, // min_epsrel
            1e-2, // min_epsabs
            1e-14, // max_epsrel
            1e-16 // max_epsabs
        );

        sum_handler.verbose = true;
        auto sum_results = sum_handler.evaluate();

        REQUIRE( typeid(sum_results) == typeid(std::vector<secdecutil::UncorrelatedDeviation<complex_t>>) );

        REQUIRE_THAT( sum_results.at(0).value.real(), Catch::Matchers::WithinAbs(2.5e-4 , 3*sum_results.at(0).uncertainty.real()) ); // expect 3 sigma agreement 99.7% of the time
        REQUIRE_THAT( sum_results.at(0).value.imag(), Catch::Matchers::WithinAbs(5e-4   , 3*sum_results.at(0).uncertainty.imag()) ); // expect 3 sigma agreement 99.7% of the time
        REQUIRE_THAT( sum_results.at(1).value.real(), Catch::Matchers::WithinAbs(-0.5   , 3*sum_results.at(1).uncertainty.real()) ); // expect 3 sigma agreement 99.7% of the time
        REQUIRE_THAT( sum_results.at(1).value.imag(), Catch::Matchers::WithinAbs( 0.25  , 3*sum_results.at(1).uncertainty.imag()) ); // expect 3 sigma agreement 99.7% of the time

        REQUIRE( std::abs(sum_results.at(0).uncertainty)/std::abs(sum_results.at(0).value) < epsrel );
        REQUIRE( std::abs(sum_results.at(1).uncertainty)/std::abs(sum_results.at(1).value) < epsrel );
    };
    SECTION("amplitude precision with Divonne, together=true") {
        using integrand_t = secdecutil::IntegrandContainer</*integrand_return_t*/ complex_t,/*x*/ double const * const,/*parameters*/ double>;
        using integral_t = secdecutil::amplitude::Integral</*integrand_return_t*/ complex_t,/*real_t*/ double>;
        
        using cuba_integrator_t = secdecutil::cuba::Divonne<complex_t>;
        using cuba_integral_t = secdecutil::amplitude::CubaIntegral</*integrand_return_t*/ complex_t,/*real_t*/ double, cuba_integrator_t, integrand_t>;
        const std::shared_ptr<cuba_integrator_t> cuba_integrator_ptr = std::make_shared<cuba_integrator_t>();
        cuba_integrator_ptr->flags = 2;
        cuba_integrator_ptr->seed = 42546;
        cuba_integrator_ptr->together = true;

        const integrand_t complex_integrand_container = integrand_t(complex_integrand.number_of_integration_variables, [](double const * const x, secdecutil::ResultInfo * result_info){return complex_integrand(x);});
        std::shared_ptr<integral_t> complex_integral_ptr = std::make_shared<cuba_integral_t>(cuba_integrator_ptr, complex_integrand_container);

        using weighted_integral_sum_t = std::vector<secdecutil::amplitude::WeightedIntegral<integral_t,/*coefficient_t*/complex_t>>;
        std::vector<weighted_integral_sum_t> integral_sums
        {
            weighted_integral_sum_t{{complex_integral_ptr, complex_t(1e-3,0.)}},
            weighted_integral_sum_t{{complex_integral_ptr, complex_t(0.,1)}}
        };
        using sum_handler_t = secdecutil::amplitude::WeightedIntegralHandler</*integrand_return_t*/ complex_t, /*real_t*/ double, /*coefficient_t*/ complex_t, /*container_t*/ std::vector>;

        double epsrel = 1e-4;

        sum_handler_t sum_handler
        (
            integral_sums,
            epsrel,
            1e-20, // epsabs
            1e6, // maxeval
            1e3, // mineval
            50., // maxincreasefac
            1e-2, // min_epsrel
            1e-2, // min_epsabs
            1e-14, // max_epsrel
            1e-16 // max_epsabs
        );

        sum_handler.verbose = true;
        auto sum_results = sum_handler.evaluate();

        REQUIRE( typeid(sum_results) == typeid(std::vector<secdecutil::UncorrelatedDeviation<complex_t>>) );

        REQUIRE_THAT( sum_results.at(0).value.real(), Catch::Matchers::WithinAbs(2.5e-4 , 3*sum_results.at(0).uncertainty.real()) ); // expect 3 sigma agreement 99.7% of the time
        REQUIRE_THAT( sum_results.at(0).value.imag(), Catch::Matchers::WithinAbs(5e-4   , 3*sum_results.at(0).uncertainty.imag()) ); // expect 3 sigma agreement 99.7% of the time
        REQUIRE_THAT( sum_results.at(1).value.real(), Catch::Matchers::WithinAbs(-0.5   , 3*sum_results.at(1).uncertainty.real()) ); // expect 3 sigma agreement 99.7% of the time
        REQUIRE_THAT( sum_results.at(1).value.imag(), Catch::Matchers::WithinAbs( 0.25    , 3*sum_results.at(1).uncertainty.imag()) ); // expect 3 sigma agreement 99.7% of the time

        REQUIRE( std::abs(sum_results.at(0).uncertainty)/std::abs(sum_results.at(0).value) < epsrel );
        REQUIRE( std::abs(sum_results.at(1).uncertainty)/std::abs(sum_results.at(1).value) < epsrel );
    };

    SECTION("amplitude precision with different errormodes") {
        using integrand_t = secdecutil::IntegrandContainer</*integrand_return_t*/ complex_t,/*x*/ double const * const,/*parameters*/ double>;
        using integral_t = secdecutil::amplitude::Integral</*integrand_return_t*/ complex_t,/*real_t*/ double>;
        
        using cuba_integrator_t = secdecutil::cuba::Vegas<complex_t>;
        using cuba_integral_t = secdecutil::amplitude::CubaIntegral</*integrand_return_t*/ complex_t,/*real_t*/ double, cuba_integrator_t, integrand_t>;
        const std::shared_ptr<cuba_integrator_t> cuba_integrator_ptr = std::make_shared<cuba_integrator_t>();
        cuba_integrator_ptr->flags = 2;
        cuba_integrator_ptr->seed = 42546;
        cuba_integrator_ptr->together = true;

        const integrand_t complex_integrand_container = integrand_t(complex_integrand.number_of_integration_variables, [](double const * const x, secdecutil::ResultInfo * result_info){return complex_integrand(x);});
        const integrand_t complex_integrand_container2 = integrand_t(complex_const_integrand.number_of_integration_variables, [](double const * const x, secdecutil::ResultInfo * result_info){return complex_const_integrand(x);});
        std::shared_ptr<integral_t> complex_integral_ptr = std::make_shared<cuba_integral_t>(cuba_integrator_ptr, complex_integrand_container);
        std::shared_ptr<integral_t> complex_integral_ptr2 = std::make_shared<cuba_integral_t>(cuba_integrator_ptr, complex_integrand_container2);

        using weighted_integral_sum_t = std::vector<secdecutil::amplitude::WeightedIntegral<integral_t,/*coefficient_t*/complex_t>>;
        std::vector<weighted_integral_sum_t> integral_sums
        {
            weighted_integral_sum_t{{complex_integral_ptr, complex_t(1.,0.)}} + 
            weighted_integral_sum_t{{complex_integral_ptr2, complex_t(0., 10000.)}}
        };
        using sum_handler_t = secdecutil::amplitude::WeightedIntegralHandler</*integrand_return_t*/ complex_t, /*real_t*/ double, /*coefficient_t*/ complex_t, /*container_t*/ std::vector>;

        double epsrel = 1e-3;

        sum_handler_t sum_handler
        (
            integral_sums,
            epsrel,
            1e-20, // epsabs
            1e6, // maxeval
            1e3, // mineval
            50., // maxincreasefac
            1e-2, // min_epsrel
            1e-2, // min_epsabs
            1e-14, // max_epsrel
            1e-16 // max_epsabs
        );

        sum_handler.verbose = true;
        sum_handler.errormode = sum_handler.all;
        //sum_handler.errormode = secdecutil::amplitude::all;
        auto sum_results = sum_handler.evaluate();

        REQUIRE( sum_results.at(0).uncertainty.real()/std::abs(sum_results.at(0).value.real()) < epsrel );
        REQUIRE( sum_results.at(0).uncertainty.imag()/std::abs(sum_results.at(0).value.imag()) < epsrel );

        std::vector<weighted_integral_sum_t> integral_sums2
        {
            weighted_integral_sum_t{{complex_integral_ptr, complex_t(0.,1.)}} + 
            weighted_integral_sum_t{{complex_integral_ptr2, complex_t(1e-4,0.)}}
        };
        sum_handler_t sum_handler2
        (
            integral_sums2,
            epsrel,
            1e-20, // epsabs
            1e6, // maxeval
            1e3, // mineval
            50., // maxincreasefac
            1e-2, // min_epsrel
            1e-2, // min_epsabs
            1e-14, // max_epsrel
            1e-16 // max_epsabs
        );

        sum_handler2.verbose = true;
        sum_handler2.errormode = sum_handler2.largest;
        sum_results = sum_handler2.evaluate();

        REQUIRE( std::max(sum_results.at(0).uncertainty.real(),sum_results.at(0).uncertainty.imag()) / 
                 std::max(std::abs(sum_results.at(0).value.real()),std::abs(sum_results.at(0).value.imag())) < epsrel );
    };

};

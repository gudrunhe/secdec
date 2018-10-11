#include "catch.hpp"
#include "../secdecutil/amplitude.hpp"

#include "../secdecutil/integrators/qmc.hpp" // secdecutil::integrators::Qmc
#include "../secdecutil/integrators/cuba.hpp" // namespace secdecutil::cuba
#include "../secdecutil/integrand_container.hpp" // secdecutil::IntegrandContainer
#include "../secdecutil/uncertainties.hpp" // secdecutil::UncorrelatedDeviation

#include <functional> // std::bind, std::placeholders
#include <memory> // std::shared_ptr, std::make_shared
#include <vector> // std::vector

#ifdef SECDEC_WITH_CUDA
    #define HOSTDEVICE __host__ __device__
#else
    #define HOSTDEVICE
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

    const static unsigned number_of_integration_variables = 3;

    HOSTDEVICE double operator()(double const * const x)
    {
        return 10. * (1. + x[0] * x[1] * x[2]); // integrates to 45/4 = 11.25
    };

} other_integrand;

TEST_CASE( "Integration with QmcIntegral", "[Integral][QmcIntegral]" ) {

    using integrator_t = secdecutil::integrators::Qmc</*integrand_return_t*/double,/*maxdim*/4,::integrators::transforms::Korobov<1>::type,simple_integrand_t>;
    auto integrator = std::make_shared<integrator_t>();
    integrator->randomgenerator.seed(42546);
    integrator->generatingvectors.clear();
    integrator->generatingvectors[65521] = {1,18303,27193,16899,31463,13841};
    integrator->generatingvectors[131071] = {1,49763,21432,15971,52704,48065};
    integrator->generatingvectors[196597] = {1,72610,13914,40202,16516,29544};
    integrator->generatingvectors[262139] = {1,76811,28708,119567,126364,5581};
    integrator->generatingvectors[327673] = {1,125075,70759,81229,99364,145331};

    secdecutil::amplitude::Integral<
                                       double, // integrand_return_t
                                       double // real_t
                                   >&&
        integral =
            secdecutil::amplitude::QmcIntegral<
                                                  double, // integrand_return_t
                                                  double, // real_t
                                                  integrator_t,
                                                  simple_integrand_t
                                              >
                (integrator, simple_integrand);

    SECTION("getter functions before compute()") {

        REQUIRE(integral.get_number_of_function_evaluations() == 0);
        REQUIRE(integral.get_next_number_of_function_evaluations() == integrator->minn);

        REQUIRE_THROWS_AS(integral.get_integral_result(), secdecutil::amplitude::integral_not_computed_error);
        REQUIRE_THROWS_WITH(integral.get_integral_result(), Equals("class Integral: get_integral_result called before compute."));

        REQUIRE_THROWS_AS(integral.get_integration_time(), secdecutil::amplitude::integral_not_computed_error);
        REQUIRE_THROWS_WITH(integral.get_integration_time(), Equals("class Integral: get_integration_time called before compute."));

        REQUIRE(integral.get_scaleexpo() > 0.5);

    };

    SECTION("setting next number of function evaluations, compute(), and getter function after compute()") {

        integral.set_next_number_of_function_evaluations(10);
        REQUIRE(integral.get_next_number_of_function_evaluations() == 10);

        integral.compute();

        // should have gone to the smallest generating vector
        REQUIRE(integral.get_number_of_function_evaluations() == 65521);
        REQUIRE(integral.get_next_number_of_function_evaluations() == 65521);

        const double value_first_estimate = integral.get_integral_result().value;
        const double uncertainty_first_estimate = integral.get_integral_result().uncertainty;
        REQUIRE(uncertainty_first_estimate < 1e-8);
        REQUIRE(uncertainty_first_estimate > 1e-9);
        REQUIRE(value_first_estimate == Approx(1./24.).epsilon(uncertainty_first_estimate));

        auto integration_time = integral.get_integration_time(); // should not throw
        REQUIRE(typeid(integration_time) == typeid(double));

        REQUIRE(integral.get_scaleexpo() > 0.5);

        // should not decrease next number of function evaluations
        integral.set_next_number_of_function_evaluations(10);
        REQUIRE(integral.get_number_of_function_evaluations() == 65521);
        REQUIRE(integral.get_next_number_of_function_evaluations() == 65521);

        // shoud increase number_of_function_evaluations
        integral.set_next_number_of_function_evaluations(200000);
        REQUIRE(integral.get_number_of_function_evaluations() == 65521);
        REQUIRE(integral.get_next_number_of_function_evaluations() == 200000);

        integral.compute();

        // should have a smaller error on the integral now
        const double value_second_estimate = integral.get_integral_result().value;
        const double uncertainty_second_estimate = integral.get_integral_result().uncertainty;
        REQUIRE(uncertainty_second_estimate < 5e-10);
        REQUIRE(value_second_estimate == Approx(1./24.).epsilon(uncertainty_second_estimate));

        // exceeding largest available QMC lattice
        integral.set_next_number_of_function_evaluations(400000);
        REQUIRE_THROWS_AS(integral.compute(), std::domain_error);
        REQUIRE_THROWS_WITH(integral.compute(), Equals("class QmcIntegral: The requested number_of_function_evaluations (400000) exceeds the largest available lattice (327673)."));

    };

};

TEST_CASE( "Integration with CubaIntegral", "[Integral][CubaIntegral]" ) {

    using integrand_t = secdecutil::IntegrandContainer</*integrand_return_t*/double,/*integrand argument type*/double const * const>; // must cast integrand to this type for cuba integrators
    integrand_t integrand(
                             simple_integrand.number_of_integration_variables,
                             std::bind(
                                           &simple_integrand_t::operator(),
                                           simple_integrand,
                                           std::placeholders::_1
                                      )
                         );

    using integrator_t = secdecutil::cuba::Vegas<double>;
    auto integrator = std::make_shared<integrator_t>();
    integrator->mineval = 12345;
    integrator->seed = 123143;

    secdecutil::amplitude::Integral<
                                       double, // integrand_return_t
                                       double // real_t
                                   >&&
        integral =
            secdecutil::amplitude::CubaIntegral<
                                                   double, // integrand_return_t
                                                   double, // real_t
                                                   integrator_t,
                                                   integrand_t
                                               >
                (integrator, integrand);

    SECTION("getter functions before compute()") {

        REQUIRE(integral.get_number_of_function_evaluations() == 0);
        REQUIRE(integral.get_next_number_of_function_evaluations() == integrator->mineval);

        REQUIRE_THROWS_AS(integral.get_integral_result(), secdecutil::amplitude::integral_not_computed_error);
        REQUIRE_THROWS_WITH(integral.get_integral_result(), Equals("class Integral: get_integral_result called before compute."));

        REQUIRE_THROWS_AS(integral.get_integration_time(), secdecutil::amplitude::integral_not_computed_error);
        REQUIRE_THROWS_WITH(integral.get_integration_time(), Equals("class Integral: get_integration_time called before compute."));

        REQUIRE(integral.get_scaleexpo() == 0.5); // Monte Carlo scaling

    };

    SECTION("setting next number of function evaluations, compute(), and getter function after compute()") {

        integral.set_next_number_of_function_evaluations(1e5);
        REQUIRE(integral.get_next_number_of_function_evaluations() == 1e5);

        integral.compute();

        // should have gone to the smallest generating vector
        REQUIRE(integral.get_number_of_function_evaluations() == 1e5);
        REQUIRE(integral.get_next_number_of_function_evaluations() == 1e5);

        const double value_first_estimate = integral.get_integral_result().value;
        const double uncertainty_first_estimate = integral.get_integral_result().uncertainty;
        REQUIRE(uncertainty_first_estimate < 1e-4);
        REQUIRE(uncertainty_first_estimate > 1e-5);
        REQUIRE(value_first_estimate == Approx(1./24.).epsilon(uncertainty_first_estimate));

        auto integration_time = integral.get_integration_time(); // should not throw
        REQUIRE(typeid(integration_time) == typeid(double));

        REQUIRE(integral.get_scaleexpo() == 0.5);

        // should not decrease next number of function evaluations
        integral.set_next_number_of_function_evaluations(10);
        REQUIRE(integral.get_number_of_function_evaluations() == 1e5);
        REQUIRE(integral.get_next_number_of_function_evaluations() == 1e5);

        // shoud increase number_of_function_evaluations
        integral.set_next_number_of_function_evaluations(5e6);
        REQUIRE(integral.get_number_of_function_evaluations() == 1e5);
        REQUIRE(integral.get_next_number_of_function_evaluations() == 5e6);

        integral.compute();

        // should have a smaller error on the integral now
        const double value_second_estimate = integral.get_integral_result().value;
        const double uncertainty_second_estimate = integral.get_integral_result().uncertainty;
        REQUIRE(uncertainty_second_estimate < 5e-6);
        REQUIRE(value_second_estimate == Approx(1./24.).epsilon(uncertainty_second_estimate));

    };

};

TEST_CASE( "Operator overloads of WeightedIntegral", "[WeightedIntegral]" ) {

    using qmc_integrator_t = secdecutil::integrators::Qmc</*integrand_return_t*/double,/*maxdim*/4,::integrators::transforms::Korobov<6>::type,simple_integrand_t>;
    auto qmc_integrator = std::make_shared<qmc_integrator_t>();
    qmc_integrator->randomgenerator.seed(42546);

    using cuba_integrand_t = secdecutil::IntegrandContainer</*integrand_return_t*/double,/*integrand argument type*/double const * const>; // must cast integrand to this type for cuba integrators
    cuba_integrand_t cuba_integrand(
                                       other_integrand.number_of_integration_variables,
                                       std::bind(
                                                     &other_integrand_t::operator(),
                                                     other_integrand,
                                                     std::placeholders::_1
                                                )
                                   );

    using cuba_integrator_t = secdecutil::cuba::Vegas<double>;
    auto cuba_integrator = std::make_shared<cuba_integrator_t>();
    cuba_integrator->mineval = 12345;
    cuba_integrator->seed = 123143;

    using integral_t =
        secdecutil::amplitude::Integral
        <
            double, // integrand_return_t
            double // real_t
        >;

    std::shared_ptr<integral_t> simple_integral =
        std::make_shared
        <
            secdecutil::amplitude::QmcIntegral
            <
                double, // integrand_return_t
                double, // real_t
                qmc_integrator_t,
                simple_integrand_t
            >
        >
        (qmc_integrator, simple_integrand);

    std::shared_ptr<integral_t> other_integral =
        std::make_shared
        <
            secdecutil::amplitude::CubaIntegral
            <
                double, // integrand_return_t
                double, // real_t
                cuba_integrator_t,
                cuba_integrand_t
            >
        >
        (cuba_integrator, cuba_integrand);

    std::vector<secdecutil::amplitude::WeightedIntegral<integral_t,/*coefficient_t*/int>> weighted_simple_integral{{simple_integral,100}};
    std::vector<secdecutil::amplitude::WeightedIntegral<integral_t,/*coefficient_t*/int>> weighted_other_integral{{other_integral}};

    SECTION(" += ") {

        weighted_simple_integral += weighted_other_integral;
        REQUIRE( weighted_simple_integral.size() == 2 );
        REQUIRE( weighted_simple_integral.at(0).coefficient == 100 );
        REQUIRE( weighted_simple_integral.at(0).integral == simple_integral );
        REQUIRE( weighted_simple_integral.at(1).coefficient == 1 );
        REQUIRE( weighted_simple_integral.at(1).integral == other_integral );

    };

    SECTION(" + ") {

        auto sum = weighted_simple_integral + weighted_other_integral;
        REQUIRE( typeid(sum) == typeid(weighted_simple_integral) );
        REQUIRE( sum.size() == 2 );
        REQUIRE( sum.at(0).coefficient == 100 );
        REQUIRE( sum.at(0).integral == simple_integral );
        REQUIRE( sum.at(1).coefficient == 1 );
        REQUIRE( sum.at(1).integral == other_integral );

    };

    SECTION(" *= ") {

        weighted_simple_integral *= 5;
        REQUIRE( weighted_simple_integral.size() == 1 );
        REQUIRE( weighted_simple_integral.at(0).coefficient == 500 );
        REQUIRE( weighted_simple_integral.at(0).integral == simple_integral );

    };

    SECTION(" * ") {

        // coefficient * integral
        auto c_times_i = 80 * weighted_simple_integral;
        REQUIRE( typeid(c_times_i) == typeid(weighted_simple_integral) );
        REQUIRE( c_times_i.size() == 1 );
        REQUIRE( c_times_i.at(0).coefficient == 8000 );
        REQUIRE( c_times_i.at(0).integral == simple_integral );

        // integral * coefficient
        auto i_times_c = weighted_other_integral * 80;
        REQUIRE( typeid(i_times_c) == typeid(weighted_other_integral) );
        REQUIRE( i_times_c.size() == 1 );
        REQUIRE( i_times_c.at(0).coefficient == 80 );
        REQUIRE( i_times_c.at(0).integral == other_integral );

    };

};

TEST_CASE( "Optimized integration with WeightedIntegralHandler", "[WeightedIntegralHandler]" ) {

    using qmc_integrator_t = secdecutil::integrators::Qmc</*integrand_return_t*/double,/*maxdim*/4,::integrators::transforms::Korobov<6>::type,simple_integrand_t>;
    auto qmc_integrator = std::make_shared<qmc_integrator_t>();
    qmc_integrator->randomgenerator.seed(42546);

    using cuba_integrand_t = secdecutil::IntegrandContainer</*integrand_return_t*/double,/*integrand argument type*/double const * const>; // must cast integrand to this type for cuba integrators
    cuba_integrand_t cuba_integrand(
                                       other_integrand.number_of_integration_variables,
                                       std::bind(
                                                     &other_integrand_t::operator(),
                                                     other_integrand,
                                                     std::placeholders::_1
                                                )
                                   );

    using cuba_integrator_t = secdecutil::cuba::Cuhre<double>;
    auto cuba_integrator = std::make_shared<cuba_integrator_t>();

    using integral_t =
        secdecutil::amplitude::Integral
        <
            double, // integrand_return_t
            double // real_t
        >;

    std::shared_ptr<integral_t> simple_integral =
        std::make_shared
        <
            secdecutil::amplitude::QmcIntegral
            <
                double, // integrand_return_t
                double, // real_t
                qmc_integrator_t,
                simple_integrand_t
            >
        >
        (qmc_integrator, simple_integrand);

    std::shared_ptr<integral_t> other_integral =
        std::make_shared
        <
            secdecutil::amplitude::CubaIntegral
            <
                double, // integrand_return_t
                double, // real_t
                cuba_integrator_t,
                cuba_integrand_t
            >
        >
        (cuba_integrator, cuba_integrand);

    using weighted_integral_sum_t = std::vector<secdecutil::amplitude::WeightedIntegral<integral_t,/*coefficient_t*/double>>;
    std::vector<weighted_integral_sum_t> integral_sums
        {
            weighted_integral_sum_t{{simple_integral, 2.5}},
            weighted_integral_sum_t{{simple_integral,12.5}} + weighted_integral_sum_t{{other_integral,1.2}},
            weighted_integral_sum_t{{simple_integral,12.5}} + weighted_integral_sum_t{{other_integral,1.2}} + weighted_integral_sum_t{{other_integral,-1.2}}
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
            1e4, // mineval
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
            REQUIRE( 1e4 == sum.mineval );
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
            1e4, // mineval
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
            REQUIRE( 1e4 == sum.mineval );
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
            1e3, // mineval
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
            REQUIRE( sum_results.at(i).value == Approx(integral_sum_solutions.at(i)).epsilon(epsrel) );
        }

    };

};

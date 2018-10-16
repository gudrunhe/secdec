#ifndef SecDecUtil_weighted_integral_hpp_included
#define SecDecUtil_weighted_integral_hpp_included

#include <algorithm> // std::max, std::min, std::sort
#include <cassert> // assert
#include <chrono> // std::chrono::steady_clock
#include <iostream> // std::cout
#include <limits> // std::numeric_limits
#include <memory> // std::shared_ptr
#include <string> // std::to_string
#include <stdexcept> // std::domain_error, std::logic_error
#include <utility> // std::declval, std::move
#include <vector> // std::vector

#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply
#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation

// wrapped integrators
#include <secdecutil/integrators/cuba.hpp>
#include <secdecutil/integrators/qmc.hpp>
// cannot set the number of points for CQuad --> not suitable for "Integral"

// TODO: documentation

namespace secdecutil {

    namespace amplitude {

        // this exception is thrown when a getter function of Integral is called before the corresponding field has been populated
        struct integral_not_computed_error : public std::logic_error { using std::logic_error::logic_error; };

        template<typename integrand_return_t, typename real_t>
        class Integral
        {
            protected:

                secdecutil::UncorrelatedDeviation<integrand_return_t> integral_result;
                real_t integration_time;

            private:

                unsigned long long int number_of_function_evaluations, next_number_of_function_evaluations;

            public:

                /*
                 * constructor
                 */
                Integral(unsigned long long int next_number_of_function_evaluations = 0) :
                    number_of_function_evaluations(0),
                    next_number_of_function_evaluations(next_number_of_function_evaluations)
                 {};

                /*
                 * destructor
                 */
                virtual ~Integral() = default; // mark the destructor virtual

                /*
                 * setter functions
                 */
                void set_next_number_of_function_evaluations(unsigned long long int new_number_of_function_evaluations)
                {
                    if(new_number_of_function_evaluations > number_of_function_evaluations) {
                        next_number_of_function_evaluations = new_number_of_function_evaluations;
                    } else {
                        next_number_of_function_evaluations = number_of_function_evaluations;
                    }
                };

                /*
                 * getter functions
                 */
                unsigned long long int get_number_of_function_evaluations() const { return number_of_function_evaluations; };
                unsigned long long int get_next_number_of_function_evaluations() const { return next_number_of_function_evaluations; };
                secdecutil::UncorrelatedDeviation<integrand_return_t> get_integral_result() const
                {
                    if(number_of_function_evaluations <= 0)
                        throw integral_not_computed_error("class Integral: get_integral_result called before compute.");
                    return integral_result;
                };
                real_t get_integration_time() const
                {
                    if(number_of_function_evaluations <= 0)
                        throw integral_not_computed_error("class Integral: get_integration_time called before compute.");
                    return integration_time;
                };
                virtual real_t get_scaleexpo() = 0;

                /*
                 * Functions to compute the integral with the given "number_of_function_evaluations".
                 */
                virtual void compute_impl() = 0; // implementation should populate "integral_result"
                void compute()
                {
                    if(next_number_of_function_evaluations > number_of_function_evaluations) {
                        auto start_time = std::chrono::steady_clock::now();
                        compute_impl();
                        auto end_time = std::chrono::steady_clock::now();
                        number_of_function_evaluations = next_number_of_function_evaluations;
                        integration_time = std::chrono::duration<real_t>(end_time - start_time).count();
                    }
                };
        };

        template<typename integrand_return_t, typename real_t, typename integrator_t, typename integrand_t>
        struct QmcIntegral : public Integral<integrand_return_t,real_t>
        {
            std::shared_ptr<integrator_t> qmc;
            integrand_t integrand;
            real_t scaleexpo;

            real_t get_scaleexpo() override { return scaleexpo; }

            /*
             * constructor
             */
            QmcIntegral(const std::shared_ptr<integrator_t>& qmc, const integrand_t& integrand) :
                Integral<integrand_return_t,real_t>(qmc->minn), qmc(qmc), integrand(integrand), scaleexpo(1.0) {};

            void compute_impl() override
            {
                using std::abs;

                unsigned long long int desired_next_n = this->get_next_number_of_function_evaluations();
                unsigned long long int next_n = qmc->get_next_n(desired_next_n);
                if(next_n < desired_next_n)
                    throw std::domain_error("class QmcIntegral: The requested number_of_function_evaluations ("
                        + std::to_string(desired_next_n) + ") exceeds the largest available lattice ("
                        + std::to_string(next_n) +").");
                this->set_next_number_of_function_evaluations( next_n ); // set number of function evluations to the next larger lattice
                qmc->minn = next_n; // set lattice size, ignore random shifts "minm"
                qmc->maxeval = 1; // make "qmc" ignore epsrel and epsabs

                secdecutil::UncorrelatedDeviation<integrand_return_t> new_result = qmc->integrate(this->integrand); // run the numerical integration

                try {
                    integrand_return_t old_error = this->get_integral_result().uncertainty;
                    if(abs(old_error) > abs(new_result.uncertainty))
                        this->integral_result = std::move(new_result);
                } catch (const integral_not_computed_error&) {
                    this->integral_result = std::move(new_result);
                };
            }
        };

        template<typename integrand_return_t, typename real_t, typename integrator_t, typename integrand_t>
        struct CubaIntegral : public Integral<integrand_return_t,real_t>
        {
            std::shared_ptr<integrator_t> integrator;
            integrand_t integrand;
            real_t scaleexpo;

            real_t get_scaleexpo() override { return scaleexpo; }

            /*
             * constructor
             */
            CubaIntegral(const std::shared_ptr<integrator_t>& integrator, const integrand_t& integrand) :
                Integral<integrand_return_t,real_t>(integrator->mineval > 1 ? integrator->mineval : 1000),
                integrator(integrator),
                integrand(integrand),
                scaleexpo(0.5)
                {};

            void compute_impl() override
            {
                unsigned long long int next_n = this->get_next_number_of_function_evaluations();
                integrator->mineval = integrator->maxeval = next_n; // set numbner of sampling points

                this->integral_result = integrator->integrate(this->integrand); // run the numerical integration
            }
        };

        template<typename integral_t, typename coefficient_t>
        struct WeightedIntegral
        {
            std::shared_ptr<integral_t> integral;
            coefficient_t coefficient;

            /*
             * constructor
             */
            WeightedIntegral(
                                const std::shared_ptr<integral_t>& integral,
                                const coefficient_t& coefficient = coefficient_t(1)
                            ) :
            integral(integral), coefficient(coefficient) {};
        };
        /*
         * required operators
         */
        template<typename integral_t, typename coefficient_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>>& operator+=
        (
            std::vector<WeightedIntegral<integral_t,coefficient_t>>& a,
            const std::vector<WeightedIntegral<integral_t,coefficient_t>>& b
        )
        {
            a.insert(a.end(), b.begin(), b.end());
            return a;
        };
        template<typename integral_t, typename coefficient_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>> operator+
        (
            std::vector<WeightedIntegral<integral_t,coefficient_t>> a,
            const std::vector<WeightedIntegral<integral_t,coefficient_t>>& b
        )
        {
            a.insert(a.end(), b.begin(), b.end());
            return a;
        };
        template<typename integral_t, typename coefficient_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>>& operator*=
        (
            std::vector<WeightedIntegral<integral_t,coefficient_t>>& a,
            const coefficient_t& b
        )
        {
            for(WeightedIntegral<integral_t,coefficient_t>& i : a)
                i.coefficient *= b;
            return a;
        };
        template<typename integral_t, typename coefficient_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>> operator*
        (
            std::vector<WeightedIntegral<integral_t,coefficient_t>> a,
            const coefficient_t& b
        )
        {
            for(WeightedIntegral<integral_t,coefficient_t>& i : a)
                i.coefficient *= b;
            return a;
        };
        template<typename integral_t, typename coefficient_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>> operator*
        (
            const coefficient_t& a,
            std::vector<WeightedIntegral<integral_t,coefficient_t>> b
        )
        {
            for(WeightedIntegral<integral_t,coefficient_t>& i : b)
                i.coefficient *= a;
            return b;
        };


        template<typename integrand_return_t, typename real_t, typename coefficient_t, template<typename...> class container_t>
        class WeightedIntegralHandler
        {
            public:

                using integral_t = Integral<integrand_return_t,real_t>;
                using term_t = WeightedIntegral<integral_t,coefficient_t>;
                using sum_return_t = decltype(std::declval<coefficient_t>()*std::declval<integral_t>().get_integral_result());

                struct sum_t
                {
                    std::vector<term_t> summands;

                    // Note: These fields are set in the constructor of WeightedIntegralHandler
                    real_t epsrel;
                    real_t epsabs;
                    unsigned long long int maxeval;
                    unsigned long long int mineval;
                    real_t maxincreasefac;
                    real_t min_epsrel;
                    real_t min_epsabs;
                    real_t max_epsrel;
                    real_t max_epsabs;

                    // construct sum_t from std::vector<term_t>
                    sum_t(const std::vector<term_t>& input) : summands(input) {};
                };

            protected:

                const std::function<sum_t(const std::vector<term_t>&)> convert_to_sum_t =
                    [] (const std::vector<term_t>& input) {  return sum_t(input);  };

            public:

                /*
                 * member fields
                 */
                bool verbose;
                real_t min_decrease_factor;
                real_t soft_wall_clock_limit;
                real_t hard_wall_clock_limit;
                container_t<sum_t> expression;

                /*
                 * constructor
                 */
                WeightedIntegralHandler
                (
                    const container_t<std::vector<term_t>>& expression,
                    real_t epsrel = 1e-2,
                    real_t epsabs = 1e-7,
                    unsigned long long int maxeval = 2147483647, // largest lattice in the default lattices of the QMC
                    unsigned long long int mineval =      50000,
                    real_t maxincreasefac = 20.,
                    real_t min_epsrel = 0.2,
                    real_t min_epsabs = 1e-4,
                    real_t max_epsrel = 1e-14,
                    real_t max_epsabs = 1e-20
                ) :
                    verbose(false),
                    min_decrease_factor(0.9),
                    hard_wall_clock_limit(std::numeric_limits<double>::infinity()),
                    soft_wall_clock_limit(std::numeric_limits<double>::infinity()),
                    expression( deep_apply(expression,convert_to_sum_t) )
                {
                    std::function<void(sum_t&)> set_parameters =
                        [&] (sum_t& sum)
                        {
                            sum.epsrel = epsrel;
                            sum.epsabs = epsabs;
                            sum.maxeval = maxeval;
                            sum.mineval = mineval;
                            sum.maxincreasefac = maxincreasefac,
                            sum.min_epsrel = min_epsrel;
                            sum.min_epsabs = min_epsabs;
                            sum.max_epsrel = max_epsrel;
                            sum.max_epsabs = max_epsabs;
                        };
                    secdecutil::deep_apply(this->expression, set_parameters);
                };

            private:

                decltype(std::chrono::steady_clock::now()) start_time;

            protected:

                /*
                 * evaluate integrals
                 */
                void evaluate_integrals()
                {
                    std::function<void(sum_t&)> compute =
                        [ this ] (sum_t& sum)
                        {
                            for (term_t& term : sum.summands)
                            {
                                const unsigned long long int curr_n = term.integral->get_number_of_function_evaluations();
                                const unsigned long long int next_n = term.integral->get_next_number_of_function_evaluations();
                                secdecutil::UncorrelatedDeviation<integrand_return_t> old_result;
                                if(verbose)
                                    try {
                                        old_result = term.integral->get_integral_result();
                                    } catch (const integral_not_computed_error&) { /* ignore */ };

                                term.integral->compute();

                                if(verbose)
                                {
                                    std::cout << "compute called" << "; sum: " << &sum << ", term: " << &term << ", integral: " << term.integral;
                                    if(next_n > curr_n)
                                        std::cout << ", res: " << old_result << " -> " << term.integral->get_integral_result()
                                                  << ", n: " << curr_n << " -> " << next_n << std::endl;
                                    else
                                        std::cout << ", res: " << term.integral->get_integral_result()
                                                  << ", n: " << next_n << std::endl;
                                }
                            }
                        };
                    secdecutil::deep_apply(expression, compute);
                };

                /*
                 * evaluate expression
                 */
                std::function<sum_return_t(const sum_t&)> compute_sum =
                    [] (const sum_t& sum)
                    {
                        sum_return_t result(0);
                        for (const term_t& term : sum.summands)
                            result += term.coefficient * term.integral->get_integral_result();
                        return result;
                    };
                container_t<sum_return_t> evaluate_expression()
                {
                    return secdecutil::deep_apply(expression, compute_sum);
                };

                /*
                 * estimate total time and try to get below wall_clock_limit
                 */
                void ensure_wall_clock_limit(bool& repeat, std::vector<integral_t*>& integrals)
                {
                    real_t elapsed_time = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - start_time).count();
                    real_t remaining_time = hard_wall_clock_limit - elapsed_time;
                    real_t time_for_next_iteration = 0.;
                    real_t decrease_factor;
                    bool can_improve_in_time = true;

                    for (integral_t* integral : integrals)
                    {
                        const unsigned long long int& curr_n = integral->get_number_of_function_evaluations();
                        const unsigned long long int& next_n = integral->get_next_number_of_function_evaluations();

                        if(next_n > curr_n) {
                            const real_t n_increase_fac = static_cast<real_t>(next_n)/static_cast<real_t>(curr_n);
                            time_for_next_iteration += integral->get_integration_time() * n_increase_fac;
                        }
                    }

                    if(verbose)
                    {
                        std::cout << std::endl;
                        std::cout << "elapsed time: " << elapsed_time << std::endl;
                        std::cout << "estimated time for integrations: " << time_for_next_iteration << std::endl;
                        std::cout << "remaining time: " << remaining_time << std::endl;
                        std::cout << "soft wall clock limit: " << soft_wall_clock_limit << std::endl;
                        std::cout << "hard wall clock limit: " << hard_wall_clock_limit << std::endl;
                    }

                    // decrease number of sampling points if sampling would run out of time
                    while(remaining_time < time_for_next_iteration)
                    {
                        if(verbose)
                            std::cout << "reducing number of samples due to time constraint ..." << std::endl;

                        decrease_factor = remaining_time / time_for_next_iteration;
                        if(decrease_factor > min_decrease_factor)
                            decrease_factor = min_decrease_factor;

                        // reduce runtime by updating with fewer samples
                        time_for_next_iteration = 0.;
                        can_improve_in_time = false;
                        for (integral_t* integral : integrals)
                        {
                            const unsigned long long int& curr_n = integral->get_number_of_function_evaluations();
                            const unsigned long long int& next_n = integral->get_next_number_of_function_evaluations();

                            if(next_n > curr_n) {
                                const real_t n_increase_fac = static_cast<real_t>(next_n)/static_cast<real_t>(curr_n) * decrease_factor;
                                unsigned long long int reduced_next_n = n_increase_fac * curr_n;
                                integral->set_next_number_of_function_evaluations( reduced_next_n );

                                if(reduced_next_n > curr_n) {
                                    can_improve_in_time = true;
                                    time_for_next_iteration += integral->get_integration_time() * n_increase_fac;
                                }
                            }
                        }

                        // stop iterating if soft limit passed
                        if(elapsed_time > soft_wall_clock_limit)
                            repeat = false;
                        else
                            // keep iterating only if improvements are possible within the time limit
                            repeat = can_improve_in_time;

                        if(verbose)
                        {
                            std::cout << "estimated time for integrations (after reduction of samples): " << time_for_next_iteration << std::endl;
                            std::cout << "can improve in time: " << (can_improve_in_time ? "true" : "false") << std::endl;
                        }
                    }

                    if(verbose)
                        std::cout << "run further refinements: " << (repeat ? "true" : "false") << std::endl;

                };

            /*
             * optimize evaluation of integrals to reach the error goals on the sums
             */
            void refine_integrals()
            {
                // ------- definition of helper subroutines - scroll down for the main part -------

                using std::abs; using std::max; using std::min; using std::sqrt; using std::isfinite;

                // end if no integral needs refinement any more
                bool repeat;

                // Ensure each integral is known to at least min_epsrel and min_epsabs
                std::function<void(sum_t&)> ensure_mineval =
                    [ &repeat ] (sum_t& sum)
                    {
                        for (term_t& term : sum.summands)
                            term.integral->set_next_number_of_function_evaluations( sum.mineval );
                    };

                std::function<void(sum_t&)> ensure_min_prec =
                    [ this, &repeat ] (sum_t& sum)
                    {
                        for (term_t& term : sum.summands)
                        {
                            if(term.integral->get_number_of_function_evaluations() < sum.maxeval)
                            {
                                auto result = term.integral->get_integral_result();
                                real_t abserr = abs(result.uncertainty);
                                if(verbose)
                                    std::cout << "sum: " << &sum << ", term: " << &term << ", integral: " << term.integral << ", current integral result: " << result << std::endl;
                                if(abserr > sum.min_epsabs)
                                {
                                    real_t relerr = abserr / abs( result.value );
                                    if(relerr > sum.min_epsrel)
                                    {
                                        unsigned long long int old_n = term.integral->get_number_of_function_evaluations();
                                        real_t new_n = old_n * pow(min(relerr/sum.min_epsrel,abserr/sum.min_epsabs), 1./term.integral->get_scaleexpo());
                                        if(new_n > old_n)
                                        {
                                            repeat = true;
                                            term.integral->set_next_number_of_function_evaluations( new_n > sum.maxeval ? sum.maxeval : static_cast<unsigned long long int>(new_n) );
                                        }
                                    }
                                }
                            }
                        }
                    };

                /*
                 *    (Q)MC scaling: err_i \propto N_i^-scaleexpo  ( N_i \propto t_i )
                 *    minimize   ERR(t1)^2 = err_1(t1)^2+err_2(T-t1)^2 wrt t1, keeping T fixed
                 *    N_i optimal <==> sqrt(t_i)/err_i =: C_i is the same for all i
                 *
                 *    achieve optimal N_i: replace t_i, err_i according to
                 *    t_i   -> t_i   * C_i^(-beta)            * fac^(+0.5/scaleexpo)
                 *    err_i -> err_i * C_i^(beta*scaleexpo)   * fac^(-0.5)
                 *    with beta=2/(2*scaleexpo+1)
                 *
                 *    optimal N_i: N_i -> N_i * C_i^(-beta) * fac^(0.5/scaleexpo)
                 */
                std::function<void(sum_t&)> ensure_error_goal =
                    [ this, &repeat ] (sum_t& sum)
                    {
                        // compute absolute error goal
                        sum_return_t current_sum = compute_sum(sum);
                        integrand_return_t current_sum_value = current_sum.value;
                        integrand_return_t current_sum_error = current_sum.uncertainty;
                        auto abs_error = abs(current_sum_error);
                        auto abs_error_goal = sum.epsrel * abs(current_sum_value);
                        abs_error_goal = max(abs_error_goal, abs_error*sum.epsrel); // If current error larger than current result set goal based on error
                        abs_error_goal = max(abs_error_goal, sum.epsabs); // Do not request an error smaller than epsabs
                        abs_error_goal *= abs_error_goal; // require variance rather than standard deviation

                        if(verbose)
                            std::cout << std::endl << "sum: " << &sum << ", current sum result : " << current_sum << std::endl;

                        // compute the C_i
                        std::vector<real_t> c; c.reserve( sum.summands.size() );
                        for(auto& term : sum.summands)
                        {
                            auto time = term.integral->get_integration_time();
                            auto integral = term.integral->get_integral_result();
                            real_t abserr = abs(integral.uncertainty);
                            real_t relerr = abserr / abs( integral.value );
                            abserr *= abs(term.coefficient); // contribution to total error of the sum

                            if(verbose)
                                std::cout << "sum: " << &sum << ", term: " << &term << ", integral: " << term.integral << ", contribution to sum error: " << abserr << std::endl;

                            if(relerr < sum.max_epsrel || abserr < sum.max_epsabs)
                                c.push_back(  -1  );
                            else
                                c.push_back(  sqrt(time) / abserr  );
                        }
                        assert(c.size() == sum.summands.size());

                        real_t oldfac, fac = 1.;
                        bool firstiter = true;
                        do {
                            oldfac = fac;

                            real_t var_estimate(0);

                            for(size_t i = 0; i < sum.summands.size(); ++i)
                            {
                                real_t& c_i = c.at(i);
                                term_t& term = sum.summands.at(i);
                                auto result = term.integral->get_integral_result();
                                real_t abserr = abs(result.uncertainty) * abs(term.coefficient); // contribution to total error of the sum
                                real_t absvar = abserr * abserr;
                                if(c_i > 0)
                                {
                                    const auto scaleexpo = term.integral->get_scaleexpo();
                                    const auto beta = 2./(2.*scaleexpo+1.);
                                    const auto beta_times_scaleexpo = beta*scaleexpo;

                                    // use all terms to obtain initial estimate of fac
                                    if( firstiter )
                                        var_estimate += absvar * pow(c_i,2*beta_times_scaleexpo);

                                    // improve fac taking into account that error of term doesn't change if n_target < n_current
                                    else
                                        var_estimate += absvar * min(pow(c_i,2*beta_times_scaleexpo)/fac, 1.);
                                } else {
                                    var_estimate += absvar;
                                }
                            }

                            fac *= var_estimate / abs_error_goal; // Factor by which we must rescale to reach absolute error goal

                            firstiter = false;

                        } while(fac < oldfac && fac != 0); // If term.integral->get_number_of_function_evaluations() > optimal N_i, use current n but update var_estimate

                        // Store target in term
                        for(size_t i = 0; i < sum.summands.size(); ++i)
                        {
                            const term_t& term = sum.summands.at(i);
                            const real_t& c_i = c.at(i);

                            if(c_i > 0)
                            {
                                const unsigned long long int& curr_n = term.integral->get_number_of_function_evaluations();
                                const unsigned long long int& next_n = term.integral->get_next_number_of_function_evaluations();
                                const auto scaleexpo = term.integral->get_scaleexpo();
                                const auto beta = 2./(2.*scaleexpo+1.);

                                unsigned long long int proposed_next_n =
                                    max(
                                            next_n,
                                            static_cast<unsigned long long int>( curr_n*pow(c_i,-beta)*pow(fac,0.5/scaleexpo) )
                                       );
                                proposed_next_n = min(sum.maxeval, proposed_next_n);

                                term.integral->set_next_number_of_function_evaluations(proposed_next_n);

                                // run again if the accuracy of an integral is increased
                                if(proposed_next_n > curr_n)
                                    repeat = true;
                            }

                        }
                    };

                // Damp very large increase in number of points
                std::function<void(sum_t&)> ensure_maxincreasefac =
                    [] (sum_t& sum)
                    {
                        for(auto& term: sum.summands)
                        {
                            // Do not increase points by more than a factor of maxincreasefac
                            const real_t curr_n = term.integral->get_number_of_function_evaluations(); // implicit cast to real_t
                            const real_t next_n = term.integral->get_next_number_of_function_evaluations(); // implicit cast to real_t
                            const unsigned long long int max_next_n = curr_n * sum.maxincreasefac; // implicit cast to unsigned long long int

                            if( next_n > max_next_n )
                                term.integral->set_next_number_of_function_evaluations( max_next_n );
                        }
                    };

                // ------------------------------- main part of the function starts here -------------------------------

                // make a unique vector of the appearing integrals
                std::vector<integral_t*> integrals;
                std::function<void(sum_t&)> populate_integrals =
                    [ &integrals ] (sum_t& sum)
                    {
                        for (term_t& term : sum.summands)
                            integrals.push_back(term.integral.get());
                    };
                secdecutil::deep_apply(expression, populate_integrals);
                std::sort(integrals.begin(), integrals.end());
                integrals.erase(std::unique(integrals.begin(), integrals.end()), integrals.end());
                integrals.shrink_to_fit();

                // initialize with minimal number of sampling points
                secdecutil::deep_apply(expression, ensure_mineval);
                if(verbose)
                    std::cout << "computing integrals to satisfy mineval" << std::endl;
                evaluate_integrals();
                if(verbose)
                    std::cout << "---------------------" << std::endl << std::endl << std::endl;

                // ensure each integral is at least known to min_epsrel and min_epsabs
                do {
                    repeat = false;

                    secdecutil::deep_apply(expression, ensure_min_prec);
                    secdecutil::deep_apply(expression, ensure_maxincreasefac);

                    if(verbose)
                        std::cout << std::endl << "computing integrals to satisfy min_epsrel and min_epsabs" << std::endl;
                    evaluate_integrals();
                    if(verbose)
                        std::cout << "---------------------" << std::endl << std::endl << std::endl;
                } while(repeat);

                // ensure the error goal of each sum is reached as good as possible under the given time constraint
                do {
                    repeat = false;

                    secdecutil::deep_apply(expression, ensure_error_goal);
                    secdecutil::deep_apply(expression, ensure_maxincreasefac);
                    ensure_wall_clock_limit(repeat, integrals);

                    if(verbose)
                        std::cout << std::endl << "computing integrals to satisfy error goals on sums" << std::endl;
                    evaluate_integrals();
                    if(verbose)
                        std::cout << "---------------------" << std::endl << std::endl << std::endl;
                } while(repeat);

            }

            public:

                /*
                 * evaluate the sum to the requested precision
                 */
                container_t<sum_return_t> evaluate()
                {
                    start_time = std::chrono::steady_clock::now();
                    refine_integrals();
                    return evaluate_expression();
                }

        };

    };

};

#endif

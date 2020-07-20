#ifndef SecDecUtil_weighted_integral_hpp_included
#define SecDecUtil_weighted_integral_hpp_included

#include <algorithm> // std::max, std::min, std::sort
#include <cassert> // assert
#include <chrono> // std::chrono::steady_clock
#include <iostream> // std::cout, std::dec
#include <limits> // std::numeric_limits
#include <memory> // std::shared_ptr
#include <string> // std::to_string
#include <stdexcept> // std::domain_error, std::logic_error, std::runtime_error
#include <thread> // std::thread
#include <utility> // std::declval, std::move
#include <vector> // std::vector
#include <fstream> // for writing changed deformation parameters to file

#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply
#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation
#include <secdecutil/sector_container.hpp> // for secdecutil::sign_check_error

// wrapped integrators
#include <secdecutil/integrators/cuba.hpp>
#include <secdecutil/integrators/qmc.hpp>
// cannot set the number of points for CQuad --> not suitable for "Integral"

// TODO: documentation

namespace secdecutil {

    namespace amplitude {

        // this exception is thrown when a getter function of Integral is called before the corresponding field has been populated
        struct integral_not_computed_error : public std::logic_error { using std::logic_error::logic_error; };

        #ifdef SECDEC_WITH_CUDA
            // cuda error handling
            struct cuda_error : public std::runtime_error { using std::runtime_error::runtime_error; };
            inline void cuda_safe_call(cudaError_t error)
            {
                if (error != cudaSuccess)
                    throw cuda_error( cudaGetErrorString(error) );
            };
        #endif

        template<typename integrand_return_t, typename real_t>
        class Integral
        {
            protected:

                secdecutil::UncorrelatedDeviation<integrand_return_t> integral_result;
                real_t integration_time;

            private:

                unsigned long long int number_of_function_evaluations, next_number_of_function_evaluations;

            public:
                typedef real_t real_t_type;
                std::string display_name = "INTEGRAL"; // used as an user readable name for the integral
                int id; // used for user output, so that the integral can be identified more easier in thousands of integrals, without going through the names

                /*
                 * constructor
                 */
                Integral(unsigned long long int next_number_of_function_evaluations = 0) :
                    number_of_function_evaluations(0),
                    next_number_of_function_evaluations(next_number_of_function_evaluations), id(0)
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
                virtual real_t get_scaleexpo() const = 0;
                virtual std::vector<std::vector<real_t*>> get_parameters() = 0;
                virtual std::vector<std::vector<real_t>> get_extra_parameters() = 0;
                virtual void clear_errors() = 0;

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

            real_t get_scaleexpo() const override { return scaleexpo; }

            std::vector<std::vector<real_t*>> get_parameters() override {return integrand.get_parameters();}
            std::vector<std::vector<real_t>> get_extra_parameters() override {return integrand.get_extra_parameters();}
            void clear_errors() override {integrand.clear_errors();}

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

            real_t get_scaleexpo() const override { return scaleexpo; }

            std::vector<std::vector<real_t*>> get_parameters() override {return integrand.get_parameters();}
            std::vector<std::vector<real_t>> get_extra_parameters() override {return integrand.get_extra_parameters();}
            void clear_errors() override {integrand.clear_errors();}

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
            std::string display_name = "WINTEGRAL";

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

        void write_map_to_file(std::map<std::string, std::vector<std::vector<double>>> map, std::string filename="changed_deformation_parameters.txt"){
            std::ofstream file;
            file.open(filename);
            file << std::setprecision(15);
            for(auto el : map){
                file << el.first << ": ";
                for(int k = 0; k <el.second.size(); k++){
                    for(int i = 0; i < el.second[k].size(); i++){
                        file << el.second[k][i];
                        if(i != el.second[k].size()-1)
                            file << ", ";
                    }
                    if(k != el.second.size()-1)
                        file << "; ";
                    else
                        file << std::endl;
                }
            }
            file.close();
        }

        std::map<std::string, std::vector<std::vector<double>>> read_map_from_file(std::string filename="changed_deformation_parameters.txt"){
            std::map<std::string, std::vector<std::vector<double>>> map;
            std::ifstream file;
            file.open(filename);
            for(std::string line; std::getline(file,line);){
                auto integral_name = line.substr(0,line.find(":"));
                auto parameters_string = line.substr(line.find(":")+1,line.size());
                std::vector<std::vector<double>> parameters_vector;
                while(parameters_string.size()){
                    parameters_vector.push_back({});
                    auto semicolon_location = std::min(parameters_string.find(";"),parameters_string.size());
                    auto parameters_k_string = parameters_string.substr(0,semicolon_location);
                    parameters_string = parameters_string.substr(std::min(semicolon_location+1,parameters_string.size()),parameters_string.size());
                    while(parameters_k_string.size()){
                        auto comma_location = std::min(parameters_k_string.find(","),parameters_k_string.size());
                        auto parameter = parameters_k_string.substr(0, comma_location);
                        parameters_k_string = parameters_k_string.substr(std::min(comma_location+1,parameters_k_string.size()), parameters_k_string.size());
                        double parameter_value = std::stod(parameter);
                        parameters_vector.back().push_back(parameter_value);
                    }
                }
                map[integral_name] = parameters_vector;
            }
            file.close();
            return map;
        }

        void print_datetime(std::string prefix = "Current time: "){
            auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << prefix << std::ctime(&t);
        }

        /*
         * evaluate a vector of integrals
         */
        template<typename integral_t>
        void evaluate_integrals(std::vector<integral_t*>& integrals, const bool& verbose, size_t number_of_threads, size_t reset_cuda_after,
                    std::map<std::string, std::vector<std::vector<double>>> changed_deformation_parameters_map)
        {
            if(number_of_threads == 0)
                ++number_of_threads;

            auto compute_integral = [ &verbose, &integrals, &changed_deformation_parameters_map ] (integral_t* integral)
                {
                    const unsigned long long int curr_n = integral->get_number_of_function_evaluations();
                    const unsigned long long int next_n = integral->get_next_number_of_function_evaluations();
                    /*secdecutil::UncorrelatedDeviation<integrand_return_t>*/ decltype(integral->get_integral_result()) old_result;
                    if(verbose)
                        try {
                            old_result = integral->get_integral_result();
                        } catch (const integral_not_computed_error&) { /* ignore */ };

                    bool failed_atleast_once = false;
                    while(true){
                        bool failed = false;
                        try{
                            integral->compute();
                        } catch(secdecutil::sign_check_error& e){
                            failed = true;
                            std::cout << "Exception: " << e.what() << std::endl;
                            integral->clear_errors();

                            std::cout << "Integral " << integral->display_name << " failed, reducing deformation parameters." << std::endl;
                            std::vector<std::vector<typename std::remove_pointer<decltype(integral)>::type::real_t_type*>> pars = integral->get_parameters();
                            std::vector<std::vector<typename std::remove_pointer<decltype(integral)>::type::real_t_type>> extra_pars = integral->get_extra_parameters();
                            failed_atleast_once = true;
                            bool changed_deformation_parameters = false;
                            for(int k = 0; k < pars.size(); k++){
                                for(int i = 0; i < pars[k].size(); i++){
                                    std::cout << "par " << k << "," << i << ": " << *pars[k][i];
                                    if(*pars[k][i] > extra_pars[k][0]){
                                        *pars[k][i] *= extra_pars[k][1];
                                        if(*pars[k][i] < extra_pars[k][0]){
                                            *pars[k][i] = extra_pars[k][0];
                                        }
                                        changed_deformation_parameters = true;
                                        std::cout << " -> " << *pars[k][i];
                                    }
                                    if(i == pars[k].size()-1 and k == pars.size()-1)
                                        std::cout << std::endl;
                                    else
                                        std::cout << ", ";
                                }
                            }
                            if(not changed_deformation_parameters){
                                throw std::runtime_error("All deformation parameters at minimum already, integral still fails.");
                            }
                        }
                        if(not failed){
                            if(failed_atleast_once){
                                std::vector<std::vector<typename std::remove_pointer<decltype(integral)>::type::real_t_type*>> pars = integral->get_parameters();
                                std::vector<std::vector<typename std::remove_pointer<decltype(integral)>::type::real_t_type>> pars_data(pars.size());
                                for(int i = 0; i < pars.size(); i++){
                                    for(auto par:pars[i]){
                                        pars_data[i].push_back(*par);
                                    }
                                }
                                changed_deformation_parameters_map[integral->display_name] = pars_data;
                                write_map_to_file(changed_deformation_parameters_map);
                            }
                            break;
                        }
                    }

                    if(verbose)
                    {
                        std::cout << "integral " << integral->id << "/" << integrals.size() << ": " << integral->display_name << ", time: ";
                        std::printf("%.1f", integral->get_integration_time());
                        std::cout << " s, ";
                        print_datetime();
                        if(next_n > curr_n){
                            std::cout << "res: " << old_result << " -> " << integral->get_integral_result()
                                      << ", n: " << curr_n << " -> " << std::dec << integral->get_number_of_function_evaluations() << std::endl;
                        }
                        else{
                            std::cout << "res: " << integral->get_integral_result()
                                      << ", n: " << std::dec << next_n << std::endl;
                        }
                        std::cout << std::endl;
                    }
                };

            std::vector<std::thread> thread_pool(number_of_threads);
            size_t idx = 0;
            #ifdef SECDEC_WITH_CUDA
                size_t job_counter = 0;
                int number_of_cuda_devices; cuda_safe_call( cudaGetDeviceCount(&number_of_cuda_devices) );
            #endif
            for (integral_t* integral : integrals)
            {
                if(thread_pool.at(idx).joinable())
                    thread_pool.at(idx).join();

                thread_pool.at(idx) = std::thread(compute_integral,integral);

                // reset cuda devices after running "reset_cuda_after" integrations
                #ifdef SECDEC_WITH_CUDA
                    if(reset_cuda_after > 0 && ++job_counter >= reset_cuda_after)
                    {
                        job_counter = 0;
                        for(std::thread& worker : thread_pool)
                            if(worker.joinable())
                                worker.join();
                        for(int device_id = 0; device_id < number_of_cuda_devices; ++device_id)
                        {
                            cuda_safe_call( cudaSetDevice(device_id) );
                            cuda_safe_call( cudaDeviceReset() );
                        }
                    }
                #endif

                if(++idx >= number_of_threads)
                    idx = 0;
            }

            for(std::thread& worker : thread_pool)
                if(worker.joinable())
                    worker.join();
        }
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
                    std::string display_name = "SUM";

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

                std::map<std::string, std::vector<std::vector<double>>> changed_deformation_parameters_map;

                // these functions give names to the sum_t sums containing the orders of the series
                void name_sum(sum_t& sum, std::string prefix){
                    sum.display_name = prefix;
                }

                template<typename T>
                void name_sum(Series<T>& series, std::string prefix){
                    for(int i = series.get_order_min(); i <= series.get_order_max(); i++){
                        name_sum(series.at(i), prefix + "_" + series.expansion_parameter + "^" + std::to_string(i));
                    }
                }

                template<typename T>
                void name_sum(std::vector<T>& vector, std::string prefix){
                    for(int i = 0; i < vector.size(); i++){
                        name_sum(vector[i], prefix);
                    }
                }

            public:

                /*
                 * member fields
                 */
                bool verbose;
                real_t min_decrease_factor;
                real_t soft_wall_clock_limit;
                real_t hard_wall_clock_limit;
                size_t number_of_threads;
                size_t reset_cuda_after;
                container_t<sum_t> expression;
                real_t epsrel;
                real_t epsabs;
                unsigned long long int maxeval;
                unsigned long long int mineval;
                real_t maxincreasefac;
                real_t min_epsrel;
                real_t min_epsabs;
                real_t max_epsrel;
                real_t max_epsabs;

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
                    number_of_threads(0),
                    reset_cuda_after(0),
                    expression( deep_apply(expression,convert_to_sum_t) ),
                    epsrel(epsrel),epsabs(epsabs),maxeval(maxeval),mineval(mineval),maxincreasefac(maxincreasefac),min_epsrel(min_epsrel),
                    min_epsabs(min_epsabs),max_epsrel(max_epsrel),max_epsabs(max_epsabs),changed_deformation_parameters_map(read_map_from_file())
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
                    name_sum(this->expression, "sum");
                };

            private:

                decltype(std::chrono::steady_clock::now()) start_time;

            protected:

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

                void print_result(){
                    std::cout << "Current result:" << std::endl;
                    auto result = evaluate_expression();
                    for (unsigned int amp_idx = 0; amp_idx < result.size(); ++amp_idx)
                        std::cout << "amplitude" << amp_idx << " = " << result.at(amp_idx) << std::endl;
                }

                /*
                 * estimate total time and try to get below wall_clock_limit
                 */
                void ensure_wall_clock_limit(bool& repeat, std::vector<integral_t*>& integrals)
                {
                    real_t elapsed_time = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - start_time).count();
                    real_t remaining_time = hard_wall_clock_limit - elapsed_time;

                    if(remaining_time <= 0.)
                    {
                        if(verbose)
                        {
                            std::cout << std::endl;
                            std::cout << "elapsed time: " << elapsed_time << " s = " << elapsed_time/60 << " min = " << elapsed_time/60/60 << " hr" << std::endl;
                            std::cout << "remaining time: " << remaining_time << " s = " << remaining_time/60 << " min = " << remaining_time/60/60 << " hr" << std::endl;
                            std::cout << "stopping due to time constraint" << std::endl;
                        }
                        return;
                    }

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
                        std::cout << "elapsed time: " << elapsed_time << " s = " << elapsed_time/60 << " min = " << elapsed_time/60/60 << " hr" << std::endl;
                        std::cout << "estimated time for integrations: " << time_for_next_iteration << " s = " << time_for_next_iteration/60 << " min = " << time_for_next_iteration/60/60 << " hr" << std::endl;
                        std::cout << "remaining time: " << remaining_time << " s = " << remaining_time/60 << " min = " << remaining_time/60/60 << " hr" << std::endl;
                        std::cout << "soft wall clock limit: " << soft_wall_clock_limit << " s = " << soft_wall_clock_limit/60 << " min = " << soft_wall_clock_limit/60/60 << " hr" << std::endl;
                        std::cout << "hard wall clock limit: " << hard_wall_clock_limit << " s = " << hard_wall_clock_limit/60 << " min = " << hard_wall_clock_limit/60/60 << " hr" << std::endl;
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
                            std::cout << "estimated time for integrations (after reduction of samples): " << time_for_next_iteration << " s = " <<
                                    time_for_next_iteration/60 << " min = " << time_for_next_iteration/60/60 << " hr" << std::endl;
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
                                    std::cout << "sum: " << sum.display_name << ", term: " << term.display_name << ", integral " << term.integral->id << ": " <<
                                                term.integral->display_name << ", current integral result: " << result << std::endl;
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
                            std::cout << std::endl << "sum: " << sum.display_name << ", current sum result: " << current_sum << std::endl;

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
                                std::cout << "sum: " << sum.display_name << ", term: " << term.display_name << ", integral " << term.integral->id << ": " <<
                                     term.integral->display_name << ", contribution to sum error: " << abserr << std::endl;

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

                // ensure at least one thread
                if(number_of_threads == 0)
                    ++number_of_threads;

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

                // read in changed deformation parameters from file
                for(int i = 0; i < integrals.size(); i++){
                    integrals[i]->id = i+1;
                    auto element = changed_deformation_parameters_map.find(integrals[i]->display_name);
                    if(element != changed_deformation_parameters_map.end()){
                        auto new_parameters = element->second;
                        auto old_parameters = integrals[i]->get_parameters();
                        for(int k = 0; k < new_parameters.size(); k++){
                            for(int i = 0; i < new_parameters[k].size(); i++){
                                *old_parameters[k][i] = new_parameters[k][i];
                                if(verbose)
                                    std::cout << "read in changed parameter " << k << "," << i << " for " << integrals[i]->display_name << ": " << *old_parameters[k][i] << std::endl;
                            }
                        }
                    }
                }

                // initialize with minimal number of sampling points
                secdecutil::deep_apply(expression, ensure_mineval);
                if(verbose){
                    print_datetime("Starting calculations: ");
                    std::cout << "computing integrals to satisfy mineval " << this->mineval << std::endl;
                }
                evaluate_integrals(integrals, verbose, number_of_threads, reset_cuda_after, changed_deformation_parameters_map);
                if(verbose){
                    std::cout << "---------------------" << std::endl << std::endl;
                    auto elapsed_time = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - start_time).count();
                    std::cout << "elapsed time: " << elapsed_time << " s = " << elapsed_time/60 << " min = " << elapsed_time/60/60 << " hr" << std::endl;
                    print_datetime();
                    print_result();
                    std::cout << std::endl;
                }

                // ensure each integral is at least known to min_epsrel and min_epsabs
                do {
                    repeat = false;

                    secdecutil::deep_apply(expression, ensure_min_prec);
                    secdecutil::deep_apply(expression, ensure_maxincreasefac);

                    if(verbose)
                        std::cout << std::endl << "computing integrals to satisfy min_epsrel " << this->min_epsrel << " or min_epsabs " << this->min_epsabs << std::endl;
                    evaluate_integrals(integrals, verbose, number_of_threads, reset_cuda_after, changed_deformation_parameters_map);
                    if(verbose){
                        std::cout << "---------------------" << std::endl << std::endl;
                        auto elapsed_time = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - start_time).count();
                        std::cout << "elapsed time: " << elapsed_time << " s = " << elapsed_time/60 << " min = " << elapsed_time/60/60 << " hr" << std::endl;
                        print_datetime();
                        print_result();
                        std::cout << std::endl;
                    }
                } while(repeat);

                // ensure the error goal of each sum is reached as good as possible under the given time constraint
                do {
                    repeat = false;

                    secdecutil::deep_apply(expression, ensure_error_goal);
                    secdecutil::deep_apply(expression, ensure_maxincreasefac);
                    ensure_wall_clock_limit(repeat, integrals);

                    if(verbose)
                        std::cout << std::endl << "computing integrals to satisfy error goals on sums: epsrel " << this->epsrel << ", epsabs " << this->epsabs << std::endl;
                    evaluate_integrals(integrals, verbose, number_of_threads, reset_cuda_after, changed_deformation_parameters_map);
                    if(verbose){
                        std::cout << "---------------------" << std::endl << std::endl;
                        print_datetime();
                        print_result();
                        std::cout << std::endl;
                    }
                } while(repeat);
                if(verbose){
                        auto elapsed_time = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - start_time).count();
                        std::cout << "elapsed time: " << elapsed_time << " s = " << elapsed_time/60 << " min = " << elapsed_time/60/60 << " hr" << std::endl;
                        print_datetime();
                }

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
#ifndef SecDecUtil_weighted_integral_hpp_included
#define SecDecUtil_weighted_integral_hpp_included

#include <algorithm> // std::max, std::min, std::sort
#include <cassert> // assert
#include <chrono> // std::chrono::steady_clock
#include <iostream> // std::cerr, std::dec
#include <iomanip> // std::fixed, std::setprecision
#include <limits> // std::numeric_limits
#include <memory> // std::shared_ptr
#include <string> // std::to_string
#include <stdexcept> // std::domain_error, std::logic_error, std::runtime_error
#include <thread> // std::thread
#include <utility> // std::declval, std::move
#include <vector> // std::vector
#include <fstream> // for writing changed deformation parameters to file
#include <sstream> // std::ostringstream 
#include <cstdio> // std::remove 
#include <unistd.h> // mkdtemp
#include <cstring> // strcpy

#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply
#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation
#include <secdecutil/sector_container.hpp> // for secdecutil::sign_check_error

// wrapped integrators
#include <secdecutil/integrators/cuba.hpp>
#include <secdecutil/integrators/qmc.hpp>
// cannot set the number of points for CQuad --> not suitable for "Integral"

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
                bool allow_refine = true; // indicates whether the integral can be improved with further iterations

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
                void set_next_number_of_function_evaluations(unsigned long long int new_number_of_function_evaluations, bool allow_decrease=false)
                {
                    if(allow_decrease)
                        next_number_of_function_evaluations = std::max({new_number_of_function_evaluations, number_of_function_evaluations});
                    else
                        next_number_of_function_evaluations = std::max({next_number_of_function_evaluations, new_number_of_function_evaluations, number_of_function_evaluations});
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
                virtual void compute_impl(const bool verbose) = 0; // implementation should populate "integral_result"
                void compute(const bool verbose)
                {
                    if(allow_refine && (next_number_of_function_evaluations > number_of_function_evaluations)) {
                        auto start_time = std::chrono::steady_clock::now();
                        compute_impl(verbose);
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

            void compute_impl(const bool verbose) override
            {
                using std::abs;

                unsigned long long int desired_next_n = this->get_next_number_of_function_evaluations();
                unsigned long long int next_n = qmc->get_next_n(desired_next_n);
                if(next_n < desired_next_n)
                {
                    this->allow_refine = false; //don't iterate with same/smaller lattice
                    // warn user that they require too many function evaluations, we return the result from the largest lattice
                    if(verbose)
                        std::cerr << "WARNING class QmcIntegral: The requested number_of_function_evaluations ("
                            + std::to_string(desired_next_n) + ") exceeds the largest available lattice ("
                            + std::to_string(next_n) +"), using largest available lattice." << std::endl;
                }
                
                // set number of function evaluations to the next larger lattice (allow decrease, e.g. if next_number_of_function_evaluations > largest lattice)
                this->set_next_number_of_function_evaluations( next_n, true );

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
            std::string tmpdir;


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
                scaleexpo(0.5),
                tmpdir("")
                {
                  // set statefile path`
                  std::string path;
                  if(const char* env_p = std::getenv("SECDEC_TMPDIR"))
                      path+= env_p;
                  else if(const char* env_p = std::getenv("TMPDIR"))
                      path+= env_p;
                  else
                      path += "/tmp";
                  if (not path.empty())
                  {
                      path  += "/pySecDec-XXXXXX" ;
                      char buf [path.length()+1];
                      strcpy(buf, path.c_str());
                      if(mkdtemp(buf))
                        tmpdir=buf;
                  }
                };

            void compute_impl(const bool verbose) override
            {
                integrator->statefiledir = tmpdir;

                if(integrator->integrator_type == 3)
                    this->allow_refine = false; //don't iterate with Divonne
                
                if(integrator->integrator_type != 3) // Divonne
                {
                  unsigned long long int next_n = this->get_next_number_of_function_evaluations();
                  
                  integrator->mineval = integrator->maxeval = next_n; // set number of sampling points


                  integrator->epsrel = 0.;
                  integrator->epsabs = 0.;
                }
                this->integral_result = integrator->integrate(this->integrand); // run the numerical integration
                this->set_next_number_of_function_evaluations( *integrator->neval, true ); // set number of function evluations to the next larger lattice
            }

            ~CubaIntegral()
            {
                for(char i : {'0','1','2'})
                  remove((tmpdir + '/' + i).c_str()); 
                rmdir(tmpdir.c_str());
            }
        };
  
        template<typename integrand_return_t, typename real_t, typename integrator_t, typename integrand_t>
        struct CQuadIntegral : public Integral<integrand_return_t,real_t>
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
            CQuadIntegral(const std::shared_ptr<integrator_t>& integrator, const integrand_t& integrand) :
                Integral<integrand_return_t,real_t>(),
                integrator(integrator),
                integrand(integrand),
                scaleexpo(0.5)
                {};

            void compute_impl(const bool verbose) override
            {
                this->allow_refine = false; //don't iterate with CQuad

                unsigned long long int next_n = this->get_next_number_of_function_evaluations();
                
                this->integral_result = integrator->integrate(this->integrand); // run the numerical integration
            }
        };
    
        template<typename integrand_return_t, typename real_t, typename integrator_t, typename integrand_t>
        struct MultiIntegratorIntegral : public Integral<integrand_return_t,real_t>
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
            MultiIntegratorIntegral(const std::shared_ptr<integrator_t>& integrator, const integrand_t& integrand) :
                Integral<integrand_return_t,real_t>(),
                integrator(integrator),
                integrand(integrand),
                scaleexpo(0.5)
                {};

            void compute_impl(const bool verbose) override
            {
                this->allow_refine = false; //don't iterate with MultiIntegrator
                
                unsigned long long int next_n = this->get_next_number_of_function_evaluations();

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
        {   // for integrals appearing in both arguments, sum the coefficients; add new integrals from b to a
            for (const auto & y: b)
            {
                bool found=false;
                for (auto & x: a)
                {
                    if(x.integral==y.integral)
                    {
                        found=true;
                        x.coefficient += y.coefficient;
                        break;
                    }
                }
                if(not found)
                    a.push_back(y);
            }
            return a;
        };
        template<typename integral_t, typename coefficient_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>> operator+
        (
            std::vector<WeightedIntegral<integral_t,coefficient_t>> a,
            const std::vector<WeightedIntegral<integral_t,coefficient_t>>& b
        )
        {
            a+=b;
            return a;
        };
        template<typename integral_t, typename coefficient_t, typename factor_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>>& operator*=
        (
            std::vector<WeightedIntegral<integral_t,coefficient_t>>& a,
            const factor_t& b
        )
        {
            for(WeightedIntegral<integral_t,coefficient_t>& i : a)
                i.coefficient *= b;
            return a;
        };
        template<typename integral_t, typename coefficient_t, typename factor_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>> operator*
        (
            std::vector<WeightedIntegral<integral_t,coefficient_t>> a,
            const factor_t& b
        )
        {
            for(WeightedIntegral<integral_t,coefficient_t>& i : a)
                i.coefficient *= b;
            return a;
        };
        template<typename integral_t, typename coefficient_t, typename factor_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>> operator*
        (
            const factor_t& a,
            std::vector<WeightedIntegral<integral_t,coefficient_t>> b
        )
        {
            for(WeightedIntegral<integral_t,coefficient_t>& i : b)
                i.coefficient *= a;
            return b;
        };
        template<typename integral_t, typename coefficient_t, typename divisor_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>>& operator/=
        (
            std::vector<WeightedIntegral<integral_t,coefficient_t>>& a,
            const divisor_t& b
        )
        {
            for(WeightedIntegral<integral_t,coefficient_t>& i : a)
                i.coefficient /= b;
            return a;
        };
        template<typename integral_t, typename coefficient_t, typename divisor_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>> operator/
        (
            std::vector<WeightedIntegral<integral_t,coefficient_t>> a,
            const divisor_t& b
        )
        {
            for(WeightedIntegral<integral_t,coefficient_t>& i : a)
                i.coefficient /=b;
            return a;
        };
        template<typename integral_t, typename coefficient_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>> operator-
        (
            std::vector<WeightedIntegral<integral_t,coefficient_t>> a
        )
        {
            a*=-1;
            return a;
        };
        template<typename integral_t, typename coefficient_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>> operator-
        (
            std::vector<WeightedIntegral<integral_t,coefficient_t>> a,
            const std::vector<WeightedIntegral<integral_t,coefficient_t>>& b
        )
        {
            a+=-b;
            return a;
        };
        template<typename integral_t, typename coefficient_t>
        std::vector<WeightedIntegral<integral_t,coefficient_t>> operator-=
        (
            std::vector<WeightedIntegral<integral_t,coefficient_t>> &a,
            const std::vector<WeightedIntegral<integral_t,coefficient_t>>& b
        )
        {
            a+=-b;
            return a;
        };

        static inline void write_map_to_file(std::map<std::string, std::vector<std::vector<double>>> map, std::string filename="changed_deformation_parameters.txt"){
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

        static inline std::map<std::string, std::vector<std::vector<double>>> read_map_from_file(std::string filename="changed_deformation_parameters.txt"){
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

        static inline void print_datetime(std::string prefix = "Current time: "){
            auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cerr << prefix << std::ctime(&t);
        }

        /*
         * evaluate a vector of integrals
         */
        template<typename integrand_return_t, typename integral_t>
        void evaluate_integrals(std::vector<integral_t*>& integrals, const bool& verbose, size_t number_of_threads, size_t reset_cuda_after,
                    std::map<std::string, std::vector<std::vector<double>>> changed_deformation_parameters_map)
        {
            if(number_of_threads == 0)
                ++number_of_threads;

            std::function<void(integral_t*)> compute_integral = [ &verbose, &integrals, &changed_deformation_parameters_map ] (integral_t* integral)
                {
                    const unsigned long long int curr_n = integral->get_number_of_function_evaluations();
                    const unsigned long long int next_n = integral->get_next_number_of_function_evaluations();
                    secdecutil::UncorrelatedDeviation<integrand_return_t> old_result;
                    if(verbose)
                        try {
                            old_result = integral->get_integral_result();
                        } catch (const integral_not_computed_error&) { /* ignore */ };

                    bool failed_atleast_once = false;
                    while(true){
                        bool failed = false;
                        try{
                            integral->compute(verbose);
                        } catch(secdecutil::sign_check_error& e){
                            failed = true;
                            std::cerr << "Exception: " << e.what() << std::endl;
                            integral->clear_errors();

                            std::cerr << "Integral " << integral->display_name << " failed, reducing deformation parameters." << std::endl;
                            std::vector<std::vector<typename std::remove_pointer<decltype(integral)>::type::real_t_type*>> pars = integral->get_parameters();
                            std::vector<std::vector<typename std::remove_pointer<decltype(integral)>::type::real_t_type>> extra_pars = integral->get_extra_parameters();
                            failed_atleast_once = true;
                            bool changed_deformation_parameters = false;
                            for(int k = 0; k < pars.size(); k++){
                                for(int i = 0; i < pars[k].size(); i++){
                                    std::cerr << "par " << k << "," << i << ": " << *pars[k][i];
                                    if(*pars[k][i] > extra_pars[k][0]){
                                        *pars[k][i] *= extra_pars[k][1];
                                        if(*pars[k][i] < extra_pars[k][0]){
                                            *pars[k][i] = extra_pars[k][0];
                                        }
                                        changed_deformation_parameters = true;
                                        std::cerr << " -> " << *pars[k][i];
                                    }
                                    if(i == pars[k].size()-1 and k == pars.size()-1)
                                        std::cerr << std::endl;
                                    else
                                        std::cerr << ", ";
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
                                //write_map_to_file(changed_deformation_parameters_map); // do not store changed lambda parameters on disk
                            }
                            break;
                        }
                    }

                    if(verbose and (next_n > curr_n))
                    {
                        std::cerr << "integral " << integral->id << "/" << integrals.size() << ": " << integral->display_name << ", time: ";
                        auto flags = std::cerr.flags();
                        std::cerr << std::fixed << std::setprecision(6) << integral->get_integration_time() << "s, ";
                        std::cerr.flags(flags);
                        print_datetime();
                        std::cerr << "res: " << old_result << " -> " << integral->get_integral_result()
                                      << ", n: " << curr_n << " -> " << std::dec << integral->get_number_of_function_evaluations() << std::endl;
                        std::cerr << std::endl;
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
            #ifdef SECDEC_WITH_CUDA
                typedef thrust::complex<real_t> complex_t;
            #else
                typedef std::complex<real_t> complex_t;
            #endif
            public:

                using integral_t = Integral<integrand_return_t,real_t>;
                using term_t = WeightedIntegral<integral_t,coefficient_t>;
                using sum_base_t = decltype(std::declval<coefficient_t>()*std::declval<integrand_return_t>());
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
                real_t decrease_to_percentage; // of remaining time
                real_t wall_clock_limit;
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

                enum ErrorMode : int { abs=0, all, largest, real, imag};
                ErrorMode errormode;

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
                    real_t max_epsabs = 1e-20,
                    ErrorMode errormode=ErrorMode::abs
                ) :
                    verbose(false),
                    min_decrease_factor(0.9),
                    wall_clock_limit(std::numeric_limits<double>::infinity()),
                    decrease_to_percentage(0.7),
                    number_of_threads(0),
                    reset_cuda_after(0),
                    expression( deep_apply(expression,convert_to_sum_t) ),
                    epsrel(epsrel),epsabs(epsabs),maxeval(maxeval),mineval(mineval),maxincreasefac(maxincreasefac),min_epsrel(min_epsrel),
                    min_epsabs(min_epsabs),max_epsrel(max_epsrel),max_epsabs(max_epsabs),changed_deformation_parameters_map(read_map_from_file()),
                    errormode(errormode)
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
                    std::cerr << "Current result:" << std::endl;
                    container_t<sum_return_t> result = evaluate_expression();
                    for (unsigned int amp_idx = 0; amp_idx < result.size(); ++amp_idx)
                        std::cerr << "amplitude" << amp_idx << " = " << result.at(amp_idx) << std::endl;
                }

                real_t apply_errormode(sum_base_t error)
                {
                    using std::abs;

                    if (errormode == this->abs)
                        return abs(error);
                    if (errormode == real)
                        return abs(complex_t(error).real());
                    if (errormode == imag)
                        return abs(complex_t(error).imag());
                    throw std::runtime_error("unexpected errormode");
                }

                /*
                 * estimate total time and try to get below wall_clock_limit
                 */
                void ensure_wall_clock_limit(bool& repeat, std::vector<integral_t*>& integrals, double decrease_to_percentage=1.0)
                {
                    real_t elapsed_time = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - start_time).count();
                    real_t remaining_time = wall_clock_limit - elapsed_time;

                    if(remaining_time <= 0.)
                    {
                        repeat = false;
                        if(verbose)
                        {
                            std::cerr << std::endl;
                            std::cerr << "elapsed time: " << elapsed_time << " s = " << elapsed_time/60 << " min = " << elapsed_time/60/60 << " hr" << std::endl;
                            std::cerr << "remaining time: " << remaining_time << " s = " << remaining_time/60 << " min = " << remaining_time/60/60 << " hr" << std::endl;
                            std::cerr << "stopping due to time constraint" << std::endl;
                        }
                        for (integral_t* integral : integrals)
                        {
                            integral->set_next_number_of_function_evaluations(  integral->get_number_of_function_evaluations(), true );
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
                        std::cerr << std::endl;
                        std::cerr << "elapsed time: " << elapsed_time << " s = " << elapsed_time/60 << " min = " << elapsed_time/60/60 << " hr" << std::endl;
                        std::cerr << "estimated time for integrations: " << time_for_next_iteration << " s = " << time_for_next_iteration/60 << " min = " << time_for_next_iteration/60/60 << " hr" << std::endl;
                        std::cerr << "remaining time: " << remaining_time << " s = " << remaining_time/60 << " min = " << remaining_time/60/60 << " hr" << std::endl;
                        std::cerr << "wall clock limit: " << wall_clock_limit << " s = " << wall_clock_limit/60 << " min = " << wall_clock_limit/60/60 << " hr" << std::endl;
                    }

                    // decrease number of sampling points if sampling would run out of time
                    while(remaining_time * decrease_to_percentage < time_for_next_iteration)
                    {
                        if(verbose)
                            std::cerr << "reducing number of samples due to time constraint ..." << std::endl;

                        decrease_factor = remaining_time * decrease_to_percentage / time_for_next_iteration ;
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
                                integral->set_next_number_of_function_evaluations( reduced_next_n, true );

                                if(reduced_next_n > curr_n) {
                                    can_improve_in_time = true;
                                    time_for_next_iteration += integral->get_integration_time() * n_increase_fac;
                                }
                            }
                        }

                        // stop iterating if wall-clock limit passed
                        if(elapsed_time > wall_clock_limit)
                        {
                            if(verbose)
                                std::cerr << "elapsed_time > wall_clock_limit: " << elapsed_time << " > " << wall_clock_limit << ", no further refinements (ensure_wall_clock_limit)" << std::endl;
                            repeat = false;
                        }
                        else
                            // keep iterating only if improvements are possible within the time limit
                            repeat = can_improve_in_time;

                        if(verbose)
                        {
                            std::cerr << "estimated time for integrations (after reduction of samples): " << time_for_next_iteration << " s = " <<
                                    time_for_next_iteration/60 << " min = " << time_for_next_iteration/60/60 << " hr" << std::endl;
                            std::cerr << "can improve in time: " << (can_improve_in_time ? "true" : "false") << std::endl;
                        }
                    }

                    if(verbose)
                        std::cerr << "run further refinements: " << (repeat ? "true" : "false") << std::endl;

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
                    [ ] (sum_t& sum)
                    {
                        for (term_t& term : sum.summands)
                            term.integral->set_next_number_of_function_evaluations( sum.mineval );
                    };

                std::function<void(sum_t&)> ensure_min_prec =
                    [ this, &repeat ] (sum_t& sum)
                    {
                        for (term_t& term : sum.summands)
                        {
                            if(term.integral->allow_refine && term.integral->get_number_of_function_evaluations() < sum.maxeval)
                            {
                                secdecutil::UncorrelatedDeviation<integrand_return_t> result = term.integral->get_integral_result();
                                real_t abserr = abs(result.uncertainty);
                                if(abserr > sum.min_epsabs)
                                {
                                    real_t relerr = abserr / abs( result.value );
                                    if(relerr > sum.min_epsrel)
                                    {
                                        unsigned long long int old_n = term.integral->get_number_of_function_evaluations();
                                        unsigned long long int new_n = static_cast<unsigned long long int>(
                                                                        min( old_n * pow(min(relerr/sum.min_epsrel,abserr/sum.min_epsabs), 1./term.integral->get_scaleexpo()),
                                                                        static_cast<real_t>(std::numeric_limits<long long>::max()))  // in some extreme cases (e.g. test cases using MC with tiny error goal), next_n could exceed LLONG_MAX
                                                                       );
                                        
                                        if(new_n > old_n)
                                        {
                                            repeat = true;
                                            term.integral->set_next_number_of_function_evaluations( new_n > sum.maxeval ? sum.maxeval : new_n );
                                            if(verbose)
                                                std::cerr << "sum: " << sum.display_name << ", term: " << term.display_name << ", integral " << term.integral->id << ": " <<
                                                        term.integral->display_name << ", current integral result: " << result << ", increase n: " << old_n << " -> " <<new_n << std::endl;
                                        } else {
                                            if(verbose)
                                                std::cerr << "sum: " << sum.display_name << ", term: " << term.display_name << ", integral " << term.integral->id << ": " <<
                                                    term.integral->display_name << ", new_n <= old_n: " << new_n << " < " << old_n << " (ensure_min_prec)" << std::endl;

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
                        if(sum.epsrel == 0. and sum.epsabs == 0.)
                            sum.epsrel=1e-50;

                        ErrorMode original_errormode = errormode; // for errormode==largest, we will temporarily set errormode = real or imag, depending on which error is larger
                        // compute absolute error goal
                        sum_return_t current_sum = compute_sum(sum);
                        sum_base_t current_sum_value = current_sum.value;
                        sum_base_t current_sum_error = current_sum.uncertainty;
                        if (errormode == largest)
                        {
                            errormode = (complex_t(current_sum.uncertainty).real() > complex_t(current_sum.uncertainty).imag()) ? real : imag;
                        }

                        real_t abs_error = apply_errormode(current_sum_error);
                        real_t abs_error_goal = sum.epsrel * apply_errormode(current_sum_value);
                        if(original_errormode == largest)
                            abs_error_goal = sum.epsrel * std::max(abs(complex_t(current_sum_value).real()), abs(complex_t(current_sum_value).imag()));
                        abs_error_goal = max(abs_error_goal, abs_error*sum.epsrel); // If current error larger than current result set goal based on error
                        abs_error_goal = max(abs_error_goal, sum.epsabs); // Do not request an error smaller than epsabs

                        if(verbose)
                            std::cerr << std::endl << "sum: " << sum.display_name << ", current sum result: " << current_sum << ",  error goal: " << abs_error_goal << " (ensure_error_goal)" << std::endl;

                        if(abs_error < abs_error_goal)
                            {
                                if(verbose)
                                    std::cerr << "sum: " << sum.display_name << ", current sum result: " << current_sum << ", abs_error < abs_error_goal: " << abs_error << " < " << abs_error_goal << ", no further refinements (ensure_error_goal)" << std::endl;
                                errormode = original_errormode;
                                return;
                            }
                        abs_error_goal *= abs_error_goal; // require variance rather than standard deviation in the following


                        // compute the C_i
                        std::vector<real_t> c; c.reserve( sum.summands.size() );
                        for(auto& term : sum.summands)
                        {
                            auto time = term.integral->get_integration_time();
                            secdecutil::UncorrelatedDeviation<integrand_return_t> integral = term.integral->get_integral_result();
                            real_t abserr = apply_errormode(integral.uncertainty * term.coefficient);
                            real_t relerr = abserr / abs( integral.value * term.coefficient);

                            if(!term.integral->allow_refine || relerr < sum.max_epsrel || abserr < sum.max_epsabs)
                                c.push_back(  -1  );
                            else
                                c.push_back(  sqrt(time) / abserr  );

                            if(verbose)
                            {
                                if(!term.integral->allow_refine)
                                    std::cerr << "sum: " << sum.display_name << ", term: " << term.display_name << ", integral " << term.integral->id << ": " <<
                                                term.integral->display_name << ", allow_refine = false, no further refinements (ensure_error_goal)" << std::endl;

                                if(relerr < sum.max_epsrel)
                                    std::cerr << "sum: " << sum.display_name << ", term: " << term.display_name << ", integral " << term.integral->id << ": " <<
                                                term.integral->display_name << ", relerr < max_epsrel: " << relerr << " < " << max_epsrel << ", no further refinements (ensure_error_goal)" << std::endl;

                                if(abserr < sum.max_epsabs)
                                    std::cerr << "sum: " << sum.display_name << ", term: " << term.display_name << ", integral " << term.integral->id << ": " <<
                                                term.integral->display_name << ", abserr < max_epsabs: " << abserr << " < " << max_epsabs << ", no further refinements (ensure_error_goal)" << std::endl;
                            }
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
                                secdecutil::UncorrelatedDeviation<integrand_return_t> result = term.integral->get_integral_result();
                                real_t abserr = apply_errormode(result.uncertainty * term.coefficient); // contribution to total error of the sum
                                real_t absvar = abserr * abserr;
                                if(c_i > 0)
                                {
                                    const real_t scaleexpo = term.integral->get_scaleexpo();
                                    const real_t beta = 2./(2.*scaleexpo+1.);
                                    const real_t beta_times_scaleexpo = beta*scaleexpo;

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

                        // reduce number of sampling points if estimated time for next iteration > decrease_to_percentage*time_remaining
                        real_t elapsed_time = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - start_time).count();
                        real_t time_remaining = (wall_clock_limit - elapsed_time)*decrease_to_percentage;
                        do {
                            oldfac = fac;

                            real_t time_estimate(0);
                            real_t scaleexpo;

                            for(size_t i = 0; i < sum.summands.size(); ++i)
                            {
                                real_t& c_i = c.at(i);
                                term_t& term = sum.summands.at(i);
                                if(c_i > 0)
                                {
                                    scaleexpo = term.integral->get_scaleexpo();
                                    const real_t beta = 2./(2.*scaleexpo+1.);

                                    real_t temp = pow(c_i,-beta)*pow(fac,0.5/scaleexpo);
                                    if (temp>0)
                                        time_estimate+= term.integral->get_integration_time()*temp;
                                }
                            }

                            if (time_estimate > time_remaining)
                            {
                                fac *= pow(time_remaining/time_estimate, 2*scaleexpo); 
                                if(verbose)
                                    std::cerr << "using reduced number of sampling points due to time limit" << std::endl;
                            }
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
                                const real_t scaleexpo = term.integral->get_scaleexpo();
                                const real_t beta = 2./(2.*scaleexpo+1.);

                                unsigned long long int proposed_next_n =
                                    max(
                                            next_n,
                                            static_cast<unsigned long long int>( min( curr_n*pow(c_i,-beta)*pow(fac,0.5/scaleexpo), 
                                                                                 static_cast<real_t>(std::numeric_limits<long long>::max()))  // in some extreme cases (e.g. test cases using MC with tiny error goal), proposed_next_n could exceed LLONG_MAX
                                                                               )
                                       );
                                proposed_next_n = min(sum.maxeval, proposed_next_n);

                                term.integral->set_next_number_of_function_evaluations(proposed_next_n);

                                // run again if the accuracy of an integral is increased
                                if(proposed_next_n > curr_n)
                                {
                                    repeat = true;
                                    if(verbose)
                                        std::cerr << "sum: " << sum.display_name << ", term: " << term.display_name << ", integral " << term.integral->id << ": " <<
                                                term.integral->display_name << ", current integral result: " << term.integral->get_integral_result() <<
                                                ", contribution to sum error: " << apply_errormode(term.integral->get_integral_result().uncertainty*term.coefficient) <<
                                                ", increase n: " << curr_n << " -> " <<proposed_next_n << " (ensure_error_goal)" << std::endl;
                                } else {
                                    if(verbose)
                                        std::cerr << "sum: " << sum.display_name << ", term: " << term.display_name << ", integral " << term.integral->id << ": " <<
                                                term.integral->display_name << ", proposed_next_n <= curr_n: " << proposed_next_n << " < " << curr_n << " (ensure_error_goal)" << std::endl;
                                }
                            }

                        }
                        errormode = original_errormode;
                    };

                // Damp very large increase in number of points
                std::function<void(sum_t&)> ensure_maxincreasefac =
                    [ this ] (sum_t& sum)
                    {
                        for(auto& term: sum.summands)
                        {
                            // Do not increase points by more than a factor of maxincreasefac
                            const real_t curr_n_as_real = term.integral->get_number_of_function_evaluations(); // implicit cast to real_t
                            const unsigned long long int next_n = term.integral->get_next_number_of_function_evaluations();
                            const unsigned long long int max_next_n = curr_n_as_real * sum.maxincreasefac; // implicit cast to unsigned long long int

                            if( next_n > max_next_n )
                            {
                                term.integral->set_next_number_of_function_evaluations( max_next_n, true );
                                if(verbose)
                                    std::cerr << "sum: " << sum.display_name << ", term: " << term.display_name << ", integral " << term.integral->id << ": " <<
                                            term.integral->display_name << ", current integral result: " << term.integral->get_integral_result() <<
                                            ", contribution to sum error: " << abs(term.integral->get_integral_result().uncertainty*term.coefficient) <<
                                            ", decreasing next_n: " << next_n << " -> max_next_n: " << max_next_n << " (ensure_maxincreasefac)" << std::endl;
                            }
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
                //    auto element = changed_deformation_parameters_map.find(integrals[i]->display_name);
                //    if(element != changed_deformation_parameters_map.end()){
                //        auto new_parameters = element->second;
                //        auto old_parameters = integrals[i]->get_parameters();
                //        for(int k = 0; k < new_parameters.size(); k++){
                //            for(int j = 0; j < new_parameters[k].size(); j++){
                //                *old_parameters[k][j] = new_parameters[k][j];
                //                if(verbose)
                //                    std::cerr << "read in changed parameter " << k << "," << j << " for " << integrals[i]->display_name << ": " << *old_parameters[k][j] << std::endl;
                //            }
                //        }
                //    }
                }

                // initialize with minimal number of sampling points
                secdecutil::deep_apply(expression, ensure_mineval);
                if(verbose){
                    print_datetime("Starting calculations: ");
                    std::cerr << "computing integrals to satisfy mineval " << this->mineval << std::endl;
                }
                evaluate_integrals<integrand_return_t>(integrals, verbose, number_of_threads, reset_cuda_after, changed_deformation_parameters_map);
                if(verbose){
                    std::cerr << "---------------------" << std::endl << std::endl;
                    auto elapsed_time = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - start_time).count();
                    std::cerr << "elapsed time: " << elapsed_time << " s = " << elapsed_time/60 << " min = " << elapsed_time/60/60 << " hr" << std::endl;
                    print_datetime();
                    print_result();
                    std::cerr << std::endl;
                }

                // ensure each integral is at least known to min_epsrel and min_epsabs
                do {
                    repeat = false;

                    if(verbose)
                        std::cerr << std::endl << "computing integrals to satisfy min_epsrel " << this->min_epsrel << " or min_epsabs " << this->min_epsabs << std::endl;

                    secdecutil::deep_apply(expression, ensure_min_prec);
                    if(verbose)
                        std::cerr << "ensure_error_goal requires further refinements: " << (repeat ? "true" : "false") << std::endl;
                    secdecutil::deep_apply(expression, ensure_maxincreasefac);
                    if(verbose)
                        std::cerr << "ensure_maxincreasefac allows/requires further refinements: " << (repeat ? "true" : "false") << std::endl;
                    ensure_wall_clock_limit(repeat, integrals, decrease_to_percentage);
                    if(verbose)
                        std::cerr << "ensure_wall_clock_limit allows/requires further refinements: " << (repeat ? "true" : "false") << std::endl;

                    evaluate_integrals<integrand_return_t>(integrals, verbose, number_of_threads, reset_cuda_after, changed_deformation_parameters_map);
                    if(verbose){
                        std::cerr << "---------------------" << std::endl << std::endl;
                        auto elapsed_time = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - start_time).count();
                        std::cerr << "elapsed time: " << elapsed_time << " s = " << elapsed_time/60 << " min = " << elapsed_time/60/60 << " hr" << std::endl;
                        print_datetime();
                        print_result();
                        std::cerr << std::endl;
                    }
                } while(repeat);

                // ensure the error goal of each sum is reached as good as possible under the given time constraint
                do {
                    repeat = false;
                    if(verbose)
                        std::cerr << std::endl << "computing integrals to satisfy error goals on sums: epsrel " << this->epsrel << ", epsabs " << this->epsabs << std::endl;

                    if (errormode == all)
                    {
                        if(verbose) std::cerr << std::endl << "ensure error goal for real part" << std::endl;
                        errormode = real;
                        secdecutil::deep_apply(expression, ensure_error_goal);
                        if(verbose) std::cerr << std::endl << "ensure error goal for imag part" << std::endl;
                        errormode = imag;
                        secdecutil::deep_apply(expression, ensure_error_goal);
                        errormode = all;
                    }
                    else 
                    {
                        secdecutil::deep_apply(expression, ensure_error_goal);
                    }
                    if(verbose)
                        std::cerr << "ensure_error_goal requires further refinements: " << (repeat ? "true" : "false") << std::endl;
                    secdecutil::deep_apply(expression, ensure_maxincreasefac);
                    if(verbose)
                        std::cerr << "ensure_maxincreasefac allows/requires further refinements: " << (repeat ? "true" : "false") << std::endl;
                    // ensure_error_goal already implements wall_clock_limit, but separately for each sum. Here we check that when considering all sums, the wall_clock_limit is not exceeded.
                    ensure_wall_clock_limit(repeat, integrals, 1.0); 
                    if(verbose)
                        std::cerr << "ensure_wall_clock_limit allows/requires further refinements: " << (repeat ? "true" : "false") << std::endl;

                    evaluate_integrals<integrand_return_t>(integrals, verbose, number_of_threads, reset_cuda_after, changed_deformation_parameters_map);
                    if(verbose){
                        std::cerr << "---------------------" << std::endl << std::endl;
                        print_datetime();
                        print_result();
                        std::cerr << std::endl;
                    }
                } while(repeat);
                if(verbose){
                        auto elapsed_time = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - start_time).count();
                        std::cerr << "elapsed time: " << elapsed_time << " s = " << elapsed_time/60 << " min = " << elapsed_time/60/60 << " hr" << std::endl;
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

namespace std {

    template <class T1, class T2, class T3>
    struct common_type<std::vector<secdecutil::amplitude::WeightedIntegral<T1,T2>>, T3> {
        using type = std::vector<secdecutil::amplitude::WeightedIntegral<T1,T2>>;
    };
        
}

#endif

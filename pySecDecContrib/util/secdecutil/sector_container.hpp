#ifndef SecDecUtil_sector_container_hpp_included
#define SecDecUtil_sector_container_hpp_included

#include <exception>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <gsl/gsl_qrng.h>
#include <secdecutil/integrand_container.hpp>

namespace secdecutil {

    #ifdef SECDEC_WITH_CUDA
        // this error is thrown if a call to a cuda function fails
        struct cuda_error : public std::runtime_error { using std::runtime_error::runtime_error; };
    #endif

    // this error is thrown if an error with the gsl occurs
    struct gsl_error : public std::runtime_error { using std::runtime_error::runtime_error; };


    /*
     * The container types pySecDec packs the integrand functions in.
     * Depending on whether we have contour deformation or not, the
     * appropriate container type is chosen.
     */
    template<typename real_t, typename complex_t, typename integrand_return_t>
    struct SectorContainerWithoutDeformation
    {
        // the call signature of an integrand
        typedef integrand_return_t IntegrandFunction
        (
         real_t const * const integration_variables,
         real_t const * const real_parameters,
         complex_t const * const complex_parameters,
         ResultInfo* result_info
        );

        #ifdef SECDEC_WITH_CUDA
            // the call signature of the get device integrand function
            typedef IntegrandFunction* GetDeviceIntegrandFunction();
        #endif

        const unsigned sector_id;
        const std::vector<int> orders;
        const unsigned number_of_integration_variables;
        IntegrandFunction * const undeformed_integrand;
        #ifdef SECDEC_WITH_CUDA
            GetDeviceIntegrandFunction * const get_device_undeformed_integrand;
            IntegrandFunction * device_undeformed_integrand;
        #endif

        std::shared_ptr<std::vector<real_t>> real_parameters;
        std::shared_ptr<std::vector<complex_t>> complex_parameters;
        // "integrand" must be a member function, otherwise we cannot bind the struct
        integrand_return_t integrand
        (
            real_t const * const integration_variables,
            real_t const * const real_parameters,
            complex_t const * const complex_parameters,
            ResultInfo* result_info
        )
        {
            return undeformed_integrand(integration_variables, real_parameters, complex_parameters, result_info);
        };

        // constructor
        SectorContainerWithoutDeformation
        (
            const unsigned sector_id,
            const std::vector<int> orders,
            const unsigned number_of_integration_variables,
            IntegrandFunction * const undeformed_integrand
            #ifdef SECDEC_WITH_CUDA
                ,GetDeviceIntegrandFunction * const get_device_undeformed_integrand
            #endif
        ) :
        sector_id(sector_id),
        orders(orders),
        number_of_integration_variables(number_of_integration_variables),
        undeformed_integrand(undeformed_integrand)
        #ifdef SECDEC_WITH_CUDA
            ,get_device_undeformed_integrand(get_device_undeformed_integrand)
        #endif
        {};
    };

    template<typename real_t, typename complex_t>
    struct SectorContainerWithDeformation
    {
        // the call signature of a deformed integrand
        typedef complex_t DeformedIntegrandFunction
        (
         real_t const * const integration_variables,
         real_t const * const real_parameters,
         complex_t const * const complex_parameters,
         real_t const * const deformation_parameters,
         ResultInfo* result_info
        );

        #ifdef SECDEC_WITH_CUDA
            // the call signature of the get device deformed integrand function
            typedef DeformedIntegrandFunction* GetDeviceDeformedIntegrandFunction();
        #endif

        // the call signature of the function that computes the maximal deformation parameters
        typedef void MaximalDeformationFunction
        (
         real_t * output_deformation_parameters,
         real_t const * const integration_variables,
         real_t const * const real_parameters,
         complex_t const * const complex_parameters,
         ResultInfo* result_info
        );

        const unsigned sector_id;
        const std::vector<int> orders;
        const unsigned number_of_integration_variables;
        DeformedIntegrandFunction * const deformed_integrand;
        #ifdef SECDEC_WITH_CUDA
            GetDeviceDeformedIntegrandFunction * const get_device_deformed_integrand;
            DeformedIntegrandFunction * device_deformed_integrand;
        #endif
        DeformedIntegrandFunction * const contour_deformation_polynomial;
        MaximalDeformationFunction * const maximal_allowed_deformation_parameters;

        // the function that optimizes the deformation_parameters
        std::vector<real_t> optimize_deformation_parameters (
                                                                real_t const * const real_parameters,
                                                                complex_t const * const complex_parameters,
                                                                const unsigned number_of_presamples = 100000,
                                                                const real_t maximum = 1.,
                                                                const real_t minimum = 1.e-5,
                                                                const real_t decrease_factor = 0.9
                                                            ) const
        {
            // define indices for the loops
            unsigned i,j;

            // if no sampling desired (number_of_presamples == 0) set the deformation parameters to the maximum
            if (number_of_presamples == 0)
                return std::vector<real_t>(number_of_integration_variables,maximum);

            // initialize the output, and temporary vectors
            std::vector<real_t> optimized_deformation_parameter_vector(number_of_integration_variables,maximum);
            std::vector<real_t> temp_deformation_parameter_vector(number_of_integration_variables,0);
            std::vector<real_t> real_sample_vector(number_of_integration_variables,0);
            real_t * optimized_deformation_parameters = optimized_deformation_parameter_vector.data();
            real_t * temp_deformation_parameters = temp_deformation_parameter_vector.data();
            real_t * real_sample = real_sample_vector.data();
            ResultInfo result_info;

            // define a Sobol sequence using the gsl
            // Restriction to at most 40 dimensions only because of the implementation in the gsl. --> Use a different Sobol implementation if higher dimensionality is needed.
            int Sobol_maxdim = 40;
            if (number_of_integration_variables > Sobol_maxdim)
                throw gsl_error("The gsl implements Sobol sequences only up to " + std::to_string(Sobol_maxdim) +" dimensions (need " +
                                std::to_string(number_of_integration_variables) + "). Please set the \"deformation_parameters\" manually.");

            // define the generator
            gsl_qrng * Sobol_generator = gsl_qrng_alloc(gsl_qrng_sobol, number_of_integration_variables);

            // find the minimum of the lambdas obtained for the different samples
            for (i=0; i<number_of_presamples; ++i)
            {
                gsl_qrng_get(Sobol_generator,real_sample);
                maximal_allowed_deformation_parameters(temp_deformation_parameters, real_sample, real_parameters, complex_parameters, &result_info);
                for (j=0; j<number_of_integration_variables; ++j)
                    if (minimum <= temp_deformation_parameters[j] && temp_deformation_parameters[j] <= maximum)
                    {
                        if (optimized_deformation_parameters[j] > temp_deformation_parameters[j])
                            optimized_deformation_parameters[j] = temp_deformation_parameters[j];
                    } else if (temp_deformation_parameters[j] < minimum) {
                        optimized_deformation_parameters[j] = minimum;
                    }
            };
            result_info.process_errors();

            // reinitialize the Sobol sequence to obtain the same samples again
            gsl_qrng_free(Sobol_generator);
            Sobol_generator = gsl_qrng_alloc(gsl_qrng_sobol, number_of_integration_variables);

            // perform the sign check for each sample; decrease the "optimized_deformation_parameters" if necessary
            for (i=0; i<number_of_presamples; ++i)
            {
                gsl_qrng_get(Sobol_generator,real_sample);
                while ( !contour_deformation_polynomial_passes_sign_check(real_sample, real_parameters, complex_parameters, optimized_deformation_parameters, &result_info) )
                {
                    for (j=0; j<number_of_integration_variables; ++j)
                        optimized_deformation_parameters[j] *= decrease_factor;
                };
            };
            result_info.process_errors();

            // delete the quasi random number generator
            gsl_qrng_free(Sobol_generator);

            return optimized_deformation_parameter_vector;
        };

        // We want to bind the real, complex, and deformation parameters to the integrand.
        // These shared pointers can be used to avoid too early deallocation.
        std::shared_ptr<std::vector<real_t>> real_parameters;
        std::shared_ptr<std::vector<complex_t>> complex_parameters;
        std::shared_ptr<std::vector<real_t>> deformation_parameters;

        private: const std::vector<real_t> zeros; public:
        // function that performs the sign check for the contour deformation
        bool contour_deformation_polynomial_passes_sign_check
        (
            real_t const * const integration_variables,
            real_t const * const real_parameters,
            complex_t const * const complex_parameters,
            real_t const * const deformation_parameters,
            ResultInfo* result_info
        ) const
        {
            return contour_deformation_polynomial(integration_variables, real_parameters, complex_parameters, deformation_parameters, result_info).imag()
                <= contour_deformation_polynomial(integration_variables, real_parameters, complex_parameters, zeros.data()          , result_info).imag();
        }

        // "integrand" must be a member function, otherwise we cannot bind the struct
        complex_t integrand (
                                real_t const * const integration_variables,
                                real_t const * const real_parameters,
                                complex_t const * const complex_parameters,
                                real_t const * const deformation_parameters,
                                ResultInfo* result_info
                            ) const
        {
            // the required sign checks are performed inside the integrand for higher performance
            try {
                return deformed_integrand(integration_variables, real_parameters, complex_parameters, deformation_parameters, result_info);
            } catch (const sign_check_error& error) {
                // rethrow but with proper error message
                auto error_message = "Contour deformation in sector \"" + std::to_string(sector_id);
                error_message += "\", order { ";
                for (auto order : orders)
                    error_message += std::to_string(order) + " ";
                error_message += "}";
                error_message += error.what();
                error_message += " yields the wrong sign of \"contour_deformation_polynomial.imag\" or";
                error_message += " \"positive_polynomial.real\". Choose a larger \"number_of_presamples\" in";
                error_message += " \"optimize_deformation_parameters\" or decrease the \"deformation_parameters\".";
                throw sign_check_error(error_message);
            }
        };

        // constructor
        SectorContainerWithDeformation
        (
            const unsigned sector_id,
            const std::vector<int> orders,
            const unsigned number_of_integration_variables,
            DeformedIntegrandFunction * const deformed_integrand,
            #ifdef SECDEC_WITH_CUDA
                GetDeviceDeformedIntegrandFunction * const get_device_deformed_integrand,
            #endif
            DeformedIntegrandFunction * const contour_deformation_polynomial,
            MaximalDeformationFunction * const maximal_allowed_deformation_parameters
        ) :
        sector_id(sector_id),
        orders(orders),
        number_of_integration_variables(number_of_integration_variables),
        deformed_integrand(deformed_integrand),
        #ifdef SECDEC_WITH_CUDA
            get_device_deformed_integrand(get_device_deformed_integrand),
        #endif
        contour_deformation_polynomial(contour_deformation_polynomial),
        maximal_allowed_deformation_parameters(maximal_allowed_deformation_parameters),
        zeros(number_of_integration_variables,0)
        {};
    };

    /*
     * Conversion functions from "secdecutil::SectorContainer"
     * to "secdecutil::IntegrandContainer".
     * These functions return a function that binds the "real_parameters" and
     * the "complex_parameters", and optimize the "deformation_parameters"
     * (if applicable).
     */
    // SectorContainerWithDeformation -> IntegrandContainer
    template<typename real_t, typename complex_t>
    std::function<secdecutil::IntegrandContainer<complex_t, real_t const * const, real_t>(secdecutil::SectorContainerWithDeformation<real_t,complex_t>)>
    SectorContainerWithDeformation_to_IntegrandContainer(const std::vector<real_t>& real_parameters, const std::vector<complex_t>& complex_parameters,
                                                         unsigned number_of_presamples = 100000, real_t deformation_parameters_maximum = 1.,
                                                         real_t deformation_parameters_minimum = 1.e-5, real_t deformation_parameters_decrease_factor = 0.9)
    {
        auto shared_real_parameters = std::make_shared<std::vector<real_t>>(real_parameters);
        auto shared_complex_parameters = std::make_shared<std::vector<complex_t>>(complex_parameters);

        return
        [ = ]
        (secdecutil::SectorContainerWithDeformation<real_t,complex_t> sector_container)
        {
            sector_container.real_parameters = shared_real_parameters;
            sector_container.complex_parameters = shared_complex_parameters;

            sector_container.deformation_parameters =
                std::make_shared<std::vector<real_t>>
                (
                    sector_container.optimize_deformation_parameters
                    (
                        sector_container.real_parameters->data(),
                        sector_container.complex_parameters->data(),
                        number_of_presamples,
                        deformation_parameters_maximum,
                        deformation_parameters_minimum,
                        deformation_parameters_decrease_factor
                    )
                );


            std::function<complex_t(real_t const * const, const real_t*, ResultInfo*)> integrand = [sector_container](real_t const * const Args, const real_t* Pars, ResultInfo* result_info){
                return sector_container.integrand(Args,sector_container.real_parameters->data(), sector_container.complex_parameters->data(), Pars, result_info);
            };

            auto deformation_parameters = std::vector<std::vector<real_t>>({*sector_container.deformation_parameters});

            auto integrand_container = secdecutil::IntegrandContainer<complex_t, real_t const * const, real_t>
                                            (sector_container.number_of_integration_variables, integrand, deformation_parameters );

            integrand_container.extra_parameters = std::vector<std::vector<real_t>>({{deformation_parameters_minimum,deformation_parameters_decrease_factor}});

            integrand_container.display_name = "sector_"+std::to_string(sector_container.sector_id)+"_order";
            for(int i = 0; i < sector_container.orders.size(); i++){
                integrand_container.display_name += "_"+std::to_string(sector_container.orders[i]);
            }

            return integrand_container;
        };
    };

    // SectorContainerWithoutDeformation -> IntegrandContainer
    template<typename integrand_return_t, typename real_t, typename complex_t>
    std::function<secdecutil::IntegrandContainer<integrand_return_t, real_t const * const>(secdecutil::SectorContainerWithoutDeformation<real_t,complex_t,integrand_return_t>)>
    SectorContainerWithoutDeformation_to_IntegrandContainer(const std::vector<real_t>& real_parameters, const std::vector<complex_t>& complex_parameters)
    {
        auto shared_real_parameters = std::make_shared<std::vector<real_t>>(real_parameters);
        auto shared_complex_parameters = std::make_shared<std::vector<complex_t>>(complex_parameters);

        return [ shared_real_parameters, shared_complex_parameters ] (secdecutil::SectorContainerWithoutDeformation<real_t,complex_t,integrand_return_t> sector_container)
        {
            sector_container.real_parameters = shared_real_parameters;
            sector_container.complex_parameters = shared_complex_parameters;

            auto integrand = std::bind(&secdecutil::SectorContainerWithoutDeformation<real_t,complex_t,integrand_return_t>::integrand, sector_container,
                                       std::placeholders::_1, sector_container.real_parameters->data(), sector_container.complex_parameters->data(), std::placeholders::_2);

            auto integrand_container = secdecutil::IntegrandContainer<integrand_return_t, real_t const * const>(sector_container.number_of_integration_variables, integrand );

            integrand_container.display_name = "sector_"+std::to_string(sector_container.sector_id)+"_order";
            for(int i = 0; i < sector_container.orders.size(); i++){
                integrand_container.display_name += "_"+std::to_string(sector_container.orders[i]);
            }

            return integrand_container;
        };
    };


    #ifdef SECDEC_WITH_CUDA

        /*
         * Integrand containter that is usable on a CUDA device.
         */

        template<typename real_t, typename complex_t, typename integrand_return_t, unsigned long long maximal_number_of_functions, size_t number_of_real_parameters, size_t number_of_complex_parameters, char... name>
        struct CudaIntegrandContainerWithoutDeformation {

            // the call signature of an integrand
            using IntegrandFunction = typename SectorContainerWithoutDeformation<real_t,complex_t,integrand_return_t>::IntegrandFunction;

            // the call signature of the get device integrand function
            using GetDeviceIntegrandFunction = typename SectorContainerWithoutDeformation<real_t,complex_t,integrand_return_t>::GetDeviceIntegrandFunction;

            // member fields
            unsigned number_of_integration_variables;
            unsigned long long number_of_functions;
            GetDeviceIntegrandFunction* get_device_functions[maximal_number_of_functions];
            IntegrandFunction* device_functions[maximal_number_of_functions];
            IntegrandFunction* host_functions[maximal_number_of_functions];
            real_t real_parameters[number_of_real_parameters];
            complex_t complex_parameters[number_of_complex_parameters];
            bool call_get_device_functions_on_copy; // must call the slow function "get_device_functions" only after all other setup is finished
            std::string display_name = "INTEGRAND";
            std::shared_ptr<ResultInfo> result_info;
            ResultInfo* result_info_device;

            std::vector<std::vector<real_t*>> get_parameters() {return std::vector<std::vector<real_t*>>();}
            std::vector<std::vector<real_t>> get_extra_parameters() {return std::vector<std::vector<real_t>>();}

            // constructor
            CudaIntegrandContainerWithoutDeformation
            (
                unsigned number_of_integration_variables=0,
                unsigned long long number_of_functions=0,
                GetDeviceIntegrandFunction *const in_get_device_functions[]=nullptr,
                IntegrandFunction *const in_host_functions[]=nullptr,
                real_t const*const in_real_parameters=nullptr,
                complex_t const*const in_complex_parameters=nullptr,
                bool in_call_get_device_functions_on_copy=false
            ):
            number_of_integration_variables(number_of_integration_variables), number_of_functions(number_of_functions),
            get_device_functions(), host_functions(), real_parameters(), complex_parameters(), call_get_device_functions_on_copy(in_call_get_device_functions_on_copy)
            {
                if (number_of_functions > maximal_number_of_functions)
                    throw std::out_of_range
                    (
                        "Constructing CudaIntegrandContainerWithoutDeformation (with maximal_number_of_functions=" +
                        std::to_string(maximal_number_of_functions) + ") with too many (" +
                        std::to_string(number_of_functions) + ") functions."
                    );
                int i;
                for (unsigned long long k=0; k<number_of_functions; ++k)
                {
                   get_device_functions[k] = in_get_device_functions[k];
                   if(call_get_device_functions_on_copy) device_functions[k] = get_device_functions[k]();
                   host_functions[k] = in_host_functions[k];
                }
                if (in_real_parameters != nullptr)
                    for (i=0; i<number_of_real_parameters; ++i)
                       real_parameters[i] = in_real_parameters[i];
                if (in_complex_parameters != nullptr)
                    for (i=0; i<number_of_complex_parameters; ++i)
                       complex_parameters[i] = in_complex_parameters[i];

                result_info = std::make_shared<ResultInfo>();

                // create a pointer for the GPU to use
                auto error = cudaMallocManaged((void**)&result_info_device,sizeof(ResultInfo));
                if (error != cudaSuccess)
                    throw cuda_error(std::string(cudaGetErrorString(error)));
                memset(result_info_device, 0, sizeof(ResultInfo));
            }

            // destructor
            ~CudaIntegrandContainerWithoutDeformation()
            {
                // if result_info will be destroyed, also destruct result_info_device
                if(result_info.use_count()==1)
                {
                    auto error = cudaFree(result_info_device);
                    //if (error != cudaSuccess)
                    //    throw cuda_error(std::string(cudaGetErrorString(error)));
                }
            }

            // copy constructor
            CudaIntegrandContainerWithoutDeformation
            (
                const CudaIntegrandContainerWithoutDeformation& other
            ):
            CudaIntegrandContainerWithoutDeformation
            (
                other.number_of_integration_variables,
                other.number_of_functions,
                other.get_device_functions,
                other.host_functions,
                other.real_parameters,
                other.complex_parameters,
                other.call_get_device_functions_on_copy
            )
            {
                display_name = other.display_name;
                result_info = other.result_info;
                result_info_device = other.result_info_device;
            }

            // converting constructor
            template
            <unsigned long long other_maximal_number_of_functions, char... other_name>
            CudaIntegrandContainerWithoutDeformation
            (
                const CudaIntegrandContainerWithoutDeformation
                      <
                          real_t, complex_t, integrand_return_t, other_maximal_number_of_functions,
                          number_of_real_parameters, number_of_complex_parameters, other_name...
                      >& other
            ) :
            number_of_integration_variables(other.number_of_integration_variables), number_of_functions(other.number_of_functions), display_name(other.display_name),
            get_device_functions(), host_functions(), real_parameters(), complex_parameters(), call_get_device_functions_on_copy(other.call_get_device_functions_on_copy),
            result_info(other.result_info), result_info_device(other.result_info_device)
            {
                if (number_of_functions > maximal_number_of_functions)
                    throw std::out_of_range
                    (
                        "Constructing CudaIntegrandContainerWithoutDeformation (with maximal_number_of_functions=" +
                        std::to_string(maximal_number_of_functions) + ") with too many (" +
                        std::to_string(number_of_functions) + ") functions."
                    );
                int i;
                for (unsigned long long k=0; k<number_of_functions; ++k)
                {
                   get_device_functions[k] = other.get_device_functions[k];
                   if(call_get_device_functions_on_copy) device_functions[k] = get_device_functions[k]();
                   host_functions[k] = other.host_functions[k];
                }
                for (i=0; i<number_of_real_parameters; ++i)
                   real_parameters[i] = other.real_parameters[i];
                for (i=0; i<number_of_complex_parameters; ++i)
                   complex_parameters[i] = other.complex_parameters[i];
            }

            // call operator
            __host__ __device__
            integrand_return_t operator()(real_t const*const integration_variables) const
            {
                integrand_return_t res = 0;
                for (unsigned long long k=0; k<number_of_functions; ++k)
                {
                    #ifdef __CUDA_ARCH__
                        res += device_functions[k](integration_variables, real_parameters, complex_parameters, result_info_device);
                    #else
                        res += host_functions[k](integration_variables, real_parameters, complex_parameters, result_info.get());
                    #endif
                }
                return res;
            }

            void process_errors() const{
                result_info_device->process_errors();
                result_info->process_errors();
            }

            void clear_errors(){
                result_info_device->clear_errors();
                result_info->clear_errors();
            }

            template
            <unsigned long long other_maximal_number_of_functions, char... other_name>
            CudaIntegrandContainerWithoutDeformation& operator+=
            (
                const CudaIntegrandContainerWithoutDeformation
                      <
                          real_t, complex_t, integrand_return_t, other_maximal_number_of_functions,
                          number_of_real_parameters, number_of_complex_parameters, other_name...
                      >& other
            )
            {
                if (number_of_functions+other.number_of_functions > maximal_number_of_functions)
                    throw std::out_of_range
                    (
                        "Constructing CudaIntegrandContainerWithoutDeformation (with maximal_number_of_functions=" +
                        std::to_string(maximal_number_of_functions) + ") with too many (" +
                        std::to_string(number_of_functions) + ") functions."
                    );
                number_of_integration_variables = std::max(this->number_of_integration_variables,other.number_of_integration_variables);
                int i;
                unsigned long long offset = this->number_of_functions;
                for (unsigned long long k=0; k<other.number_of_functions; ++k)
                {
                    get_device_functions[k+offset] = other.get_device_functions[k];
                    host_functions[k+offset] = other.host_functions[k];
                }
                for (i=0; i<number_of_real_parameters; ++i)
                    real_parameters[i] = other.real_parameters[i];
                for (i=0; i<number_of_complex_parameters; ++i)
                    complex_parameters[i] = other.complex_parameters[i];
                number_of_functions = this->number_of_functions+other.number_of_functions;
                return *this;
            }
            template
            <unsigned long long other_maximal_number_of_functions, char... other_name>
            friend CudaIntegrandContainerWithoutDeformation operator+
            (
                const CudaIntegrandContainerWithoutDeformation& ic1,
                const CudaIntegrandContainerWithoutDeformation
                      <
                          real_t, complex_t, integrand_return_t, other_maximal_number_of_functions,
                          number_of_real_parameters, number_of_complex_parameters, other_name...
                      >& ic2
            )
            {
                CudaIntegrandContainerWithoutDeformation out(ic1);
                out += ic2;
                return out;
            }

        };

        template
        <typename real_t, typename complex_t, unsigned long long maximal_number_of_functions,
        size_t maximal_number_of_integration_variables, size_t number_of_real_parameters,
        size_t number_of_complex_parameters, char... name>
        struct CudaIntegrandContainerWithDeformation {

            // the call signature of an integrand
            using IntegrandFunction = typename SectorContainerWithDeformation<real_t,complex_t>::DeformedIntegrandFunction;

            // the call signature of the get device deformed integrand function
            using GetDeviceIntegrandFunction = typename SectorContainerWithDeformation<real_t,complex_t>::GetDeviceDeformedIntegrandFunction;

            using integrand_return_t = complex_t;

            // member fields
            unsigned number_of_integration_variables;
            unsigned long long number_of_functions;
            GetDeviceIntegrandFunction* get_device_functions[maximal_number_of_functions];
            IntegrandFunction* device_functions[maximal_number_of_functions];
            IntegrandFunction* host_functions[maximal_number_of_functions];
            real_t real_parameters[number_of_real_parameters];
            complex_t complex_parameters[number_of_complex_parameters];
            real_t deformation_parameters[maximal_number_of_functions][maximal_number_of_integration_variables];
            std::vector<std::vector<real_t>> extra_parameters;
            std::vector<std::vector<real_t*>> deformation_parameters_ptrs;
            bool call_get_device_functions_on_copy;
            std::string display_name = "INTEGRAND";
            std::shared_ptr<ResultInfo> result_info;
            ResultInfo* result_info_device;

            auto get_parameters() -> decltype(deformation_parameters_ptrs){return deformation_parameters_ptrs;}
            std::vector<std::vector<real_t>> get_extra_parameters(){return extra_parameters;}

            // constructor
            CudaIntegrandContainerWithDeformation
            (
                unsigned number_of_integration_variables=0,
                unsigned long long number_of_functions=0,
                GetDeviceIntegrandFunction *const in_get_device_functions[]=nullptr,
                IntegrandFunction *const in_host_functions[]=nullptr,
                real_t const*const in_real_parameters=nullptr,
                complex_t const*const in_complex_parameters=nullptr,
                real_t const*const in_deformation_parameters[maximal_number_of_integration_variables]=nullptr,
                bool in_call_get_device_functions_on_copy=false
            ):
            number_of_integration_variables(number_of_integration_variables), number_of_functions(number_of_functions),
            get_device_functions(), host_functions(), real_parameters(), complex_parameters(), deformation_parameters(), call_get_device_functions_on_copy(in_call_get_device_functions_on_copy),
            deformation_parameters_ptrs(), extra_parameters()
            {
                if (number_of_functions > maximal_number_of_functions)
                    throw std::out_of_range
                    (
                        "Constructing CudaIntegrandContainerWithDeformation (with maximal_number_of_functions=" +
                        std::to_string(maximal_number_of_functions) + ") with too many (" +
                        std::to_string(number_of_functions) + ") functions."
                    );
                int i;
                for (unsigned long long k=0; k<number_of_functions; ++k)
                {
                    get_device_functions[k] = in_get_device_functions[k];
                    if(call_get_device_functions_on_copy) device_functions[k] = get_device_functions[k]();
                    host_functions[k] = in_host_functions[k];
                    if (in_deformation_parameters != nullptr){
                        deformation_parameters_ptrs.push_back(std::vector<real_t*>());
                        for (i=0; i<number_of_integration_variables; ++i){
                            deformation_parameters[k][i] = in_deformation_parameters[k][i];
                            deformation_parameters_ptrs.back().push_back(&deformation_parameters[k][i]);
                        }
                    }
                }
                if (in_real_parameters != nullptr)
                    for (i=0; i<number_of_real_parameters; ++i)
                        real_parameters[i] = in_real_parameters[i];
                if (in_complex_parameters != nullptr)
                    for (i=0; i<number_of_complex_parameters; ++i)
                        complex_parameters[i] = in_complex_parameters[i];

                result_info = std::make_shared<ResultInfo>();

                // create a pointer for the GPU to use
                auto error = cudaMallocManaged((void**)&result_info_device,sizeof(ResultInfo));
                if (error != cudaSuccess)
                    throw cuda_error(std::string(cudaGetErrorString(error)));
                memset(result_info_device, 0, sizeof(ResultInfo));
            }
            
            // destructor
            ~CudaIntegrandContainerWithDeformation()
            {
                // if result_info will be destroyed, also destruct result_info_device
                if(result_info.use_count()==1)
                {
                    auto error = cudaFree(result_info_device);
                    //if (error != cudaSuccess)
                    //    throw cuda_error(std::string(cudaGetErrorString(error)));
                }
            }

            // copy constructor
            // Note: can not just call constructor as this would force us to convert [][] to **
            CudaIntegrandContainerWithDeformation
            (
                const CudaIntegrandContainerWithDeformation& other
            ):
            number_of_integration_variables(other.number_of_integration_variables), number_of_functions(other.number_of_functions), display_name(other.display_name),
            get_device_functions(), host_functions(), real_parameters(), complex_parameters(), deformation_parameters(), call_get_device_functions_on_copy(other.call_get_device_functions_on_copy),
            deformation_parameters_ptrs(), extra_parameters(other.extra_parameters), result_info(other.result_info), result_info_device(other.result_info_device)
            {
                if (number_of_functions > maximal_number_of_functions)
                    throw std::out_of_range
                    (
                        "Constructing CudaIntegrandContainerWithDeformation (with maximal_number_of_functions=" +
                         std::to_string(maximal_number_of_functions) + ") with too many (" +
                         std::to_string(number_of_functions) + ") functions."
                    );
                int i;
                for (unsigned long long k=0; k<number_of_functions; ++k)
                {
                    get_device_functions[k] = other.get_device_functions[k];
                    if(call_get_device_functions_on_copy) device_functions[k] = get_device_functions[k]();
                    host_functions[k] = other.host_functions[k];
                    if (other.deformation_parameters != nullptr){
                        deformation_parameters_ptrs.push_back(std::vector<real_t*>());
                        for (i=0; i<number_of_integration_variables; ++i){
                            deformation_parameters[k][i] = other.deformation_parameters[k][i];
                            deformation_parameters_ptrs.back().push_back(&deformation_parameters[k][i]);
                        }
                    }
                }
                if (other.real_parameters != nullptr)
                    for (i=0; i<number_of_real_parameters; ++i)
                        real_parameters[i] = other.real_parameters[i];
                if (other.complex_parameters != nullptr)
                    for (i=0; i<number_of_complex_parameters; ++i)
                        complex_parameters[i] = other.complex_parameters[i];
            }

            // converting constructor
            template
            <unsigned long long other_maximal_number_of_functions, char... other_name>
            CudaIntegrandContainerWithDeformation
            (
                const CudaIntegrandContainerWithDeformation
                      <
                          real_t, complex_t, other_maximal_number_of_functions,
                          maximal_number_of_integration_variables, number_of_real_parameters,
                          number_of_complex_parameters, other_name...
                      >& other
            ) :
            number_of_integration_variables(other.number_of_integration_variables), number_of_functions(other.number_of_functions), display_name(other.display_name),
            get_device_functions(), host_functions(), real_parameters(), complex_parameters(), call_get_device_functions_on_copy(other.call_get_device_functions_on_copy),
            deformation_parameters_ptrs(), extra_parameters(other.extra_parameters), result_info(other.result_info), result_info_device(other.result_info_device)
            {
                if (number_of_functions > maximal_number_of_functions)
                    throw std::out_of_range
                    (
                        "Constructing CudaIntegrandContainerWithDeformation (with maximal_number_of_functions=" +
                        std::to_string(maximal_number_of_functions) + ") with too many (" +
                        std::to_string(number_of_functions) + ") functions."
                    );
                int i;
                for (unsigned long long k=0; k<number_of_functions; ++k)
                {
                    get_device_functions[k] = other.get_device_functions[k];
                    if(call_get_device_functions_on_copy) device_functions[k] = get_device_functions[k]();
                    host_functions[k] = other.host_functions[k];
                    deformation_parameters_ptrs.push_back(std::vector<real_t*>());
                    for (i=0; i<number_of_integration_variables; ++i){
                        deformation_parameters[k][i] = other.deformation_parameters[k][i];
                        deformation_parameters_ptrs.back().push_back(&deformation_parameters[k][i]);
                    }
                }
                for (i=0; i<number_of_real_parameters; ++i)
                    real_parameters[i] = other.real_parameters[i];
                for (i=0; i<number_of_complex_parameters; ++i)
                    complex_parameters[i] = other.complex_parameters[i];
            }

            // call operator
            __host__ __device__
            integrand_return_t operator()(real_t const*const integration_variables) const
            {
                integrand_return_t res = 0;
                for (unsigned long long k=0; k<number_of_functions; ++k)
                {
                    #ifdef __CUDA_ARCH__
                        res += device_functions[k](integration_variables, real_parameters, complex_parameters, deformation_parameters[k], result_info_device);
                    #else
                        res += host_functions[k](integration_variables, real_parameters, complex_parameters, deformation_parameters[k], result_info.get());
                    #endif
                }
                return res;
            }

            void process_errors() const{
                result_info_device->process_errors();
                result_info->process_errors();
            }

            void clear_errors(){
                result_info_device->clear_errors();
                result_info->clear_errors();
            }

            template
            <unsigned long long other_maximal_number_of_functions, char... other_name>
            CudaIntegrandContainerWithDeformation& operator+=
            (
                const CudaIntegrandContainerWithDeformation
                      <
                          real_t,complex_t,other_maximal_number_of_functions,
                          maximal_number_of_integration_variables, number_of_real_parameters,
                          number_of_complex_parameters, other_name...
                      >& other
            )
            {
                if (number_of_functions+other.number_of_functions > maximal_number_of_functions)
                    throw std::out_of_range
                    (
                        "Constructing CudaIntegrandContainerWithDeformation (with maximal_number_of_functions=" +
                        std::to_string(maximal_number_of_functions) + ") with too many (" +
                        std::to_string(number_of_functions) + ") functions."
                    );
                number_of_integration_variables = std::max(this->number_of_integration_variables,other.number_of_integration_variables);
                int i;
                unsigned long long offset = this->number_of_functions;
                for (unsigned long long k=0; k<other.number_of_functions; ++k)
                {
                    get_device_functions[k+offset] = other.get_device_functions[k];
                    host_functions[k+offset] = other.host_functions[k];
                    deformation_parameters_ptrs.push_back(std::vector<real_t*>());
                    for (i=0; i<other.number_of_integration_variables; ++i){
                        deformation_parameters[k+offset][i] = other.deformation_parameters[k][i];
                        deformation_parameters_ptrs.back().push_back(&deformation_parameters[k+offset][i]);
                    }
                }
                for (i=0; i<number_of_real_parameters; ++i)
                    real_parameters[i] = other.real_parameters[i];
                for (i=0; i<number_of_complex_parameters; ++i)
                    complex_parameters[i] = other.complex_parameters[i];
                number_of_functions = this->number_of_functions+other.number_of_functions;
                return *this;
            }
            template
            <unsigned long long other_maximal_number_of_functions, char... other_name>
            friend CudaIntegrandContainerWithDeformation operator+
            (
                const CudaIntegrandContainerWithDeformation& ic1,
                const CudaIntegrandContainerWithDeformation
                      <
                          real_t,complex_t,other_maximal_number_of_functions,
                          maximal_number_of_integration_variables, number_of_real_parameters,
                          number_of_complex_parameters, other_name...
                      >& ic2
            )
            {
                CudaIntegrandContainerWithDeformation out(ic1);
                out += ic2;
                return out;
            }

        };

        /*
         * Conversion functions from "secdecutil::SectorContainer"
         * to "secdecutil::CudaIntegrandContainer".
         * These functions return a function that "binds" the "real_parameters" and
         * the "complex_parameters", and optimize the "deformation_parameters"
         * (if applicable) without using "std::bind".
         */
        // SectorContainerWithDeformation -> CudaIntegrandContainer
        template<size_t maximal_number_of_integration_variables, size_t number_of_real_parameters, size_t number_of_complex_parameters, char... name, typename real_t, typename complex_t>
        std::function
        <
            CudaIntegrandContainerWithDeformation<real_t,complex_t,1/*maximal_number_of_functions*/,maximal_number_of_integration_variables,number_of_real_parameters,number_of_complex_parameters,name...>
            (secdecutil::SectorContainerWithDeformation<real_t,complex_t>)
        >
        SectorContainerWithDeformation_to_CudaIntegrandContainer(const std::vector<real_t>& real_parameters, const std::vector<complex_t>& complex_parameters,
                                                                 unsigned number_of_presamples = 100000, real_t deformation_parameters_maximum = 1.,
                                                                 real_t deformation_parameters_minimum = 1.e-5, real_t deformation_parameters_decrease_factor = 0.9)
        {
            return
            [ = ]
            (secdecutil::SectorContainerWithDeformation<real_t,complex_t> sector_container)
            {
                std::vector<real_t> optimized_deformation_parameters =
                    sector_container.optimize_deformation_parameters
                        (
                            real_parameters.data(),
                            complex_parameters.data(),
                            number_of_presamples,
                            deformation_parameters_maximum,
                            deformation_parameters_minimum,
                            deformation_parameters_decrease_factor
                        );
                real_t const*const optimized_deformation_parameters_ptr = optimized_deformation_parameters.data();

                auto integrand_container =
                CudaIntegrandContainerWithDeformation<real_t,complex_t,1/*maximal_number_of_functions*/,maximal_number_of_integration_variables,number_of_real_parameters,number_of_complex_parameters,name...>
                {
                    sector_container.number_of_integration_variables,
                    1, // number_of_functions
                    &sector_container.get_device_deformed_integrand,
                    &sector_container.deformed_integrand, // host function
                    real_parameters.data(),
                    complex_parameters.data(),
                    &optimized_deformation_parameters_ptr
                };

                integrand_container.extra_parameters = std::vector<std::vector<real_t>>{{deformation_parameters_minimum,deformation_parameters_decrease_factor}};

                integrand_container.display_name = "sector_"+std::to_string(sector_container.sector_id)+"_order";
                for(int i = 0; i < sector_container.orders.size(); i++){
                    integrand_container.display_name += "_"+std::to_string(sector_container.orders[i]);
                }

                return integrand_container;
            };
        };

        // SectorContainerWithoutDeformation -> CudaIntegrandContainer
        template<typename integrand_return_t, size_t number_of_real_parameters, size_t number_of_complex_parameters, char... name, typename real_t, typename complex_t>
        std::function
        <
            CudaIntegrandContainerWithoutDeformation<real_t,complex_t,integrand_return_t,1/*maximal_number_of_functions*/,number_of_real_parameters,number_of_complex_parameters,name...>
            (secdecutil::SectorContainerWithoutDeformation<real_t,complex_t,integrand_return_t>)
        >
        SectorContainerWithoutDeformation_to_CudaIntegrandContainer(const std::vector<real_t>& real_parameters, const std::vector<complex_t>& complex_parameters)
        {

            return [ real_parameters, complex_parameters ] (secdecutil::SectorContainerWithoutDeformation<real_t,complex_t,integrand_return_t> sector_container)
            {
                auto integrand_container =
                CudaIntegrandContainerWithoutDeformation<real_t,complex_t,integrand_return_t,1/*maximal_number_of_functions*/,number_of_real_parameters,number_of_complex_parameters,name...>
                {
                    sector_container.number_of_integration_variables,
                    1, // number_of_functions
                    &sector_container.get_device_undeformed_integrand,
                    &sector_container.undeformed_integrand, // host function
                    real_parameters.data(),
                    complex_parameters.data()
                };

                integrand_container.display_name = "sector_"+std::to_string(sector_container.sector_id)+"_order";
                for(int i = 0; i < sector_container.orders.size(); i++){
                    integrand_container.display_name += "_"+std::to_string(sector_container.orders[i]);
                }

                return integrand_container;
            };
        };

    #endif

}

#endif

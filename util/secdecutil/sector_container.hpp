#ifndef SecDecUtil_sector_container_hpp_included
#define SecDecUtil_sector_container_hpp_included

#include <exception>
#include <functional>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <secdecutil/integrand_container.hpp>

namespace secdecutil {

    // the return type of the integral transformation (contour deformation)
    template<typename complex_t>
    struct integral_transformation_t
    {
        std::vector<complex_t> transformed_variables;
        complex_t Jacobian_determinant;
    };

    // this error is thrown if the sign check of the deformation (contour_deformation_polynomial.imag() <= 0) fails
    struct sign_check_error : public std::runtime_error { using std::runtime_error::runtime_error; };


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
         complex_t const * const complex_parameters
         );

        const unsigned sector_id;
        const std::vector<int> orders;
        const unsigned number_of_integration_variables;
        IntegrandFunction * const undeformed_integrand;

        std::shared_ptr<std::vector<real_t>> real_parameters;
        std::shared_ptr<std::vector<complex_t>> complex_parameters;
        // "integrand" must be a member function, otherwise we cannot bind the struct
        integrand_return_t integrand
        (
            real_t const * const integration_variables,
            real_t const * const real_parameters,
            complex_t const * const complex_parameters
        )
        {
            return undeformed_integrand(integration_variables, real_parameters, complex_parameters);
        };

        // constructor
        SectorContainerWithoutDeformation
        (
            const unsigned sector_id,
            const std::vector<int> orders,
            const unsigned number_of_integration_variables,
            IntegrandFunction * const undeformed_integrand
        ) :
        sector_id(sector_id),
        orders(orders),
        number_of_integration_variables(number_of_integration_variables),
        undeformed_integrand(undeformed_integrand)
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
         real_t const * const deformation_parameters
         );

        // the call signature of the function that computes the maximal deformation parameters
        typedef void MaximalDeformationFunction
        (
         real_t * output_deformation_parameters,
         real_t const * const integration_variables,
         real_t const * const real_parameters,
         complex_t const * const complex_parameters
         );

        const unsigned sector_id;
        const std::vector<int> orders;
        const unsigned number_of_integration_variables;
        DeformedIntegrandFunction * const deformed_integrand;
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

            // Define a lambda function that generates 'number_of_integration_variables'-dimensional
            // uniformly-distributed points.
            // We use the c++11-builtin 'minstd_rand' and seed it with the 'number_of_presamples'.
            std::minstd_rand random_number_generator(/* seed  = */ number_of_presamples);
            std::uniform_real_distribution<real_t> uniform_distribution(0,1);
            auto generate_sample =
            [ &random_number_generator , &uniform_distribution , this ]
            (real_t * output)
            {
                for (unsigned k = 0 ; k < number_of_integration_variables ; ++k)
                    output[k] = uniform_distribution(random_number_generator);
            };

            // find the minimum of the lambdas obtained for the different samples
            for (i=0; i<number_of_presamples; ++i)
            {
                generate_sample(real_sample);
                maximal_allowed_deformation_parameters(temp_deformation_parameters, real_sample, real_parameters, complex_parameters);
                for (j=0; j<number_of_integration_variables; ++j)
                    if (minimum <= temp_deformation_parameters[j] && temp_deformation_parameters[j] <= maximum)
                    {
                        if (optimized_deformation_parameters[j] > temp_deformation_parameters[j])
                            optimized_deformation_parameters[j] = temp_deformation_parameters[j];
                    } else if (temp_deformation_parameters[j] < minimum) {
                        optimized_deformation_parameters[j] = minimum;
                    }
            };

            // reseed the random number generator to obtain the same samples again
            random_number_generator.seed(number_of_presamples);

            // perform the sign check for each sample; decrease the "optimized_deformation_parameters" if necessary
            integral_transformation_t<complex_t> deformation;
            for (i=0; i<number_of_presamples; ++i)
            {
                generate_sample(real_sample);
                while ( !contour_deformation_polynomial_passes_sign_check(real_sample, real_parameters, complex_parameters, optimized_deformation_parameters) )
                {
                    for (j=0; j<number_of_integration_variables; ++j)
                        optimized_deformation_parameters[j] *= decrease_factor;
                };
            };

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
            real_t const * const deformation_parameters
        ) const
        {
            return contour_deformation_polynomial(integration_variables, real_parameters, complex_parameters, deformation_parameters).imag()
                <= contour_deformation_polynomial(integration_variables, real_parameters, complex_parameters, zeros.data()          ).imag();
        }

        // "integrand" must be a member function, otherwise we cannot bind the struct
        complex_t integrand (
                                real_t const * const integration_variables,
                                real_t const * const real_parameters,
                                complex_t const * const complex_parameters,
                                real_t const * const deformation_parameters
                            ) const
        {
            // the required sign checks are performed inside the integrand for higher performance
            try {
                return deformed_integrand(integration_variables, real_parameters, complex_parameters, deformation_parameters);
            } catch (const sign_check_error& error) {
                // rethrow but with proper error message
                auto error_message = "Contour deformation in sector \"" + std::to_string(sector_id);
                error_message += "\", order { ";
                for (auto order : orders)
                    error_message += std::to_string(order) + " ";
                error_message += "} yields the wrong sign of \"contour_deformation_polynomial.imag\" or";
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
            DeformedIntegrandFunction * const contour_deformation_polynomial,
            MaximalDeformationFunction * const maximal_allowed_deformation_parameters
        ) :
        sector_id(sector_id),
        orders(orders),
        number_of_integration_variables(number_of_integration_variables),
        deformed_integrand(deformed_integrand),
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
    std::function<secdecutil::IntegrandContainer<complex_t, real_t const * const>(secdecutil::SectorContainerWithDeformation<real_t,complex_t>)>
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

            auto integrand = std::bind(&secdecutil::SectorContainerWithDeformation<real_t,complex_t>::integrand, sector_container,
                                       std::placeholders::_1, sector_container.real_parameters->data(), sector_container.complex_parameters->data(),
                                       sector_container.deformation_parameters->data());

            return secdecutil::IntegrandContainer<complex_t, real_t const * const>(sector_container.number_of_integration_variables, integrand);
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
                                       std::placeholders::_1, sector_container.real_parameters->data(), sector_container.complex_parameters->data());

            return secdecutil::IntegrandContainer<integrand_return_t, real_t const * const>(sector_container.number_of_integration_variables, integrand );
        };
    };

}

#endif

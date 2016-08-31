#ifndef SecDecUtil_sector_container_hpp_included
#define SecDecUtil_sector_container_hpp_included

#include <exception>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <gsl/gsl_qrng.h>
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
         complex_t const * const complex_parameters
         );

        const unsigned sector_id;
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
            const unsigned number_of_integration_variables,
            IntegrandFunction * const undeformed_integrand
        ) :
        sector_id(sector_id),
        number_of_integration_variables(number_of_integration_variables),
        undeformed_integrand(undeformed_integrand)
        {};
    };

    template<typename real_t, typename complex_t>
    struct SectorContainerWithDeformation
    {
        // the call signature of an integrand to be deformed
        typedef complex_t DeformableIntegrandFunction
        (
         complex_t const * const transformed_integration_variables,
         real_t const * const real_parameters,
         complex_t const * const complex_parameters
         );

        // the call signature of the integral transformation (contour deformation)
        typedef integral_transformation_t<complex_t> ContourDeformationFunction
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
        const unsigned number_of_integration_variables;
        DeformableIntegrandFunction * const undeformed_integrand;
        ContourDeformationFunction * const contour_deformation;
        DeformableIntegrandFunction * const contour_deformation_polynomial;
        MaximalDeformationFunction * const maximal_allowed_deformation_parameters;

        // the function that optimizes the deformation_parameters
        std::vector<real_t> optimize_deformation_parameters (
                                                                real_t const * const real_parameters,
                                                                complex_t const * const complex_parameters,
                                                                const unsigned number_of_samples = 100000,
                                                                const real_t maximum = 1.,
                                                                const real_t minimum = 1.e-5,
                                                                const real_t decrease_factor = 0.9
                                                            ) const
        {
            // define indices for the loops
            unsigned i,j;

            // if no sampling desired (number_of_samples == 0) set the deformation parameters to the maximum
            if (number_of_samples == 0)
                return std::vector<real_t>(number_of_integration_variables,maximum);

            // initialize the output, and temporary vectors
            std::vector<real_t> optimized_deformation_parameter_vector(number_of_integration_variables,0);
            std::vector<real_t> temp_deformation_parameter_vector(number_of_integration_variables,0);
            std::vector<real_t> real_sample_vector(number_of_integration_variables,0);
            std::vector<complex_t> complex_sample_vector(number_of_integration_variables,0);
            real_t * optimized_deformation_parameters = optimized_deformation_parameter_vector.data();
            real_t * temp_deformation_parameters = temp_deformation_parameter_vector.data();
            real_t * real_sample = real_sample_vector.data();
            complex_t * complex_sample = complex_sample_vector.data();

            // define a Sobol sequence using the gsl
            // Restriction to at most 40 dimensions only because of the implementation in the gsl. --> Use a different Sobol implementation if higher dimensionality is needed.
            int Sobol_maxdim = 40;
            if (number_of_integration_variables > Sobol_maxdim)
                throw gsl_error("The gsl implements Sobol sequences only up to " + std::to_string(Sobol_maxdim) +" dimensions (need " +
                                std::to_string(number_of_integration_variables) + "). Please set the \"deformation_parameters\" manually.");

            // define the generator
            gsl_qrng * Sobol_generator = gsl_qrng_alloc(gsl_qrng_sobol, number_of_integration_variables);

            // average over the lambdas obtained for the different samples
            for (i=0; i<number_of_samples; ++i)
            {
                gsl_qrng_get(Sobol_generator,real_sample);
                maximal_allowed_deformation_parameters(temp_deformation_parameters, real_sample, real_parameters, complex_parameters);
                for (j=0; j<number_of_integration_variables; ++j)
                    if (minimum <= temp_deformation_parameters[j] && temp_deformation_parameters[j] <= maximum)
                        optimized_deformation_parameters[j] += temp_deformation_parameters[j] / number_of_samples;
                    else if (temp_deformation_parameters[j] > maximum)
                        optimized_deformation_parameters[j] += maximum / number_of_samples;
                    else
                        optimized_deformation_parameters[j] += minimum / number_of_samples;
            };

            // reinitialize the Sobol sequence to obtain the same samples again
            gsl_qrng_free(Sobol_generator);
            Sobol_generator = gsl_qrng_alloc(gsl_qrng_sobol, number_of_integration_variables);

            // perform the sign check for each sample; decrease the "optimized_deformation_parameters" if necessary
            integral_transformation_t<complex_t> deformation;
            for (i=0; i<number_of_samples; ++i)
            {
                gsl_qrng_get(Sobol_generator,real_sample);

                // the "contour_deformation_polynomial" takes complex --> must cast the "sample"
                for (j=0; j<number_of_integration_variables; ++j)
                    complex_sample[j] = real_sample[j];

                deformation = contour_deformation(real_sample, real_parameters, complex_parameters, optimized_deformation_parameters);
                while (contour_deformation_polynomial(deformation.transformed_variables.data(), real_parameters, complex_parameters).imag() >
                       contour_deformation_polynomial(complex_sample, real_parameters, complex_parameters).imag())
                {
                    for (j=0; j<number_of_integration_variables; ++j)
                        optimized_deformation_parameters[j] *= decrease_factor;
                    deformation = contour_deformation(real_sample, real_parameters, complex_parameters, optimized_deformation_parameters);
                };
            };

            // delete the quasi random number generator
            gsl_qrng_free(Sobol_generator);

            return optimized_deformation_parameter_vector;
        };

        // We want to bind the real, complex, and deformation parameters to the integrand.
        // These shared pointers can be used to avoid too early deallocation.
        std::shared_ptr<std::vector<real_t>> real_parameters;
        std::shared_ptr<std::vector<complex_t>> complex_parameters;
        std::shared_ptr<std::vector<real_t>> deformation_parameters;
        complex_t integrand (
                                real_t const * const integration_variables,
                                real_t const * const real_parameters,
                                complex_t const * const complex_parameters,
                                real_t const * const deformation_parameters
                            ) const
        {
            auto deformation = contour_deformation(integration_variables, real_parameters, complex_parameters, deformation_parameters);

            auto untransformed_integration_variable_vector = std::vector<complex_t>(number_of_integration_variables);
            auto untransformed_integration_variables = untransformed_integration_variable_vector.data();
            for (unsigned i=0; i<number_of_integration_variables; ++i)
                untransformed_integration_variables[i] = integration_variables[i];

            if (contour_deformation_polynomial(deformation.transformed_variables.data(), real_parameters, complex_parameters).imag() >
                contour_deformation_polynomial(untransformed_integration_variables, real_parameters, complex_parameters).imag())
                throw sign_check_error("Contour deformation in sector \"" + std::to_string(sector_id) + "\" yields the wrong sign of \"contour_deformation_polynomial.imag\". Choose a larger \"number_of_samples\" in \"optimize_deformation_parameters\" (recommended) or decrease \"deformation_parameters\".");

            return deformation.Jacobian_determinant * undeformed_integrand(deformation.transformed_variables.data(), real_parameters, complex_parameters);
        };

        // constructor
        SectorContainerWithDeformation
        (
            const unsigned sector_id,
            const unsigned number_of_integration_variables,
            DeformableIntegrandFunction * const undeformed_integrand,
            ContourDeformationFunction * const contour_deformation,
            DeformableIntegrandFunction * const contour_deformation_polynomial,
            MaximalDeformationFunction * const maximal_allowed_deformation_parameters
        ) :
        sector_id(sector_id),
        number_of_integration_variables(number_of_integration_variables),
        undeformed_integrand(undeformed_integrand),
        contour_deformation(contour_deformation),
        contour_deformation_polynomial(contour_deformation_polynomial),
        maximal_allowed_deformation_parameters(maximal_allowed_deformation_parameters)
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
                                                         unsigned number_of_samples = 100000, real_t deformation_parameters_maximum = 1.,
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
                        number_of_samples,
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

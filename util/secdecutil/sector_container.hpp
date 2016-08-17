#ifndef SecDecUtil_sector_container_hpp_included
#define SecDecUtil_sector_container_hpp_included

#include <exception>
#include <functional>
#include <memory>
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

        // the call signature of the function to optimize deformation parameters (contour deformation) // TODO: do we need this?
        typedef void OptimizeDeformationFunction
        (
         real_t const * const initial_guess,
         real_t const * const real_parameters,
         complex_t const * const complex_parameters,
         const size_t number_of_samples
         );

        const unsigned sector_id;
        const unsigned number_of_integration_variables;
        DeformableIntegrandFunction * const undeformed_integrand;
        ContourDeformationFunction * const contour_deformation;
        DeformableIntegrandFunction * const contour_deformation_polynomial;

        std::vector<real_t> optimize_deformation_parameters (
                                                             real_t const * const initial_guess,
                                                             real_t const * const real_parameters,
                                                             complex_t const * const complex_parameters,
                                                             const size_t number_of_samples
                                                             ) const;

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
            if (contour_deformation_polynomial(deformation.transformed_variables.data(), real_parameters, complex_parameters).imag() > 0.)
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
            DeformableIntegrandFunction * const contour_deformation_polynomial
        ) :
        sector_id(sector_id),
        number_of_integration_variables(number_of_integration_variables),
        undeformed_integrand(undeformed_integrand),
        contour_deformation(contour_deformation),
        contour_deformation_polynomial(contour_deformation_polynomial)
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
                                          unsigned number_of_samples = 10000, real_t deformation_parameters_initial_guess = 1.)
    {
        auto shared_real_parameters = std::make_shared<std::vector<real_t>>(real_parameters);
        auto shared_complex_parameters = std::make_shared<std::vector<complex_t>>(complex_parameters);

        return
        [shared_real_parameters,shared_complex_parameters,number_of_samples,deformation_parameters_initial_guess]
        (secdecutil::SectorContainerWithDeformation<real_t,complex_t> sector_container)
        {
            sector_container.real_parameters = shared_real_parameters;
            sector_container.complex_parameters = shared_complex_parameters;

            sector_container.deformation_parameters = std::make_shared<std::vector<real_t>>(sector_container.number_of_integration_variables,deformation_parameters_initial_guess);
            // TODO call optimize lambda with "number_of_samples";

            auto integrand = std::bind(&secdecutil::SectorContainerWithDeformation<real_t,complex_t>::integrand, sector_container,
                                       std::placeholders::_1, sector_container.real_parameters->data(), sector_container.complex_parameters->data(),
                                       sector_container.deformation_parameters->data());

            return secdecutil::IntegrandContainer<complex_t, real_t const * const>(sector_container.number_of_integration_variables, integrand );
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

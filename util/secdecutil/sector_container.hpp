#ifndef SecDecUtil_sector_container_hpp_included
#define SecDecUtil_sector_container_hpp_included

#include <exception>
#include <functional>

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
        IntegrandFunction * const integrand;
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
    };
}

#endif

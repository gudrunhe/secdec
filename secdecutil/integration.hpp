#ifndef SecDecUtil_integration_hpp_included
#define SecDecUtil_integration_hpp_included

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

    // container class to collect the integrand functions
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

        // the call signature of the function to optimize deformation parameters (contour deformation)
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

    
    template <typename T, typename ...Args>
    class IntegrandContainer {

    private:
        enum Operation { add, subtract, multiply, divide };

    public:
        
        int number_of_integration_variables;
        std::function<T(Args...)> integrand;
        
        /*
         *  Helper functions
         */
        template<int operation>
        static IntegrandContainer add_subtract_multiply_or_divide(const IntegrandContainer& ic1, const IntegrandContainer& ic2)
        {
            int number_of_integration_variables = std::max(ic1.number_of_integration_variables, ic2.number_of_integration_variables);
            std::function<T(Args...)> integrand;
            
            if (operation == add) {
                integrand = [ic1, ic2] (Args... x) { return ic1.integrand(x...) + ic2.integrand(x...); };
            } else if (operation == subtract ) {
                integrand = [ic1, ic2] (Args... x) { return ic1.integrand(x...) - ic2.integrand(x...); };
            } else if (operation == multiply ) {
                integrand = [ic1, ic2] (Args... x) { return ic1.integrand(x...) * ic2.integrand(x...); };
            } else if ( operation == divide ) {
                integrand = [ic1, ic2] (Args... x) { return ic1.integrand(x...) / ic2.integrand(x...); };
            }

            return IntegrandContainer(number_of_integration_variables, integrand);
        };

        /*
         *  Compound assignment operators
         */
        IntegrandContainer& operator-=(const IntegrandContainer& ic1)
        {
            *this = *this - ic1;
            return *this;
        };
        
        IntegrandContainer& operator+=(const IntegrandContainer& ic1)
        {
            *this = *this + ic1;
            return *this;
        };

        // TODO - bof add tests
        IntegrandContainer& operator*=(const IntegrandContainer& ic1)
        {
            *this = *this * ic1;
            return *this;
        };
        IntegrandContainer& operator/=(const IntegrandContainer& ic1)
        {
            *this = *this / ic1;
            return *this;
        };
        // TODO - eof add tests

        /*
         *  Binary operators
         */
        friend IntegrandContainer operator-(const IntegrandContainer& ic1, const IntegrandContainer& ic2)
        {
            return add_subtract_multiply_or_divide<subtract>(ic1,ic2);
        };
        
        friend IntegrandContainer operator+(const IntegrandContainer& ic1, const IntegrandContainer& ic2)
        {
            return add_subtract_multiply_or_divide<add>(ic1,ic2);
        };

        // TODO - bof add tests
        friend IntegrandContainer operator*(const IntegrandContainer& ic1, const IntegrandContainer& ic2)
        {
            return add_subtract_multiply_or_divide<multiply>(ic1,ic2);
        };
        friend IntegrandContainer operator/(const IntegrandContainer& ic1, const IntegrandContainer& ic2)
        {
            return add_subtract_multiply_or_divide<divide>(ic1,ic2);
        };
        // TODO - eof add tests

        IntegrandContainer(const int number_of_integration_variables, const std::function<T(Args...)>& integrand):
        number_of_integration_variables (number_of_integration_variables), integrand(integrand)
        {};

        IntegrandContainer()
        {};
        
    };
    
}

#endif
#ifndef SecDecUtil_integrand_container_hpp_included
#define SecDecUtil_integrand_container_hpp_included

#include <complex>
#include <functional>

namespace secdecutil {

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
         * unary operators
         */
        IntegrandContainer operator+() const
        {
            return *this;
        };

        IntegrandContainer operator-() const
        {
            return IntegrandContainer
            (
                number_of_integration_variables,
                [this] (Args... x) { return - integrand(x...); }
            );
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

        friend IntegrandContainer operator*(const IntegrandContainer& ic1, const IntegrandContainer& ic2)
        {
            return add_subtract_multiply_or_divide<multiply>(ic1,ic2);
        };
        friend IntegrandContainer operator/(const IntegrandContainer& ic1, const IntegrandContainer& ic2)
        {
            return add_subtract_multiply_or_divide<divide>(ic1,ic2);
        };

        /*
         * Constructors
         */
        IntegrandContainer(const int number_of_integration_variables, const std::function<T(Args...)>& integrand):
        number_of_integration_variables (number_of_integration_variables), integrand(integrand)
        {};

        // default constructor (the "zero-integrand")
        IntegrandContainer() :
        number_of_integration_variables(0),integrand([](...){return T();})
        {};

    };

  template<typename T, typename... Args>
  IntegrandContainer<T, Args...> complex_to_real
  (IntegrandContainer<std::complex<T>, Args...> ic, T (*operation) (const std::complex<T>& ))
  {
    std::function<T(Args...)> new_integrand = [&ic, operation] (const Args... integration_variables){ return operation(ic.integrand(integration_variables...)); };
    return IntegrandContainer<T, Args...>(ic.number_of_integration_variables, new_integrand);
  }
  
}

#endif

#ifndef SecDecUtil_uncertainties_hpp_included
#define SecDecUtil_uncertainties_hpp_included

#include <cmath>
#include <complex>
#include <ostream>

#ifdef SECDEC_WITH_CUDA
    #include <thrust/complex.h>
#endif

namespace secdecutil {

    /*!
     * Implement the propagation of !!UNCORRELATED!!
     * standard deviations via the class
     * "UncorrelatedDeviation".
     *
     */

    /****************************************************************************************************************************************
     *                                                                                                                                      *
     *     .. note::                                                                                                                        *
     *     Division by :cpp:class:`UncorrelatedDeviation` is not implemented as this operation is not always well defined.                  *
     *     Specifically, it is ill defined in the case that the errors are Gaussian distributed as the expectation value,                   *
     *                                                                                                                                      *
     *     .. math::                                                                                                                        *
     *     \mathrm{E}\left[\frac{1}{X}\right] = \int_{-\infty}^{\infty} \frac{1}{X} p(X)\ \mathrm{d}X,                                      *
     *                                                                                                                                      *
     *     where                                                                                                                            *
     *                                                                                                                                      *
     *     .. math::                                                                                                                        *
     *     p(X) = \frac{1}{\sqrt{2 \pi \sigma^2 }} \exp\left( -\frac{(x-\mu)^2}{2\sigma^2} \right),                                         *
     *                                                                                                                                      *
     *     is undefined in the Riemann or Lebesgue sense. The rule :math:`\delta(a/b) = |a/b| \sqrt{ (\delta a/a)^2 + (\delta b/b)^2 }`     *
     *     can not be derived from the first principles of probability theory.                                                              *
     *                                                                                                                                      *
     ****************************************************************************************************************************************/

    template <typename T>
    class UncorrelatedDeviation {

    public:

        T value;
        T uncertainty;

        /*
         * check for nan
         */

        bool isnan(){
            return value!=value;
        }

        /*
         *  Unary Operators
         */
        UncorrelatedDeviation operator+() const
        {
            return *this;
        }

        UncorrelatedDeviation operator-() const
        {
            return UncorrelatedDeviation( -value, uncertainty);
        }

        /*
         *  Compound assignment operators
         */
        UncorrelatedDeviation& operator+=(const T& number)
        {
            this->value += number;
            return *this;
        };
        UncorrelatedDeviation& operator+=(const UncorrelatedDeviation& gu1)
        {
            this->value += gu1.value;
            this->uncertainty = std::sqrt( this->uncertainty*this->uncertainty + gu1.uncertainty*gu1.uncertainty );
            return *this;
        };

        UncorrelatedDeviation& operator-=(const T& number)
        {
            this->value -= number;
            return *this;
        };
        UncorrelatedDeviation& operator-=(const UncorrelatedDeviation& gu1)
        {
            this->value -= gu1.value;
            this->uncertainty = std::sqrt( this->uncertainty*this->uncertainty + gu1.uncertainty*gu1.uncertainty );
            return *this;
        };

        UncorrelatedDeviation& operator*=(const T& number)
        {
            this->value *= number;
            this->uncertainty = std::sqrt( this->uncertainty*number*this->uncertainty*number );
            return *this;
        };
        UncorrelatedDeviation& operator*=(const UncorrelatedDeviation& gu1)
        {
            T old_value = this->value;
            this->value *= gu1.value;
            this->uncertainty = std::sqrt( this->uncertainty*gu1.value*this->uncertainty*gu1.value + gu1.uncertainty*old_value*gu1.uncertainty*old_value + this->uncertainty*this->uncertainty*gu1.uncertainty*gu1.uncertainty );
            return *this;
        };

        UncorrelatedDeviation& operator/=(const T& number)
        {
            this->value /= number;
            this->uncertainty = std::sqrt( this->uncertainty*this->uncertainty/(number*number) );
            return *this;
        };

        /*
         *  Binary operators
         */
        friend UncorrelatedDeviation operator+(const UncorrelatedDeviation& uc1, const UncorrelatedDeviation& uc2)
        {
            UncorrelatedDeviation output = uc1; // copy
            output += uc2;
            return output;
        };
        friend UncorrelatedDeviation operator+(const UncorrelatedDeviation& uc1, const T& value)
        {
            UncorrelatedDeviation output = uc1; // copy
            output += value;
            return output;
        };
        friend UncorrelatedDeviation operator+(const T& value, const UncorrelatedDeviation& uc)
        {
            UncorrelatedDeviation output = uc; // copy
            output += value;
            return output;
        };

        friend UncorrelatedDeviation operator-(const UncorrelatedDeviation& uc1, const UncorrelatedDeviation& uc2)
        {
            UncorrelatedDeviation output = uc1; // copy
            output -= uc2;
            return output;
        };
        friend UncorrelatedDeviation operator-(const UncorrelatedDeviation& uc1, const T& value)
        {
            UncorrelatedDeviation output = uc1; // copy
            output -= value;
            return output;
        };
        friend UncorrelatedDeviation operator-(const T& value, const UncorrelatedDeviation& uc)
        {
            UncorrelatedDeviation output = -uc; // copy
            output += value;
            return output;
        };

        friend UncorrelatedDeviation operator*(const UncorrelatedDeviation& uc1, const UncorrelatedDeviation& uc2)
        {
            UncorrelatedDeviation output = uc1; // copy
            output *= uc2;
            return output;
        };
        friend UncorrelatedDeviation operator*(const UncorrelatedDeviation& uc1, const T& value)
        {
            UncorrelatedDeviation output = uc1; // copy
            output *= value;
            return output;
        };
        friend UncorrelatedDeviation operator*(const T& value, const UncorrelatedDeviation& uc2)
        {
            UncorrelatedDeviation output = uc2; // copy
            output *= value;
            return output;
        };

        friend UncorrelatedDeviation operator/(const UncorrelatedDeviation& uc1, const T& value)
        {
            UncorrelatedDeviation output = uc1; // copy
            output /= value;
            return output;
        };

        /*
         *  Printing
         */
        friend std::ostream& operator<< (std::ostream& os, const UncorrelatedDeviation& gu1)
        {
            os << gu1.value << " +/- " << gu1.uncertainty;
            return os;
        };

        /*
         *  Constructors
         */
        UncorrelatedDeviation(T value, T uncertainty):
        value(value), uncertainty(uncertainty)
        {};

        // construct with zero uncertainty
        UncorrelatedDeviation(T value):
        value(value), uncertainty(0)
        {};

        // default constructor
        UncorrelatedDeviation():
        value(0), uncertainty(0)
        {};

    };


    // specialization for complex
    #define UNCORRELATEDDEVIATIONCOMPLEX(complex_template) template <typename T> \
    class UncorrelatedDeviation<complex_template<T>> { \
 \
    public: \
 \
        complex_template<T> value; \
        complex_template<T> uncertainty; \
 \
        /* \
         * check for nan \
         */ \
 \
        bool isnan(){ \
            return value!=value; \
        } \
 \
        /* \
         *  Unary Operators \
         */ \
        UncorrelatedDeviation operator+() const \
        { \
            return *this; \
        } \
 \
        UncorrelatedDeviation operator-() const \
        { \
            return UncorrelatedDeviation( -value, uncertainty); \
        } \
 \
        /* \
         *  Compound assignment operators \
         */ \
        UncorrelatedDeviation& operator+=(const T& number) \
        { \
            this->value += number; \
            return *this; \
        }; \
        UncorrelatedDeviation& operator+=(const complex_template<T>& number) \
        { \
            this->value += number; \
            return *this; \
        }; \
        UncorrelatedDeviation& operator+=(const UncorrelatedDeviation& gu1) \
        { \
            this->value += gu1.value; \
            this->uncertainty = complex_template<T> \
                                {std::sqrt( this->uncertainty.real()*this->uncertainty.real() + gu1.uncertainty.real()*gu1.uncertainty.real() ),  /* real part */ \
                                 std::sqrt( this->uncertainty.imag()*this->uncertainty.imag() + gu1.uncertainty.imag()*gu1.uncertainty.imag() )}; /* imaginary part */ \
            return *this; \
        }; \
 \
        UncorrelatedDeviation& operator-=(const T& number) \
        { \
            this->value -= number; \
            return *this; \
        }; \
        UncorrelatedDeviation& operator-=(const complex_template<T>& number) \
        { \
            this->value -= number; \
            return *this; \
        }; \
        UncorrelatedDeviation& operator-=(const UncorrelatedDeviation& gu1) \
        { \
            this->value -= gu1.value; \
            this->uncertainty = complex_template<T> \
                                {std::sqrt( this->uncertainty.real()*this->uncertainty.real() + gu1.uncertainty.real()*gu1.uncertainty.real() ),  /* real part */ \
                                 std::sqrt( this->uncertainty.imag()*this->uncertainty.imag() + gu1.uncertainty.imag()*gu1.uncertainty.imag() )}; /* imaginary part */ \
            return *this; \
        }; \
 \
        UncorrelatedDeviation& operator*=(const T& number) \
        { \
            auto real0 = UncorrelatedDeviation<T>(this->value.real(), this->uncertainty.real()); \
            auto imag0 = UncorrelatedDeviation<T>(this->value.imag(), this->uncertainty.imag()); \
 \
            auto new_real_part = real0*number; \
            auto new_imag_part = number*imag0; \
 \
            this->value = complex_template<T>{new_real_part.value,new_imag_part.value}; \
            this->uncertainty = complex_template<T>{new_real_part.uncertainty,new_imag_part.uncertainty}; \
 \
            return *this; \
        }; \
        UncorrelatedDeviation& operator*=(const complex_template<T>& number) \
        { \
            auto real0 = UncorrelatedDeviation<T>(this->value.real(), this->uncertainty.real()); \
            auto imag0 = UncorrelatedDeviation<T>(this->value.imag(), this->uncertainty.imag()); \
 \
            auto new_real_part = real0*number.real() - imag0*number.imag(); \
            auto new_imag_part = real0*number.imag() + number.real()*imag0; \
 \
            this->value = complex_template<T>{new_real_part.value,new_imag_part.value}; \
            this->uncertainty = complex_template<T>{new_real_part.uncertainty,new_imag_part.uncertainty}; \
 \
            return *this; \
        }; \
        UncorrelatedDeviation& operator*=(const UncorrelatedDeviation& gu1) \
        { \
            auto real0 = UncorrelatedDeviation<T>(this->value.real(), this->uncertainty.real()); \
            auto imag0 = UncorrelatedDeviation<T>(this->value.imag(), this->uncertainty.imag()); \
 \
            auto real1 = UncorrelatedDeviation<T>(gu1.value.real(), gu1.uncertainty.real()); \
            auto imag1 = UncorrelatedDeviation<T>(gu1.value.imag(), gu1.uncertainty.imag()); \
 \
            auto new_real_part = real0*real1 - imag0*imag1; \
            auto new_imag_part = real0*imag1 + real1*imag0; \
 \
            this->value = complex_template<T>{new_real_part.value,new_imag_part.value}; \
            this->uncertainty = complex_template<T>{new_real_part.uncertainty,new_imag_part.uncertainty}; \
 \
            return *this; \
        }; \
 \
        UncorrelatedDeviation& operator/=(const T& number) \
        { \
            auto real0 = UncorrelatedDeviation<T>(this->value.real(), this->uncertainty.real()); \
            auto imag0 = UncorrelatedDeviation<T>(this->value.imag(), this->uncertainty.imag()); \
 \
            auto denominator = number*number; \
 \
            auto new_real_part = (real0*number) / denominator; \
            auto new_imag_part = (number*imag0) / denominator; \
 \
            this->value = complex_template<T>{new_real_part.value,new_imag_part.value}; \
            this->uncertainty = complex_template<T>{new_real_part.uncertainty,new_imag_part.uncertainty}; \
 \
            return *this; \
        }; \
        UncorrelatedDeviation& operator/=(const complex_template<T>& number) \
        { \
            auto real0 = UncorrelatedDeviation<T>(this->value.real(), this->uncertainty.real()); \
            auto imag0 = UncorrelatedDeviation<T>(this->value.imag(), this->uncertainty.imag()); \
 \
            auto denominator = number.real()*number.real() + number.imag()*number.imag(); \
 \
            auto new_real_part = (real0*number.real() + imag0*number.imag()) / denominator; \
            auto new_imag_part = (number.real()*imag0 - real0*number.imag()) / denominator; \
 \
            this->value = complex_template<T>{new_real_part.value,new_imag_part.value}; \
            this->uncertainty = complex_template<T>{new_real_part.uncertainty,new_imag_part.uncertainty}; \
 \
            return *this; \
        }; \
 \
        /* \
         *  Binary operators \
         */ \
        friend UncorrelatedDeviation operator+(const UncorrelatedDeviation& uc1, const UncorrelatedDeviation& uc2) \
        { \
            UncorrelatedDeviation output = uc1; /* copy */ \
            output += uc2; \
            return output; \
        }; \
        friend UncorrelatedDeviation operator+(const UncorrelatedDeviation& uc1, const T& value) \
        { \
            UncorrelatedDeviation output = uc1; /* copy */ \
            output += value; \
            return output; \
        }; \
        friend UncorrelatedDeviation operator+(const UncorrelatedDeviation& uc1, const complex_template<T>& value) \
        { \
            UncorrelatedDeviation output = uc1; /* copy */ \
            output += value; \
            return output; \
        }; \
        friend UncorrelatedDeviation operator+(const T& value, const UncorrelatedDeviation& uc) \
        { \
            UncorrelatedDeviation output = value; /* convert */ \
            output += uc; \
            return output; \
        }; \
        friend UncorrelatedDeviation operator+(const complex_template<T>& value, const UncorrelatedDeviation& uc) \
        { \
            UncorrelatedDeviation output = value; /* convert */ \
            output += uc; \
            return output; \
        }; \
 \
        friend UncorrelatedDeviation operator-(const UncorrelatedDeviation& uc1, const UncorrelatedDeviation& uc2) \
        { \
            UncorrelatedDeviation output = uc1; /* copy */ \
            output -= uc2; \
            return output; \
        }; \
        friend UncorrelatedDeviation operator-(const UncorrelatedDeviation& uc1, const T& value) \
        { \
            UncorrelatedDeviation output = uc1; /* copy */ \
            output -= value; \
            return output; \
        }; \
        friend UncorrelatedDeviation operator-(const UncorrelatedDeviation& uc1, const complex_template<T>& value) \
        { \
            UncorrelatedDeviation output = uc1; /* copy */ \
            output -= value; \
            return output; \
        }; \
        friend UncorrelatedDeviation operator-(const T& value, const UncorrelatedDeviation& uc) \
        { \
            UncorrelatedDeviation output = value; /* convert */ \
            output -= uc; \
            return output; \
        }; \
        friend UncorrelatedDeviation operator-(const complex_template<T>& value, const UncorrelatedDeviation& uc) \
        { \
            UncorrelatedDeviation output = value; /* convert */ \
            output -= uc; \
            return output; \
        }; \
 \
        friend UncorrelatedDeviation operator*(const UncorrelatedDeviation& uc1, const UncorrelatedDeviation& uc2) \
        { \
            UncorrelatedDeviation output = uc1; /* copy */ \
            output *= uc2; \
            return output; \
        }; \
        friend UncorrelatedDeviation operator*(const UncorrelatedDeviation& uc1, const T& value) \
        { \
            UncorrelatedDeviation output = uc1; /* copy*/ \
            output *= value; \
            return output; \
        }; \
        friend UncorrelatedDeviation operator*(const UncorrelatedDeviation& uc1, const complex_template<T>& value) \
        { \
            UncorrelatedDeviation output = uc1; /* copy */ \
            output *= value; \
            return output; \
        }; \
        friend UncorrelatedDeviation operator*(const T& value, const UncorrelatedDeviation& uc2) \
        { \
            UncorrelatedDeviation output = value; /* convert */ \
            output *= uc2; \
            return output; \
        }; \
        friend UncorrelatedDeviation operator*(const complex_template<T>& value, const UncorrelatedDeviation& uc2) \
        { \
            UncorrelatedDeviation output = value; /* convert */ \
            output *= uc2; \
            return output; \
        }; \
 \
        friend UncorrelatedDeviation operator/(const UncorrelatedDeviation& uc1, const T& value) \
        { \
            UncorrelatedDeviation output = uc1; /* copy */ \
            output /= value; \
            return output; \
        }; \
        friend UncorrelatedDeviation operator/(const UncorrelatedDeviation& uc1, const complex_template<T>& value) \
        { \
            UncorrelatedDeviation output = uc1; /* copy */ \
            output /= value; \
            return output; \
        }; \
 \
        /* \
         *  Printing \
         */ \
        friend std::ostream& operator<< (std::ostream& os, const UncorrelatedDeviation& gu1) \
        { \
            os << gu1.value << " +/- " << gu1.uncertainty; \
            return os; \
        }; \
 \
        /* \
         *  Constructors \
         */ \
        UncorrelatedDeviation(complex_template<T> value, complex_template<T> uncertainty): \
        value(value), uncertainty(uncertainty) \
        {}; \
 \
        /* construct with zero uncertainty */ \
        UncorrelatedDeviation(T value): \
        value(value), uncertainty(0) \
        {}; \
        UncorrelatedDeviation(complex_template<T> value): \
        value(value), uncertainty(0) \
        {}; \
        UncorrelatedDeviation(): \
        value(0), uncertainty(0) \
        {}; \
 \
        /* converting constructor "T -> complex<T>" */ \
        UncorrelatedDeviation(UncorrelatedDeviation<T> ud): \
        value(ud.value), uncertainty(ud.uncertainty) \
        {}; \
 \
    };

    UNCORRELATEDDEVIATIONCOMPLEX(std::complex)
    #ifdef SECDEC_WITH_CUDA
        UNCORRELATEDDEVIATIONCOMPLEX(thrust::complex)
    #endif

    #undef UNCORRELATEDDEVIATIONCOMPLEX

}

#endif

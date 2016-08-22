#ifndef SecDecUtil_uncertainties_hpp_included
#define SecDecUtil_uncertainties_hpp_included

#include <cmath>
#include <complex>
#include <ostream>

namespace secdecutil {

    /*!
     * Implement the propagation of !!UNCORRELATED!!
     * standard deviations via the class
     * "UncorrelatedDeviation".
     *
     */

    template <typename T>
    class UncorrelatedDeviation {

    public:

        T value;
        T uncertainty;

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
        template <typename Tp>
        UncorrelatedDeviation<Tp>& operator+=(const UncorrelatedDeviation<Tp>& gu1)
        {
            this->value += gu1.value;
            this->uncertainty = std::sqrt( this->uncertainty*this->uncertainty + gu1.uncertainty*gu1.uncertainty );
            return *this;
        };
        template<typename Tinner>
        UncorrelatedDeviation<std::complex<Tinner>>& operator+=(const UncorrelatedDeviation<std::complex<Tinner>>& gu1)
        {
            this->value += gu1.value;
            this->uncertainty = {std::sqrt( this->uncertainty.real()*this->uncertainty.real() + gu1.uncertainty.real()*gu1.uncertainty.real() ),  // real part
                                 std::sqrt( this->uncertainty.imag()*this->uncertainty.imag() + gu1.uncertainty.imag()*gu1.uncertainty.imag() )}; // imaginary part
            return *this;
        };

        template <typename Tp>
        UncorrelatedDeviation<Tp>& operator-=(const UncorrelatedDeviation<Tp>& gu1)
        {
            this->value -= gu1.value;
            this->uncertainty = std::sqrt( this->uncertainty*this->uncertainty + gu1.uncertainty*gu1.uncertainty );
            return *this;
        };
        template<typename Tinner>
        UncorrelatedDeviation<std::complex<Tinner>>& operator-=(const UncorrelatedDeviation<std::complex<Tinner>>& gu1)
        {
            this->value -= gu1.value;
            this->uncertainty = {std::sqrt( this->uncertainty.real()*this->uncertainty.real() + gu1.uncertainty.real()*gu1.uncertainty.real() ),  // real part
                                 std::sqrt( this->uncertainty.imag()*this->uncertainty.imag() + gu1.uncertainty.imag()*gu1.uncertainty.imag() )}; // imaginary part
            return *this;
        };

        template <typename Tp>
        UncorrelatedDeviation<Tp>& operator*=(const UncorrelatedDeviation<Tp>& gu1)
        {
            T old_value = this->value;
            this->value *= gu1.value;
            this->uncertainty = std::sqrt( this->uncertainty*gu1.value*this->uncertainty*gu1.value + gu1.uncertainty*old_value*gu1.uncertainty*old_value );
            return *this;
        };
        template<typename Tinner>
        UncorrelatedDeviation<std::complex<Tinner>>& operator*=(const UncorrelatedDeviation<std::complex<Tinner>>& gu1)
        {
            auto real0 = UncorrelatedDeviation<Tinner>(this->value.real(), this->uncertainty.real());
            auto imag0 = UncorrelatedDeviation<Tinner>(this->value.imag(), this->uncertainty.imag());

            auto real1 = UncorrelatedDeviation<Tinner>(gu1.value.real(), gu1.uncertainty.real());
            auto imag1 = UncorrelatedDeviation<Tinner>(gu1.value.imag(), gu1.uncertainty.imag());

            auto new_real_part = real0*real1 - imag0*imag1;
            auto new_imag_part = real0*imag1 + real1*imag0;

            this->value = {new_real_part.value,new_imag_part.value};
            this->uncertainty = {new_real_part.uncertainty,new_imag_part.uncertainty};

            return *this;
        };

        template<typename Tp>
        UncorrelatedDeviation<Tp>& operator/=(const UncorrelatedDeviation<Tp>& gu1)
        {
            T old_value = this->value;
            this->value /= gu1.value;
            this->uncertainty = std::sqrt( this->uncertainty*this->uncertainty/(gu1.value*gu1.value) + this->value*this->value*gu1.uncertainty*gu1.uncertainty/(gu1.value*gu1.value) );
            return *this;
        };
        template<typename Tinner>
        UncorrelatedDeviation<std::complex<Tinner>>& operator/=(const UncorrelatedDeviation<std::complex<Tinner>>& gu1)
        {
            auto real0 = UncorrelatedDeviation<Tinner>(this->value.real(), this->uncertainty.real());
            auto imag0 = UncorrelatedDeviation<Tinner>(this->value.imag(), this->uncertainty.imag());

            auto real1 = UncorrelatedDeviation<Tinner>(gu1.value.real(), gu1.uncertainty.real());
            auto imag1 = UncorrelatedDeviation<Tinner>(gu1.value.imag(), gu1.uncertainty.imag());

            auto denominator = real1*real1 + imag1*imag1;

            auto new_real_part = (real0*real1 + imag0*imag1) / denominator;
            auto new_imag_part = (real1*imag0 - real0*imag1) / denominator;

            this->value = {new_real_part.value,new_imag_part.value};
            this->uncertainty = {new_real_part.uncertainty,new_imag_part.uncertainty};

            return *this;
        };

        /*
         *  Binary operators
         */
        friend UncorrelatedDeviation operator+(const UncorrelatedDeviation& gu1, const UncorrelatedDeviation& gu2)
        {
            auto output = gu1; // copy
            output += gu2;
            return output;
        };

        friend UncorrelatedDeviation operator-(const UncorrelatedDeviation& gu1, const UncorrelatedDeviation& gu2)
        {
            auto output = gu1; // copy
            output -= gu2;
            return output;
        };

        friend UncorrelatedDeviation operator*(const UncorrelatedDeviation& gu1, const UncorrelatedDeviation& gu2)
        {
            auto output = gu1; // copy
            output *= gu2;
            return output;
        };

        friend UncorrelatedDeviation operator/(const UncorrelatedDeviation& gu1, const UncorrelatedDeviation& gu2)
        {
            auto output = gu1; // copy
            output /= gu2;
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
         *  Constructor
         */
        UncorrelatedDeviation(T value, T uncertainty):
        value(value), uncertainty(uncertainty)
        {};

        // construct with zero uncertainty
        UncorrelatedDeviation(T value):
        value(value), uncertainty(0)
        {};

    };

}

#endif

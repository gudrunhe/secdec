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
        UncorrelatedDeviation& operator+=(const Tp& number)
        {
            this->value += number;
            return *this;
        };
        template <typename Tp>
        UncorrelatedDeviation& operator+=(const UncorrelatedDeviation<Tp>& gu1)
        {
            this->value += gu1.value;
            this->uncertainty = std::sqrt( this->uncertainty*this->uncertainty + gu1.uncertainty*gu1.uncertainty );
            return *this;
        };
        template<typename Tinner>
        UncorrelatedDeviation& operator+=(const UncorrelatedDeviation<std::complex<Tinner>>& gu1)
        {
            this->value += gu1.value;
            this->uncertainty = {std::sqrt( this->uncertainty.real()*this->uncertainty.real() + gu1.uncertainty.real()*gu1.uncertainty.real() ),  // real part
                                 std::sqrt( this->uncertainty.imag()*this->uncertainty.imag() + gu1.uncertainty.imag()*gu1.uncertainty.imag() )}; // imaginary part
            return *this;
        };

        template <typename Tp>
        UncorrelatedDeviation& operator-=(const Tp& number)
        {
            this->value -= number;
            return *this;
        };
        template <typename Tp>
        UncorrelatedDeviation& operator-=(const UncorrelatedDeviation<Tp>& gu1)
        {
            this->value -= gu1.value;
            this->uncertainty = std::sqrt( this->uncertainty*this->uncertainty + gu1.uncertainty*gu1.uncertainty );
            return *this;
        };
        template<typename Tinner>
        UncorrelatedDeviation& operator-=(const UncorrelatedDeviation<std::complex<Tinner>>& gu1)
        {
            this->value -= gu1.value;
            this->uncertainty = {std::sqrt( this->uncertainty.real()*this->uncertainty.real() + gu1.uncertainty.real()*gu1.uncertainty.real() ),  // real part
                                 std::sqrt( this->uncertainty.imag()*this->uncertainty.imag() + gu1.uncertainty.imag()*gu1.uncertainty.imag() )}; // imaginary part
            return *this;
        };

        template <typename Tp>
        UncorrelatedDeviation& operator*=(const Tp& number)
        {
            this->value *= number;
            this->uncertainty = std::sqrt( this->uncertainty*number*this->uncertainty*number );
            return *this;
        };
        template<typename Tinner>
        UncorrelatedDeviation& operator*=(const std::complex<Tinner>& number)
        {
            auto real0 = UncorrelatedDeviation<Tinner>(this->value.real(), this->uncertainty.real());
            auto imag0 = UncorrelatedDeviation<Tinner>(this->value.imag(), this->uncertainty.imag());

            auto new_real_part = real0*number.real() - imag0*number.imag();
            auto new_imag_part = real0*number.imag() + number.real()*imag0;

            this->value = {new_real_part.value,new_imag_part.value};
            this->uncertainty = {new_real_part.uncertainty,new_imag_part.uncertainty};

            return *this;
        };
        template <typename Tp>
        UncorrelatedDeviation& operator*=(const UncorrelatedDeviation<Tp>& gu1)
        {
            T old_value = this->value;
            this->value *= gu1.value;
            this->uncertainty = std::sqrt( this->uncertainty*gu1.value*this->uncertainty*gu1.value + gu1.uncertainty*old_value*gu1.uncertainty*old_value );
            return *this;
        };
        template<typename Tinner>
        UncorrelatedDeviation& operator*=(const UncorrelatedDeviation<std::complex<Tinner>>& gu1)
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
        UncorrelatedDeviation& operator/=(const Tp& number)
        {
            this->value /= number;
            this->uncertainty = std::sqrt( this->uncertainty*this->uncertainty/(number*number) );
            return *this;
        };
        template<typename Tinner>
        UncorrelatedDeviation& operator/=(const std::complex<Tinner>& number)
        {
            auto real0 = UncorrelatedDeviation<Tinner>(this->value.real(), this->uncertainty.real());
            auto imag0 = UncorrelatedDeviation<Tinner>(this->value.imag(), this->uncertainty.imag());

            auto denominator = number.real()*number.real() + number.imag()*number.imag();

            auto new_real_part = (real0*number.real() + imag0*number.imag()) / denominator;
            auto new_imag_part = (number.real()*imag0 - real0*number.imag()) / denominator;

            this->value = {new_real_part.value,new_imag_part.value};
            this->uncertainty = {new_real_part.uncertainty,new_imag_part.uncertainty};

            return *this;
        };
        template<typename Tp>
        UncorrelatedDeviation& operator/=(const UncorrelatedDeviation<Tp>& gu1)
        {
            this->value /= gu1.value;
            this->uncertainty = std::sqrt( this->uncertainty*this->uncertainty/(gu1.value*gu1.value) + this->value*this->value*gu1.uncertainty*gu1.uncertainty/(gu1.value*gu1.value) );
            return *this;
        };
        template<typename Tinner>
        UncorrelatedDeviation& operator/=(const UncorrelatedDeviation<std::complex<Tinner>>& gu1)
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
        template<typename Tin>
        UncorrelatedDeviation(Tin value):
        value(value), uncertainty(0)
        {};

        // converting constructor
        template<typename Tother>
        UncorrelatedDeviation(UncorrelatedDeviation<Tother> ud):
        value(ud.value), uncertainty(ud.uncertainty)
        {};

    };

    /*
     *  Binary operators
     */
    template<typename T1, typename T2>
    inline auto operator+(const UncorrelatedDeviation<T1>& uc1, const UncorrelatedDeviation<T2>& uc2)
    -> UncorrelatedDeviation<decltype(uc1.value + uc2.value)>
    {
        UncorrelatedDeviation<decltype(uc1.value + uc2.value)> output = uc1; // copy
        output += uc2;
        return output;
    };
    template<typename T1, typename T2>
    inline auto operator+(const UncorrelatedDeviation<T1>& uc1, const T2& value)
    -> UncorrelatedDeviation<decltype(uc1.value + value)>
    {
        UncorrelatedDeviation<decltype(uc1.value + value)> output = uc1; // copy
        output += value;
        return output;
    };
    template<typename T1, typename T2>
    inline auto operator+(const T1& value, const UncorrelatedDeviation<T2>& uc)
    -> UncorrelatedDeviation<decltype(value + uc.value)>
    {
        UncorrelatedDeviation<decltype(value + uc.value)> output = uc; // copy
        output += value;
        return output;
    };

    template<typename T1, typename T2>
    inline auto operator-(const UncorrelatedDeviation<T1>& uc1, const UncorrelatedDeviation<T2>& uc2)
    -> UncorrelatedDeviation<decltype(uc1.value - uc2.value)>
    {
        UncorrelatedDeviation<decltype(uc1.value - uc2.value)> output = uc1; // copy
        output -= uc2;
        return output;
    };
    template<typename T1, typename T2>
    inline auto operator-(const UncorrelatedDeviation<T1>& uc1, const T2& value)
    -> UncorrelatedDeviation<decltype(uc1.value - value)>
    {
        UncorrelatedDeviation<decltype(uc1.value - value)> output = uc1; // copy
        output -= value;
        return output;
    };
    template<typename T1, typename T2>
    inline auto operator-(const T1& value, const UncorrelatedDeviation<T2>& uc)
    -> UncorrelatedDeviation<decltype(value + uc.value)>
    {
        UncorrelatedDeviation<decltype(value - uc.value)> output = -uc; // copy
        output += value;
        return output;
    };

    template<typename T1, typename T2>
    inline auto operator*(const UncorrelatedDeviation<T1>& uc1, const UncorrelatedDeviation<T2>& uc2)
    -> UncorrelatedDeviation<decltype(uc1.value * uc2.value)>
    {
        UncorrelatedDeviation<decltype(uc1.value * uc2.value)> output = uc1; // copy
        output *= uc2;
        return output;
    };
    template<typename T1, typename T2>
    inline auto operator*(const UncorrelatedDeviation<T1>& uc1, const T2& value)
    -> UncorrelatedDeviation<decltype(uc1.value * value)>
    {
        UncorrelatedDeviation<decltype(uc1.value * value)> output = uc1; // copy
        output *= value;
        return output;
    };
    template<typename T1, typename T2>
    inline auto operator*(const T1& value, const UncorrelatedDeviation<T2>& uc2)
    -> UncorrelatedDeviation<decltype(value * uc2.value)>
    {
        UncorrelatedDeviation<decltype(value * uc2.value)> output = uc2; // copy
        output *= value;
        return output;
    };

    template<typename T1, typename T2>
    inline auto operator/(const UncorrelatedDeviation<T1>& uc1, const UncorrelatedDeviation<T2>& uc2)
    -> UncorrelatedDeviation<decltype(uc1.value / uc2.value)>
    {
        UncorrelatedDeviation<decltype(uc1.value / uc2.value)> output = uc1; // copy
        output /= uc2;
        return output;
    };
    template<typename T1, typename T2>
    inline auto operator/(const UncorrelatedDeviation<T1>& uc1, const T2& value)
    -> UncorrelatedDeviation<decltype(uc1.value / value)>
    {
        UncorrelatedDeviation<decltype(uc1.value / value)> output = uc1; // copy
        output /= value;
        return output;
    };
    template<typename T1, typename T2>
    inline auto operator/(const T1& value, const UncorrelatedDeviation<T2>& uc)
    -> UncorrelatedDeviation<decltype(value / uc.value)>
    {
        UncorrelatedDeviation<decltype(value / uc.value)> output = value; // copy
        output /= uc;
        return output;
    };

}

#endif

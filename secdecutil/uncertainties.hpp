#ifndef SecDecUtil_uncertainties_hpp_included
#define SecDecUtil_uncertainties_hpp_included

#include <cmath>

namespace secdecutil {

    // Warning - assumes uncorrelated errors, neglects covariance matrix
    // TODO specialisations for complex numbers

    template <typename T>
    class GaussianUncertainty {

    private:
        enum Operation { add, subtract, multiply, divide };

    public:
        
        T value;
        T uncertainty;

        /*
         *  Helper functions
         */
        template<int operation>
        static GaussianUncertainty add_subtract_multiply_or_divide(const GaussianUncertainty& gu1, const GaussianUncertainty& gu2)
        {
            T value;
            T uncertainty;

            if (operation == add) {

                value = gu1.value + gu2.value;
                uncertainty = sqrt( pow( gu1.uncertainty, 2) + pow( gu2.uncertainty, 2) );

            } else if (operation == subtract ) {

                value = gu1.value - gu2.value;
                uncertainty = sqrt( pow( gu1.uncertainty, 2) + pow( gu2.uncertainty, 2) );

            } else if (operation == multiply ) {

                value = gu1.value * gu2.value;
                uncertainty = fabs(value) * sqrt( pow(gu1.uncertainty/gu1.value, 2) + pow(gu2.uncertainty/gu2.value, 2) );

            } else if ( operation == divide ) {

                value = gu1.value / gu2.value;
                uncertainty = fabs(value) * sqrt( pow(gu1.uncertainty/gu1.value, 2) + pow(gu2.uncertainty/gu2.value, 2) );

            }

            return GaussianUncertainty(value, uncertainty);
        };

        /*
         *  Unary Operators
         */
        GaussianUncertainty operator+() const
        {
            return *this;
        }

        GaussianUncertainty operator-() const
        {
            return GaussianUncertainty( -value, uncertainty);
        }

        /*
         *  Compound assignment operators
         */
        GaussianUncertainty& operator+=(const GaussianUncertainty& gu1)
        {
            *this = *this + gu1;
            return *this;
        };

        GaussianUncertainty& operator-=(const GaussianUncertainty& gu1)
        {
            *this = *this - gu1;
            return *this;
        };

        GaussianUncertainty& operator*=(const GaussianUncertainty& gu1)
        {
            *this = *this * gu1;
            return *this;
        };
        GaussianUncertainty& operator/=(const GaussianUncertainty& gu1)
        {
            *this = *this / gu1;
            return *this;
        };

        /*
         *  Binary operators
         */
        friend GaussianUncertainty operator+(const GaussianUncertainty& gu1, const GaussianUncertainty& gu2)
        {
            return add_subtract_multiply_or_divide<add>(gu1,gu2);
        };

        friend GaussianUncertainty operator-(const GaussianUncertainty& gu1, const GaussianUncertainty& gu2)
        {
            return add_subtract_multiply_or_divide<subtract>(gu1,gu2);
        };

        friend GaussianUncertainty operator*(const GaussianUncertainty& gu1, const GaussianUncertainty& gu2)
        {
            return add_subtract_multiply_or_divide<multiply>(gu1,gu2);
        };

        friend GaussianUncertainty operator/(const GaussianUncertainty& gu1, const GaussianUncertainty& gu2)
        {
            return add_subtract_multiply_or_divide<divide>(gu1,gu2);
        };

        /*
         *  Printing
         */
        friend std::ostream& operator<< (std::ostream& os, const GaussianUncertainty& gu1)
        {
            os << "( " << gu1.value << " +/- " << gu1.uncertainty << " )";
            return os;
        };

        /*
         *  Constructor
         */
        GaussianUncertainty(T value, T uncertainty):
        value (value), uncertainty(uncertainty)
        {};

    };
    
}

#endif
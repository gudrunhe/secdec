#ifndef __SecDec_include_guard_series
#define __SecDec_include_guard_series

#include <vector>
#include <stdexcept>
#include <iostream> // std::ostream

namespace secdecutil {

    // Store \Sum_{i=order_min}^order_max c_i * eps^i
    template <typename T>
    class Series {

    protected:

        int order_min;
        int order_max;
        std::vector<T> content;
        bool truncated_above;

        /*
         *  Helper functions
         */

        // Checks if series has a term at order order
        bool hasTerm(int order) const
        {
            if ( (order >= order_min) && (order <= order_max) )
                return true;

            return false;
        }
        // Performs subtraction U=true or addition U=false of two series
        // \Sum_i=a1^b1 c_i * eps^i - \Sum_j=a2^b2 c'_j * eps^j = \Sum_i=min(a1,a2)^max(b1,b2) (c_i - c'_i) * eps^i
        // or
        // \Sum_i=a1^b1 c_i * eps^i + \Sum_j=a2^b2 c'_j * eps^j = \Sum_i=min(a1,a2)^max(b1,b2) (c_i + c'_i) * eps^i
        template<bool U>
        static Series sub_or_add(const Series& s1, const Series& s2)
        {
            // TODO Check expansion parameter

            // Assume series are not truncated
            bool truncated_above = false;
            int order_min = std::min( s1.order_min, s2.order_min );
            int order_max = std::max( s1.order_max, s2.order_max );

            // Check if new series must be truncated from above
            if ( s1.truncated_above  )
            {
                truncated_above = true;
                if (s1.order_max < order_max)
                    order_max = s1.order_max;
            }

            if ( s2.truncated_above )
            {
                truncated_above = true;
                if (s2.order_max < order_max)
                    order_max = s2.order_max;
            }

            std::vector<T> content;
            content.reserve(order_max-order_min+1);
            if ( U ) {
                // Perform subtraction
                for ( int i = order_min; i < order_max + 1; i++)
                {
                    if ( s1.hasTerm(i) && s2.hasTerm(i) )
                    {
                        content.push_back(s1.at(i) - s2.at(i));
                    } else if ( s1.hasTerm(i) )
                    {
                        content.push_back(s1.at(i));
                    } else // s2.hasTerm(i)
                    {
                        content.push_back(-s2.at(i));
                    }
                }
            } else {
                // Perform addition
                for ( int i = order_min; i < order_max + 1; i++)
                {
                    if ( s1.hasTerm(i) && s2.hasTerm(i) )
                    {
                        content.push_back(s1.at(i) + s2.at(i));
                    } else if ( s1.hasTerm(i) )
                    {
                        content.push_back(s1.at(i));
                    } else // s2.hasTerm(i)
                    {
                        content.push_back(s2.at(i));
                    }
                }
            }
            return Series(order_min, order_max, content, truncated_above);
        }

    public:

        const int get_order_min() { return order_min; }
        const int get_order_max() { return order_max; }
        const bool get_truncated_above() { return truncated_above; }

        const T& operator[](const int i) const { return content[i-order_min]; }
        T& operator[](const int i)             { return content[i-order_min]; }
        const T& at(const int i) const         { return content.at(i-order_min); }
        T& at(const int i)                     { return content.at(i-order_min); }

        /*
         *  Comparator Operators
         */
        friend bool operator==(const Series& s1, const Series& s2)
        {
            // TODO Check expansion parameter

            if ( s1.order_min != s2.order_min )
                return false;
            if ( s1.order_max != s2.order_max )
                return false;
            if ( s1.truncated_above != s2.truncated_above )
                return false;
            if ( s1.content != s2.content )
                return false;
            return true;
        }

        friend bool operator!= (const Series& s1, const Series& s2)
        {
            return !( s1 == s2 );
        }

        /*
         *  Unary Operators
         */
        Series operator-() const
        {
            return *this * (-1);
        }

        Series operator+() const
        {
            return *this;
        }

        /*
         *  Compound assignment operators
         */
        Series& operator-=(const Series& s1)
        {
            *this = *this - s1;
            return *this;
        }

        Series& operator+=(const Series& s1)
        {
            *this = *this + s1;
            return *this;
        }

        Series& operator*=(const Series& s1)
        {
            *this = *this * s1;
            return *this;
        }

        Series& operator*=(const T& val)
        {
            *this = *this * val;
            return *this;
        }

        // TODO operator/=

        /*
         *  Binary operators
         */
        friend Series operator-(const Series& s1, const Series& s2)
        {
            return sub_or_add<true>(s1,s2);
        }

        friend Series operator+(const Series& s1, const Series& s2)
        {
            return sub_or_add<false>(s1,s2);
        };

        // (\Sum_i=a1^b1 c_i * eps^i) * (\Sum_j=a2^b2 c'_j * eps^j)
        friend Series operator*(const Series& s1, const Series& s2)
        {
            // TODO Check expansion parameter

            // Assume series are not truncated
            bool truncated_above = false;
            int order_min = s1.order_min + s2.order_min;
            int order_max = s1.order_max + s2.order_max;
            size_t current_index;

            if ( s1.truncated_above )
            {
                truncated_above = true;
                order_max = std::min( order_max , s1.order_max + s2.order_min );
            }
            if ( s2.truncated_above )
            {
                truncated_above = true;
                order_max = std::min( order_max , s1.order_min + s2.order_max );
            }

            // Perform multiplication
            std::vector<T> content;
            content.reserve(order_max-order_min+1);
            for ( int i = s1.order_min; i < s1.order_max + 1; i++ )
            {
                for ( int j = s2.order_min; j < s2.order_max + 1; j++ )
                {
                    current_index = i+j-order_min;
                    if ( ( (i+j) >= order_min ) && ( (i+j) <= order_max ) && s1.hasTerm(i) && s2.hasTerm(j) )
                        if ( current_index < content.size() ) // term exists
                            content.at(current_index) += ( s1.at(i) * s2.at(j) );
                        else // term must be created
                            content.push_back( s1.at(i) * s2.at(j) );
                }
            }

            return Series(order_min, order_max, content, truncated_above);
        }

        // d * \Sum_i=a1^b1 c_i * eps^i = \Sum_i=a1^b1 d* c_i * eps^i
        friend Series operator*( const Series& s1, const T& val)
        {
            std::vector<T> content;
            content.reserve(s1.order_max-s1.order_min+1);
            for ( int i = s1.order_min; i < s1.order_max + 1; i++ )
                content.push_back(val * s1.at(i));

            return Series(s1.order_min, s1.order_max, content, s1.truncated_above);
        }

        friend Series operator*( const T& val, const Series& s1)
        {
            return s1*val;
        }

        // TODO friend operator/

        friend std::ostream& operator<< (std::ostream& os, const Series& s1)
        {
            for ( int i = s1.order_min; i < s1.order_max + 1; i++)
                os << " + (" << s1.at(i) << ")*ep^(" << i << ")";

            if ( s1.truncated_above )
                os << " + O(ep^(" << (s1.order_max+1) << "))";

            return os;
        };

        // Constructor
        Series(
               int order_min,
               int order_max,
               std::vector<T> content,
               bool truncated_above = true
               ) :
        order_min(order_min), order_max(order_max),
        content(content), truncated_above(truncated_above)
        {
            if ( content.size() != (order_max-order_min+1) )
                throw std::invalid_argument("Incorrect number of series coefficients. got: " + std::to_string(content.size()) + ", expected: " + std::to_string(order_max-order_min+1));
        }

    };

}

#endif

#ifndef SecDecUtil_series_hpp_included
#define SecDecUtil_series_hpp_included

#include <iostream>
#include <string>
#include <vector>

namespace secdecutil {

    /*!
     * Store an expression of the form
     * "\Sum_{i=order_min}^order_max c_i * eps^i"
     *
     * The Series class overlads arithmetic operators.
     *
     * The expansion parameter can be set by changing
     * the public member variable ``expansion_parameter``.
     */
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

        // Performs subtraction (``subtract=true``) or addition (``subtract=false``) of two series
        // \Sum_i=a1^b1 c_i * eps^i - \Sum_j=a2^b2 c'_j * eps^j = \Sum_i=min(a1,a2)^max(b1,b2) (c_i - c'_i) * eps^i
        // or
        // \Sum_i=a1^b1 c_i * eps^i + \Sum_j=a2^b2 c'_j * eps^j = \Sum_i=min(a1,a2)^max(b1,b2) (c_i + c'_i) * eps^i
        template<bool subtract>
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
            if ( subtract ) {
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

        typedef typename std::vector<T>::iterator iterator;
        typedef typename std::vector<T>::const_iterator const_iterator;
        typedef typename std::vector<T>::reverse_iterator reverse_iterator;
        typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;
        typedef T& reference;
        typedef const T& const_reference;
        typedef T value_type;

        std::string expansion_parameter; // default value "x" set in constructor
        const int get_order_min() const { return order_min; }
        const int get_order_max() const { return order_max; }
        const bool get_truncated_above() const { return truncated_above; }

        iterator begin() noexcept { return content.begin(); }
        const_iterator begin() const noexcept { return content.begin(); }
        iterator end() noexcept { return content.end(); }
        const_iterator end() const noexcept { return content.end(); }
        reverse_iterator rbegin() noexcept { return content.rbegin(); }
        const_reverse_iterator rbegin() const noexcept  { return content.rbegin(); }
        reverse_iterator rend() noexcept { return content.rend(); }
        const_reverse_iterator rend() const noexcept  { return content.rend(); }
        const_iterator cbegin() const noexcept { return content.cbegin(); }
        const_iterator cend() const noexcept { return content.cend(); }
        const_reverse_iterator crbegin() const noexcept { return content.crbegin(); }
        const_reverse_iterator crend() const noexcept { return content.crend(); }

        reference operator[](const int i)             { return content[i-order_min]; }
        const_reference operator[](const int i) const { return content[i-order_min]; }
        reference at(const int i)                     { return content.at(i-order_min); }
        const_reference at(const int i) const         { return content.at(i-order_min); }
        reference front() { return content.front(); }
        const_reference front() const { return content.front(); }
        reference back() { return content.back(); }
        const_reference back() const { return content.back(); }
        value_type* data() noexcept { return content.data(); }

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
                    {
                        if ( current_index < content.size() ) // term exists
                            content.at(current_index) += ( s1.at(i) * s2.at(j) );
                        else // term must be created
                            content.push_back( s1.at(i) * s2.at(j) );
                    }
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
            int i;

            for ( i = s1.order_min; i < s1.order_max + 1; i++)
            {
                os << " + (" << s1.at(i) << ")";
                if ( i != 0 )
                {
                    os << "*" << s1.expansion_parameter;
                    if (i != 1)
                        os << "^" << i;
                }
            }

            if ( s1.truncated_above )
            {
                os << " + O(";
                os << s1.expansion_parameter;
                if ( i != 1 )
                    os << "^" << (i);
                os << ")";
            }

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
        content(content), truncated_above(truncated_above),
        expansion_parameter("x")
        {
            if ( content.size() != (order_max-order_min+1) )
                throw std::invalid_argument("Incorrect number of series coefficients. got: " + std::to_string(content.size()) + ", expected: " + std::to_string(order_max-order_min+1));
        }

    };

}

#endif

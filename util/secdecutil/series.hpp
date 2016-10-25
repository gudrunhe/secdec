#ifndef SecDecUtil_series_hpp_included
#define SecDecUtil_series_hpp_included

#include <exception>
#include <iostream>
#include <string>
#include <vector>

namespace secdecutil {

    // this exception is thrown when a binary operator on Series with different expansion parameters is called
    struct expansion_parameter_mismatch_error : public std::runtime_error { using std::runtime_error::runtime_error; };

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

        // Performs subtraction (``subtract=true``) or addition (``subtract=false``) of two series
        // \Sum_i=a1^b1 c_i * eps^i - \Sum_j=a2^b2 c'_j * eps^j = \Sum_i=min(a1,a2)^max(b1,b2) (c_i - c'_i) * eps^i
        // or
        // \Sum_i=a1^b1 c_i * eps^i + \Sum_j=a2^b2 c'_j * eps^j = \Sum_i=min(a1,a2)^max(b1,b2) (c_i + c'_i) * eps^i
        template<bool subtract, typename T1, typename T2>
        static auto sub_or_add(const Series<T1>& s1, const Series<T2>& s2)
        -> Series<decltype(s1.at(s1.get_order_min()) + s2.at(s2.get_order_min()))>
        {
            using Tout = decltype(s1.at(s1.get_order_min()) + s2.at(s2.get_order_min()));

            if (s1.expansion_parameter != s2.expansion_parameter)
                throw expansion_parameter_mismatch_error("\"" + s1.expansion_parameter + "\" != \"" + s2.expansion_parameter + "\"");

            // Assume series are not truncated
            bool truncated_above = false;
            int order_min = std::min( s1.get_order_min(), s2.get_order_min() );
            int order_max = std::max( s1.get_order_max(), s2.get_order_max() );

            // Check if new series must be truncated from above
            if ( s1.get_truncated_above() )
            {
                truncated_above = true;
                if (s1.get_order_max() < order_max)
                    order_max = s1.get_order_max();
            }

            if ( s2.get_truncated_above() )
            {
                truncated_above = true;
                if (s2.get_order_max() < order_max)
                    order_max = s2.get_order_max();
            }

            std::vector<Tout> content;
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
            return Series<Tout>(order_min, order_max, content, truncated_above, s1.expansion_parameter /* equality of the expansion parameter was checked earlier */);
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
        int get_order_min() const { return order_min; }
        int get_order_max() const { return order_max; }
        bool get_truncated_above() const { return truncated_above; }

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
         *  Helper functions
         */
        // Checks if series has a term at order order
        bool hasTerm(int order) const
        {
            if ( (order >= order_min) && (order <= order_max) )
                return true;

            return false;
        }

        /*
         *  Comparator Operators
         */
        template<typename Tother>
        friend bool operator==(const Series& s1, const Series<Tother>& s2)
        {
            if (s1.expansion_parameter != s2.expansion_parameter)
                return false;
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

        template<typename Tother>
        friend bool operator!= (const Series& s1, const Series<Tother>& s2)
        {
            return !( s1 == s2 );
        }

        /*
         *  Unary Operators
         */
        Series operator-() const
        {
            std::vector<T> content;
            content.reserve(this->order_max-this->order_min+1);
            for ( int i = this->order_min; i < this->order_max + 1; i++ )
                content.push_back(-this->at(i));

            return Series(this->order_min, this->order_max, content, this->truncated_above, this->expansion_parameter);
        }

        Series operator+() const
        {
            return *this;
        }

        /*
         *  Compound assignment operators
         */
        template<typename Tother>
        Series& operator-=(const Tother& other)
        {
            *this = *this - other;
            return *this;
        }

        template<typename Tother>
        Series& operator+=(const Tother& other)
        {
            *this = *this + other;
            return *this;
        }

        template<typename Tother>
        Series& operator*=(const Tother& other)
        {
            *this = *this * other;
            return *this;
        }

        // TODO operator/=

        /*
         *  Binary operators
         */
        template<typename T1, typename T2>
        friend auto operator+(const Series<T1>& s1, const Series<T2>& s2)
        -> Series<decltype(s1.at(s1.get_order_min()) + s2.at(s2.get_order_min()))>;
        template<typename T1, typename T2>
        friend auto operator+(const Series<T1>& series, const T2& scalar)
        -> Series<decltype(series.at(series.get_order_min()) + scalar)>;
        template<typename T1, typename T2>
        friend auto operator+(const T1& scalar, const Series<T2>& series)
        -> Series<decltype(scalar + series.at(series.get_order_min()))>;

        template<typename T1, typename T2>
        friend auto operator-(const Series<T1>& s1, const Series<T2>& s2)
        -> Series<decltype(s1.at(s1.get_order_min()) + s2.at(s2.get_order_min()))>;
        template<typename T1, typename T2>
        friend auto operator-(const Series<T1>& series, const T2& scalar)
        -> Series<decltype(series.at(series.get_order_min()) + scalar)>;
        template<typename T1, typename T2>
        friend auto operator-(const T1& scalar, const Series<T2>& series)
        -> Series<decltype(scalar + series.at(series.get_order_min()))>;

        // (\Sum_i=a1^b1 c_i * eps^i) * (\Sum_j=a2^b2 c'_j * eps^j)
        template<typename T2>
        friend auto operator*(const Series& s1, const Series<T2>& s2)
        -> Series<decltype(s1.at(s1.get_order_min()) * s2.at(s2.get_order_min()))>
        {
            using Tout = decltype(s1.at(s1.get_order_min()) * s2.at(s2.get_order_min()));

            if (s1.expansion_parameter != s2.expansion_parameter)
                throw expansion_parameter_mismatch_error("\"" + s1.expansion_parameter + "\" != \"" + s2.expansion_parameter + "\"");

            // Assume series are not truncated
            bool truncated_above = false;
            int order_min = s1.get_order_min() + s2.get_order_min();
            int order_max = s1.get_order_max() + s2.get_order_max();
            size_t current_index;

            if ( s1.get_truncated_above() )
            {
                truncated_above = true;
                order_max = std::min( order_max , s1.get_order_max() + s2.get_order_min() );
            }
            if ( s2.get_truncated_above() )
            {
                truncated_above = true;
                order_max = std::min( order_max , s1.get_order_min() + s2.get_order_max() );
            }

            // Perform multiplication
            std::vector<Tout> content;
            content.reserve(order_max-order_min+1);
            for ( int i = s1.order_min; i < s1.order_max + 1; i++ )
            {
                for ( int j = s2.get_order_min(); j < s2.get_order_max() + 1; j++ )
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

            return Series<Tout>(order_min, order_max, content, truncated_above, s1.expansion_parameter /* equality of the expansion parameter was checked earlier */);
        }

        // d * \Sum_i=a1^b1 c_i * eps^i = \Sum_i=a1^b1 d* c_i * eps^i
        // Note: "val*Series" is defined below assuming that is is equal to "Series*val"
        template<typename T2>
        friend auto operator*(const Series& s1, const T2& val)
        -> Series<decltype(s1.at(s1.get_order_min()) * val)>
        {
            using Tout = decltype(s1.at(s1.get_order_min()) * val);
            std::vector<Tout> content;
            content.reserve(s1.order_max-s1.order_min+1);
            for ( int i = s1.order_min; i < s1.order_max + 1; i++ )
                content.push_back(val * s1.at(i));

            return Series<Tout>(s1.order_min, s1.order_max, content, s1.truncated_above, s1.expansion_parameter);
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
               bool truncated_above = true,
               const std::string expansion_parameter = "x"
               ) :
        order_min(order_min), order_max(order_max),
        content(content), truncated_above(truncated_above),
        expansion_parameter(expansion_parameter)
        {
            if (order_min > order_max)
                throw std::invalid_argument("\"order_min\" (" + std::to_string(order_min) + ") > \"order_max\" (" + std::to_string(order_max) + ") in series constructor.");
            if ( content.size() != (order_max-order_min+1) )
                throw std::invalid_argument("Incorrect number of series coefficients. got: " + std::to_string(content.size()) + ", expected: " + std::to_string(order_max-order_min+1));
        }

        // TODO: converting constructor

    };

    /*
     *  Binary operators
     */
    template<typename T1, typename T2>
    inline auto operator+(const Series<T1>& s1, const Series<T2>& s2)
    -> Series<decltype(s1.at(s1.get_order_min()) + s2.at(s2.get_order_min()))>
    {
        return s1.template sub_or_add<false>(s1,s2);
    };
    template<typename T1, typename T2>
    inline auto operator+(const Series<T1>& series, const T2& scalar)
    -> Series<decltype(series.at(series.get_order_min()) + scalar)>
    {
        return series.template sub_or_add<false>(series,Series<T2>{0,0,{scalar},false,series.expansion_parameter});
    };
    template<typename T1, typename T2>
    inline auto operator+(const T1& scalar, const Series<T2>& series)
    -> Series<decltype(scalar + series.at(series.get_order_min()))>
    {
        return series.template sub_or_add<false>(Series<T1>{0,0,{scalar},false,series.expansion_parameter}, series);
    };

    template<typename T1, typename T2>
    inline auto operator-(const Series<T1>& s1, const Series<T2>& s2)
    -> Series<decltype(s1.at(s1.get_order_min()) + s2.at(s2.get_order_min()))>
    {
        return s1.template sub_or_add<true>(s1,s2);
    };
    template<typename T1, typename T2>
    inline auto operator-(const Series<T1>& series, const T2& scalar)
    -> Series<decltype(series.at(series.get_order_min()) + scalar)>
    {
        return series.template sub_or_add<true>(series,Series<T2>{0,0,{scalar},false,series.expansion_parameter});
    };
    template<typename T1, typename T2>
    inline auto operator-(const T1& scalar, const Series<T2>& series)
    -> Series<decltype(scalar + series.at(series.get_order_min()))>
    {
        return series.template sub_or_add<true>(Series<T1>{0,0,{scalar},false,series.expansion_parameter}, series);
    };

    /*
     * This operator must be defined outside of the class, otherwise the
     * compiler generates multiple definitions for Series<T1> * Series<T2>
     * (one from "T=T1" and one from "T=T2").
     */
    template<typename T1, typename T2>
    auto operator*(const T1& val, const Series<T2>& s1)
    -> Series<decltype(s1.at(s1.get_order_min()) * val)>
    {
        return s1*val;
    }

}

#endif

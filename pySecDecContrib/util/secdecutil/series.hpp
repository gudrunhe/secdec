#ifndef SecDecUtil_series_hpp_included
#define SecDecUtil_series_hpp_included

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
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
        
        template<typename U>
        struct is_series {
          static bool const value = false;
        };
        template<typename U>
        struct is_series<Series<U>> {
          static bool const value = true;
        };

        template<typename U, typename V = U>
        struct CreateContent
        {
            static V create(const Series<U>& series, const V& scalar = V())
            {
                return scalar;
            }
        };

        template<typename U, typename V>
        struct CreateContent<Series<U>,V>
        {
            static Series<U> create(const Series<Series<U>>& series)
            {
                return {0,0,{CreateContent<U>::create(series.get_content().at(0))},false,series.get_content().at(0).expansion_parameter};
            }
            static auto create(const Series<Series<U>>& series, const V& scalar)
            -> Series<decltype(CreateContent<U,V>::create(series.get_content().at(0), scalar))>
            {
                return {0,0,{CreateContent<U,V>::create(series.get_content().at(0),scalar)},false,series.get_content().at(0).expansion_parameter};
            }
        };

        int order_min;
        int order_max;
        std::vector<T> content;
        bool truncated_above;

        /*
         *  Helper functions
         */

        /*
         * Performs subtraction of two series
         * \Sum_i=a1^b1 c_i * eps^i - \Sum_j=a2^b2 c'_j * eps^j = \Sum_i=min(a1,a2)^max(b1,b2) (c_i - c'_i) * eps^i
         */
        template<typename T1, typename T2>
        static auto subtract(const Series<T1>& s1, const Series<T2>& s2)
        {
            using Tout = typename std::remove_cv<std::common_type_t<T1,T2>>::type;

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
            // Perform subtraction
            for ( int i = order_min; i < order_max + 1; i++)
            {
                if ( s1.has_term(i) && s2.has_term(i) )
                {
                    content.push_back(s1.at(i) - s2.at(i));
                } else if ( s1.has_term(i) )
                {
                    content.push_back(s1.at(i));
                } else if ( s2.has_term(i) )
                {
                    content.push_back(-s2.at(i));
                } else // construct default
                {
                    content.push_back(CreateContent<Tout>::create(s1));
                }
            }
            return Series<Tout>(order_min, order_max, content, truncated_above, s1.expansion_parameter /* equality of the expansion parameter was checked earlier */);
        }

        /*
         * Performs addition of two series
         * \Sum_i=a1^b1 c_i * eps^i + \Sum_j=a2^b2 c'_j * eps^j = \Sum_i=min(a1,a2)^max(b1,b2) (c_i + c'_i) * eps^i
         */
        template<typename T1, typename T2>
        static auto add(const Series<T1>& s1, const Series<T2>& s2)
        {
            using Tout = typename std::remove_cv<std::common_type_t<T1,T2>>::type;

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
            // Perform addition
            for ( int i = order_min; i < order_max + 1; i++)
            {
                if ( s1.has_term(i) && s2.has_term(i) )
                {
                    content.push_back(s1.at(i) + s2.at(i));
                } else if ( s1.has_term(i) )
                {
                    content.push_back(s1.at(i));
                } else if ( s2.has_term(i) )
                {
                    content.push_back(s2.at(i));
                } else // construct default
                {
                    content.push_back(CreateContent<Tout>::create(s1));
                }
            }
            return Series<Tout>(order_min, order_max, content, truncated_above, s1.expansion_parameter /* equality of the expansion parameter was checked earlier */);
        }

        // (\Sum_i=a1^b1 c_i * eps^i) * (\Sum_j=a2^b2 c'_j * eps^j)
        template<typename T1, typename T2>
        static auto multiply_series(const Series<T1>& s1, const Series<T2>& s2)
        {
            using Tout = typename std::remove_cv<std::common_type_t<T1,T2>>::type;

            if (s1.expansion_parameter != s2.expansion_parameter)
                throw expansion_parameter_mismatch_error("\"" + s1.expansion_parameter + "\" != \"" + s2.expansion_parameter + "\"");

            // Determine limits of resulting series
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
                    if ( ( (i+j) >= order_min ) && ( (i+j) <= order_max ) && s1.has_term(i) && s2.has_term(j) )
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
        template<bool from_left, typename T1, typename T2>
        static auto multiply_scalar(const Series<T1>& s1, const T2& val)
        {
            using Tout = typename std::remove_cv<std::common_type_t<T1,T2>>::type;
            std::vector<Tout> content;
            content.reserve(s1.order_max-s1.order_min+1);
            for ( int i = s1.order_min; i < s1.order_max + 1; i++ )
                if (from_left)
                    content.push_back(val * s1.at(i));
                else
                    content.push_back(s1.at(i) * val);
            return Series<Tout>(s1.order_min, s1.order_max, content, s1.truncated_above, s1.expansion_parameter);
        }

        /*
         * (\Sum_i=a1^b1 c_i * eps^i) / (\Sum_j=a2^b2 c'_j * eps^j)
         * when multiplying out "(\Sum_i=a1^b1 c_i * eps^i) * (\Sum_k=a3^b3 x_k * eps^k) = (\Sum_j=a2^b2 c'_j * eps^j)"
         * and comparing coefficients of "eps^<power>" this comes down to solving a triangular system of equations:
         *
         * Suppose the x_k in "\Sum_k=a3^b3 x_k * eps^k = (\Sum_i=a1^b1 c_i * eps^i) / (\Sum_j=a2^b2 c'_j * eps^j)" are unknown, then:
         *
         *   (c_0  * exp^0 + c_1  * eps^1 + c_2  * eps^2 + ...) * (x_0 * exp^0 + x_1 * eps^1 + x_2 * eps^2 + ...)
         * =  c'_0 * exp^0 + c'_1 * eps^1 + c'_2 * eps^2 + ...
         *
         * by comparing coefficients:
         *
         * c_0 * x_0                         = c'_0
         * c_1 * x_0 + c_0 * x_1             = c'_1
         * c_2 * x_0 + c_1 * x_1 + c_0 * x_2 = c'_2
         *
         */
        template<typename T1, typename T2>
        static auto division_series_by_series(const Series<T1>& s1, const Series<T2>& s2)
        {
            using Tout = typename std::remove_cv<std::common_type_t<T1,T2>>::type;

            if (s1.expansion_parameter != s2.expansion_parameter)
                throw expansion_parameter_mismatch_error("\"" + s1.expansion_parameter + "\" != \"" + s2.expansion_parameter + "\"");

            auto s1_content = s1.get_content();
            auto s2_content = s2.get_content();

            // result can only be not truncated if only one order in denominator
            bool truncated_above = s2.get_order_min() == s2.get_order_max() ? s1.get_truncated_above() || s2.get_truncated_above() : true;

            // Determine limits of resulting series
            int order_min = s1.get_order_min() - s2.get_order_min();
            int order_max = order_min + std::max(s1_content.size(),s2_content.size()) - 1;
            if ( s1.get_truncated_above() )
                order_max = std::min( order_max , order_min + int( s1_content.size() ) - 1 );
            if ( s2.get_truncated_above() )
                order_max = std::min( order_max , order_min + int( s2_content.size() ) - 1 );

            // Perform backsubstitution avoiding an explicit zero, since
            // we do not know if "Tout(0)" or "Tout()" is defined. The
            // nested if statements drop terms that have a multiplicative
            // zero.
            std::vector<Tout> content;
            T2 denominator( s2_content.at(0) );
            content.reserve(order_max-order_min+1);
            Tout tmp( s1_content.at(0)/ s2_content.at(0) );
            bool tmp_initialized;

            // solve  "c_i * x_0 + ... + c_0 * x_i = c'_i" for x_i
            for ( int i = 0; i <= order_max-order_min; ++i )
            {
                if (i == 0)
                {
                    // c_0 * x_0 = c'0  <=>  x_0 = c'_0 / c_0
                    // that is how "tmp" is initialized above
                    content.push_back(tmp);

                } else {

                    // "x_i = (  c'_i - c_i * x_0 - ... - c_1 * x_(i-1)  )    /    c_0"

                    tmp_initialized = false;

                    // if "c'_i" available (assume "c'_i = 0" otherwise)
                    if (s1_content.size() > i)
                    {
                        tmp = s1_content.at(i);
                        tmp_initialized = true;
                    }

                    for ( int j = 0; j < i; ++j )
                    {

                        // if "c_(i-j)" available (assume "c_(i-j) = 0" otherwise)
                        if (s2_content.size() > i-j)
                        {
                            if (tmp_initialized)
                            {
                                tmp -= content.at(j) * s2_content.at(i-j);
                            } else {
                                tmp = - content.at(j) * s2_content.at(i-j);
                                tmp_initialized = true;
                            }
                        }
                    }

                    tmp /= denominator; // final division by "c_0"
                    content.push_back(tmp);

                }

            }

            return Series<Tout>(order_min, order_max, content, truncated_above, s1.expansion_parameter /* equality of the expansion parameter was checked earlier */);

        }

        template<typename T1, typename T2>
        static auto division_series_by_scalar(const Series<T1>& s1, const T2& value)
        {
            using Tout = typename std::remove_cv<std::common_type_t<T1,T2>>::type;

            // copy content
            std::vector<Tout> content = s1.get_content();

            // divide every series coefficient by the scalar
            for ( auto& item : content )
                item /= value;

            return Series<Tout>(s1.order_min, s1.order_max, content, s1.truncated_above, s1.expansion_parameter);
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
        const std::vector<T>& get_content() const { return content; }

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

        void push_back(const_reference val) { content.push_back(val); ++order_max; };
        void push_back(T&& val) { content.push_back(std::forward<T>(val)); ++order_max; };
        void pop_back() { content.pop_back(); --order_max; };

        /*
         *  Helper functions
         */
        // Checks if series has a term at order order
        bool has_term(int order) const
        {
            if ( (order >= order_min) && (order <= order_max) )
                return true;

            return false;
        }

        /*
         *  Comparator Operators
         */
        template<typename T2>
        friend bool operator==(const Series& s1, const Series<T2>& s2)
        {
            if (s1.expansion_parameter != s2.expansion_parameter)
                return false;
            if ( s1.get_order_min() != s2.get_order_min() )
                return false;
            if ( s1.get_order_max() != s2.get_order_max() )
                return false;
            if ( s1.get_truncated_above() != s2.get_truncated_above() )
                return false;
            for ( size_t idx = 0 ; idx < s1.get_order_max()-s1.get_order_min()+1 ; ++idx )
                if ( s1.get_content().at(idx) != s2.get_content().at(idx) )
                    return false;
            return true;
        }

        template<typename T2>
        friend bool operator!= (const Series& s1, const Series<T2>& s2)
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

        template<typename Tother>
        Series& operator/=(const Tother& other)
        {
            *this = *this / other;
            return *this;
        }

        /*
         *  Binary operators
         */
        template<typename T2>
        friend auto operator+(const Series& s1, const Series<T2>& s2)
        {
            return s1.template add(s1,s2);
        }
        template<typename T2, typename = typename std::enable_if<!is_series<T2>::value>::type>
        friend auto operator+(const Series& series, const T2& scalar)
        {
            auto content = Series<T>::CreateContent<T,T2>::create(series,scalar);
            return series.template add(series,Series<decltype(content)>{0,0,{content},false,series.expansion_parameter});
        };
        template<typename T2, typename = typename std::enable_if<!is_series<T2>::value>::type>
        friend auto operator+(const T2& scalar, const Series& series)
        {
            auto content = Series<T>::CreateContent<T,T2>::create(series,scalar);
            return series.template add(Series<decltype(content)>{0,0,{content},false,series.expansion_parameter}, series);
        };

        template<typename T2>
        friend auto operator-(const Series& s1, const Series<T2>& s2)
        {
            return s1.template subtract(s1,s2);
        };
        template<typename T2, typename = typename std::enable_if<!is_series<T2>::value>::type>
        friend auto operator-(const Series& series, const T2& scalar)
        {
            auto content = Series<T>::CreateContent<T,T2>::create(series,scalar);
            return series.template subtract(series,Series<decltype(content)>{0,0,{content},false,series.expansion_parameter});
        };
        template<typename T2, typename = typename std::enable_if<!is_series<T2>::value>::type>
        friend auto operator-(const T2& scalar, const Series& series)
        {
            auto content = Series<T>::CreateContent<T,T2>::create(series,scalar);
            return series.template subtract(Series<decltype(content)>{0,0,{content},false,series.expansion_parameter}, series);
        };

        template<typename T2>
        friend auto operator*(const Series& s1, const Series<T2>& s2)
        {
            return s1.template multiply_series(s1,s2);
        }
        template<typename T2, typename = typename std::enable_if<!is_series<T2>::value>::type>
        friend auto operator*(const Series& s1, const T2& val)
        {
            return s1.template multiply_scalar</* from_left = */ false>(s1, val);
        }
        template<typename T2, typename = typename std::enable_if<!is_series<T2>::value>::type>
        friend auto operator*(const T2& val, const Series& s1)
        {
            return s1.template multiply_scalar</* from_left = */ true>(s1, val);
        }

        template<typename T2>
        friend auto operator/(const Series& s1, const Series<T2>& s2)
        {
            return s1.template division_series_by_series(s1,s2);
        }
        template<typename T2, typename = typename std::enable_if<!is_series<T2>::value>::type>
        friend auto operator/(const Series& s1, const T2& val)
        {
            return s1.template division_series_by_scalar(s1,val);
        }
        template<typename T2, typename = typename std::enable_if<!is_series<T2>::value>::type>
        friend auto operator/(const T2& val, const Series& s1)
        {
            auto content = Series<T>::CreateContent<T,T2>::create(s1,val);
            return s1.template division_series_by_series(Series<decltype(content)>{0,0,{content},false,s1.expansion_parameter},s1);
        }

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

        // converting constructor
        template<typename U>
        Series(const Series<U>& s) :
        order_min(s.get_order_min()), order_max(s.get_order_max()),
        truncated_above(s.get_truncated_above()),
        expansion_parameter(s.expansion_parameter)
        {
            content.reserve(order_max-order_min+1);
            for ( const auto& item : s.get_content() )
                content.push_back(item);
        }

    };

}

namespace std {

    template <typename T1, typename T2>
    struct common_type<secdecutil::Series<T1>, T2> {
        using type = secdecutil::Series<typename std::common_type<T1, T2>::type>;
    };

    template <typename T1, typename T2>
    struct common_type<secdecutil::Series<T1>, secdecutil::Series<T2>> {
        using type = secdecutil::Series<typename std::common_type<T1, T2>::type>;
    };

}


#endif

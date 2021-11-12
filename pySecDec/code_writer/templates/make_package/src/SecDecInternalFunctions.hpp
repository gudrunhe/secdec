#ifndef %(name)s_src_SecDecInternalFunctions_hpp_included
#define %(name)s_src_SecDecInternalFunctions_hpp_included

#include "%(name)s.hpp"

#include <cmath>
#ifdef SECDEC_WITH_CUDA
    #include <thrust/complex.h>
#else
    #include <complex>
#endif
#include <stdexcept>
#include <string>

#ifdef SECDEC_WITH_CUDA
    #define STD thrust
#else
    #define STD std
#endif

namespace %(name)s
{
    // imaginary unit
    #ifdef SECDEC_WITH_CUDA
        // FORM notation
        #define i_ complex_t{0,1}
        // sympy notation
        #define I complex_t{0,1}
    #else
        constexpr complex_t i_{0,1}; // FORM notation
        constexpr complex_t I = i_; // sympy notation
    #endif


    // required functions
    // --{

    /*
     * We need the nonstandard continuation on the negative x-axis;
     * i.e. log(-1) = -i*pi.
     */
    // use std::log(real_t) but override log(complex_t)
    #define %(name)s_contour_deformation %(contour_deformation)i
    #define %(name)s_has_complex_parameters %(have_complex_parameters)i
    #define %(name)s_enforce_complex_return_type %(enforce_complex_return_type)i
    #if %(name)s_has_complex_parameters || %(name)s_contour_deformation || %(name)s_enforce_complex_return_type
      #ifdef SECDEC_WITH_CUDA
        __host__ __device__
      #endif
        inline complex_t log(complex_t arg)
        {
            if (arg.imag() == 0)
                arg = complex_t(arg.real(),-0.);
            return STD::log(arg);
        }
    #else
      #ifdef SECDEC_WITH_CUDA
        __host__ __device__
      #endif
        inline real_t log(real_t arg)
        {
            if (arg < 0)
            {
                #ifdef SECDEC_WITH_CUDA
                    return std::nan("");
                #else
                    std::string error_message;
                    error_message += "Encountered \"log(<negative real>)\" in a real-valued integrand function of \"%(name)s\". ";
                    error_message += "Try to enforce complex return values for the generated integrands; i.e. set ";
                    error_message += "\"enforce_complex=True\" in the corresponding call to \"loop_package\" or \"make_package\".";
                    throw std::domain_error(error_message);
                #endif
            }
            return std::log(arg);
        }
    #endif

    /*
     * We do not want to use "std::pow(double, int)" because the g++ compiler
     * casts it to "pow(double, double)" which is extremely inefficient but
     * demanded by the c++ standard.
     *
     * Note: Using std::pow and the flags "-O2" and "-ffast-math" with g++ is even faster,
     *       but "-ffast-math" is g++ specific and allows "unsafe math optimizations".
     *       The intel compiler produces code that runs faster when using std::pow and
     *       "-O2" than with this function.
     *       However, the c++ standard requires that the second argument of
     *       "std::pow(double, int)" is casted to double. To comply with the standard, we
     *       decided to implement our own optimized power function rather than relying on
     *       the compiler to perform optimizations possibly disallowed by the c++ standard.
     *       Playing around with "std::pow" and the aforementioned switches is nevertheless
     *       worth a try in practical applications where high performance is needed.
     */
    template <typename Tbase> inline
    #ifdef SECDEC_WITH_CUDA
      __host__ __device__
    #endif
    Tbase SecDecInternalPow(Tbase base, int exponent)
    {

        #ifndef SECDEC_WITH_CUDA
            if (exponent > 1024 or exponent < -1024)
                return std::pow(base, exponent);
        #endif

        if (exponent < 0)
            return Tbase(1)/SecDecInternalPow(base, -exponent);

        else if (exponent == 0)
            return Tbase(1);

        else if (exponent == 1)
            return base;

        else if (exponent == 2)
            return base * base;

        else if (exponent == 3)
            return base * base * base;

        else if (exponent == 4)
        {
            Tbase result = base;
            result *= result;
            result *= result;
            return result;
        }

        else if (exponent == 5)
        {
            Tbase result = base;
            result *= result;
            result *= result;
            return result * base;
        }

        else if (exponent == 6)
        {
            Tbase result = base * base * base;
            return result * result;
        }

        else if (exponent == 7)
        {
            Tbase result = base * base * base;
            return result * result * base;
        }

        else if (exponent == 8)
        {
            Tbase result = base;
            result *= result;
            result *= result;
            return result * result;
        }

        else if (exponent == 16)
        {
            Tbase tmp = base * base;
            tmp *= tmp;
            tmp *= tmp;
            tmp *= tmp;
            return tmp;
        }

        unsigned half_exponent = exponent / 2;
        Tbase out = SecDecInternalPow(base, half_exponent);

        out *= out;
        if (2 * half_exponent == exponent) // exponent is even
            return out;
        else // exponent is odd --> need another factor of the base due to integer division above
            return out * base;
    }

    #ifdef SECDEC_WITH_CUDA
      __host__ __device__
    #endif
    real_t inline pow(real_t x, int y)
    {
        return SecDecInternalPow(x, y);
    }
    #ifdef SECDEC_WITH_CUDA
      __host__ __device__
    #endif
    complex_t inline pow(complex_t x, int y)
    {
        return SecDecInternalPow(x, y);
    }
    #if %(name)s_has_complex_parameters || %(name)s_contour_deformation || %(name)s_enforce_complex_return_type
        template <typename Tbase, typename Texponent>
      #ifdef SECDEC_WITH_CUDA
        __host__ __device__
      #endif
        complex_t pow(Tbase base, Texponent exponent)
        {
            #ifdef SECDEC_WITH_CUDA
                complex_t cbase = base;
                if (cbase.imag() == 0)
                    return thrust::pow( complex_t(cbase.real(),-0.) , exponent );
                return thrust::pow(cbase, exponent);
            #else
                if (std::imag(base) == 0)
                    return std::pow( complex_t(STD::real(base),-0.) , exponent );
                return std::pow(base, exponent);
            #endif
        }
    #else
      #ifdef SECDEC_WITH_CUDA
        __host__ __device__
      #endif
        inline real_t pow(real_t base, real_t exponent)
        {
            if (base < 0)
            {
                #ifdef SECDEC_WITH_CUDA
                    return std::nan("");
                #else
                    std::string error_message;
                    error_message += "Encountered \"pow(<negative real>, <rational>)\" in a real-valued integrand function of ";
                    error_message += "\"%(name)s\". Try to enforce complex return values for the generated integrands; i.e. set ";
                    error_message += "\"enforce_complex=True\" in the corresponding call to \"loop_package\" or \"make_package\".";
                    throw std::domain_error(error_message);
                #endif
            }
            return std::pow(base, exponent);
        }
    #endif

    /*
     * Overload binary arithmetic operators between "int" and "complex_t"
     * by converting "int" to "real_t".
     */
    #ifdef SECDEC_WITH_CUDA
      #define INT_COMPLEX(OPERATOR) \
          inline __host__ __device__ complex_t operator OPERATOR (const int x, const complex_t y) \
          { \
              return static_cast<real_t>(x) OPERATOR y; \
          } \
          inline __host__ __device__ complex_t operator OPERATOR (const complex_t x, const int y) \
          { \
              return x OPERATOR static_cast<real_t>(y); \
          }
    #else
      #define INT_COMPLEX(OPERATOR) \
          inline complex_t operator OPERATOR (const int x, const complex_t y) \
          { \
              return static_cast<real_t>(x) OPERATOR y; \
          } \
          inline complex_t operator OPERATOR (const complex_t x, const int y) \
          { \
              return x OPERATOR static_cast<real_t>(y); \
          }
    #endif
    INT_COMPLEX(+)
    INT_COMPLEX(-)
    INT_COMPLEX(*)
    INT_COMPLEX(/)
    #undef INT_COMPLEX

    #ifdef SECDEC_WITH_CUDA
        #define SecDecFn __host__ __device__ static inline
     #else
        #define SecDecFn static inline
    #endif
    SecDecFn real_t SecDecInternalRealPart(const real_t x) { return x; }
    SecDecFn real_t SecDecInternalRealPart(const complex_t x) { return x.real(); }
    SecDecFn real_t SecDecInternalImagPart(const real_t x) { return 0; }
    SecDecFn real_t SecDecInternalImagPart(const complex_t x) { return x.imag(); }
    SecDecFn complex_t SecDecInternalI(const real_t x) { return complex_t{0, x}; }
    SecDecFn complex_t SecDecInternalI(const complex_t x) { return complex_t{-x.imag(), x.real()}; }

    #undef %(name)s_contour_deformation
    #undef %(name)s_has_complex_parameters
    #undef %(name)s_enforce_complex_return_type
    #undef SecDecFn

    // --}

};

#undef STD

#define SecDecInternalDenominator(x) 1.0/(x)
#define SecDecInternalLog(x) log(x)
#define SecDecInternalPow(x, n) pow(x, n)

#ifdef SECDEC_WITH_CUDA
#define SecDecInternalAbs(x) thrust::abs(complex_t{x})
#else
#define SecDecInternalAbs(x) std::abs(x)
#endif

#define likely(x) (x)
#define unlikely(x) (x)

#define SecDecInternalOutputDeformationParameters(i, v) output_deformation_parameters[i] = (v);
#define SecDecInternalSignCheckErrorPositivePolynomial(id) \
    { \
      secdecutil::ResultInfo current_result; \
      current_result.return_value = secdecutil::ResultInfo::ReturnValue::sign_check_error_positive_polynomial; \
      current_result.signCheckId = (id); \
      result_info->fill_if_empty_threadsafe(current_result); \
    };
#define SecDecInternalSignCheckErrorContourDeformation(id) \
    { \
      secdecutil::ResultInfo current_result; \
      current_result.return_value = secdecutil::ResultInfo::ReturnValue::sign_check_error_contour_deformation; \
      current_result.signCheckId = (id); \
      result_info->fill_if_empty_threadsafe(current_result); \
    };

#if defined(__GNUC__) || defined(__NVCC__)
#define restrict __restrict__
#else
#define restrict
#endif

#endif

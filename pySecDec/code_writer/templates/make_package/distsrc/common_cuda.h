#pragma nv_diag_suppress 177 // disable the "function was declared but never referenced" warning

#include <cinttypes>
#include <cub/block/block_reduce.cuh>
#include <math_constants.h>
#include <thrust/complex.h>

#define likely(x) __builtin_expect((x), 1)
#define unlikely(x) __builtin_expect((x), 0)
#define mathfn __device__ static inline

typedef double real_t;
typedef thrust::complex<real_t> complex_t;
#define REAL_NAN CUDART_NAN

#if SECDEC_RESULT_IS_COMPLEX 
    typedef complex_t result_t;
    mathfn complex_t SecDecInternalLog(const real_t x)
    { return (x >= 0) ? log(x) : complex_t{log(-x), -M_PI}; };
    mathfn complex_t SecDecInternalLog(const complex_t x)
    { return (x.imag() == 0) ? SecDecInternalLog(x.real()) : log(x); };
    mathfn complex_t SecDecInternalPow(real_t x, real_t n)
    { return (x >= 0) ? pow(x, n) : conj(pow(complex_t{x}, n)); }
    mathfn complex_t SecDecInternalPow(complex_t x, real_t n)
    { return conj(pow(conj(x), n)); }
#else
    typedef real_t result_t;
    mathfn real_t SecDecInternalLog(real_t x)
    { return log(x); };
    mathfn real_t SecDecInternalExp(real_t x)
    { return exp(x); };
    mathfn real_t SecDecInternalPow(real_t x, real_t n)
    { return pow(x, n); }
#endif

mathfn double SecDecInternalSqr(const double x) { return x*x; }
mathfn complex_t SecDecInternalSqr(const complex_t x) { return x*x; }
mathfn real_t SecDecInternalRealPart(const real_t x) { return x; }
mathfn real_t SecDecInternalRealPart(const complex_t x) { return x.real(); }
mathfn real_t SecDecInternalImagPart(const real_t x) { return 0; }
mathfn real_t SecDecInternalImagPart(const complex_t x) { return x.imag(); }
mathfn complex_t SecDecInternalI(const real_t x) { return complex_t{0, x}; }
mathfn complex_t SecDecInternalI(const complex_t x) { return complex_t{-x.imag(), x.real()}; }

mathfn real_t exp(int n) { return exp(real_t(n)); }

#define SecDecInternalNPow(x, n) SecDecInternalNPowTemplate<n>(x)
template<unsigned n, typename T> inline T
SecDecInternalNPowTemplate(const T &x)
{
    if constexpr ((n%%2) == 0) {
        return SecDecInternalSqr(SecDecInternalNPowTemplate<n/2,T>(x));
    } else {
        return SecDecInternalSqr(SecDecInternalNPowTemplate<n/2,T>(x))*x;
    }
}
template<> inline real_t SecDecInternalNPowTemplate<1,real_t>(const real_t &x) { return x; };
template<> inline complex_t SecDecInternalNPowTemplate<1,complex_t>(const complex_t &x) { return x; };

#define SecDecInternalQuo(n, d) (((real_t)(n))/(d))
#define SecDecInternalDenominator(x) ((real_t)1/(x))
#define i_ (complex_t{0,1})

mathfn real_t none_f(real_t x) { return x; }
mathfn real_t none_w(real_t x) { return 1; }
mathfn real_t baker_f(real_t x) { auto a = 2*x; auto b = 2-a; return (a <= b) ? a : b; }
mathfn real_t baker_w(real_t x) { return 1; }
mathfn real_t korobov1x1_w(const real_t x) { auto u = (1-x)*x; return 6*u; }
mathfn real_t korobov1x1_f(const real_t x) {
    auto y = (x <= 0.5) ? x : 1-x;
    auto y2 = SecDecInternalSqr(y);
    auto h = y2*(3 - 2*y);
    return (x <= 0.5) ? h : 1-h;
}
mathfn real_t korobov2x2_w(const real_t x) { auto u = (1-x)*x; return 30*SecDecInternalSqr(u); }
mathfn real_t korobov2x2_f(const real_t x) {
    auto y = (x <= 0.5) ? x : 1-x;
    auto y3 = y*SecDecInternalSqr(y);
    auto h = y3*(10 + y*(-15 + 6*y));
    return (x <= 0.5) ? h : 1-h;
}
mathfn real_t korobov3x3_w(const real_t x) { auto u = (1-x)*x; return 140*u*SecDecInternalSqr(u); }
mathfn real_t korobov3x3_f(const real_t x) {
    auto y = (x <= 0.5) ? x : 1-x;
    auto y4 = SecDecInternalSqr(SecDecInternalSqr(y));
    auto h = y4*(35 + y*(-84 + (70 - 20*y)*y));
    return (x <= 0.5) ? h : 1-h;
}
mathfn real_t korobov4x4_w(const real_t x) { auto u = (1-x)*x; return 630*SecDecInternalSqr(SecDecInternalSqr(u)); }
mathfn real_t korobov4x4_f(const real_t x) {
    auto y = (x <= 0.5) ? x : 1-x;
    auto y5 = y*SecDecInternalSqr(SecDecInternalSqr(y));
    auto h = y5*(126 + y*(-420 + y*(540 + y*(-315 + 70*y))));
    return (x <= 0.5) ? h : 1-h;
}
mathfn real_t korobov5x5_w(const real_t x) { auto u = (1-x)*x; return 2772*u*SecDecInternalSqr(SecDecInternalSqr(u)); }
mathfn real_t korobov5x5_f(const real_t x) {
    auto y = (x <= 0.5) ? x : 1-x;
    auto y6 = SecDecInternalSqr(y*SecDecInternalSqr(y));
    auto h = y6*(462 + y*(-1980 + y*(3465 + y*(-3080 + (1386 - 252*y)*y))));
    return (x <= 0.5) ? h : 1-h;
}
mathfn real_t korobov6x6_w(const real_t x) { auto u = (1-x)*x; return 12012*SecDecInternalSqr(u*SecDecInternalSqr(u)); }
mathfn real_t korobov6x6_f(const real_t x) {
    auto y = (x <= 0.5) ? x : 1-x;
    auto y7 = y*SecDecInternalSqr(y*SecDecInternalSqr(y));
    auto h = y7*(1716 + y*(-9009 + y*(20020 + y*(-24024 + y*(16380 + y*(-6006 + 924*y))))));
    return (x <= 0.5) ? h : 1-h;
}

mathfn uint64_t mulmod(uint64_t a, uint64_t b, uint64_t k) {
    // assume 0 <= a,b <= k < 2^53
    if (k <= 3037000499) { // floor(Int, sqrt(2^63-1))
        return (a*b) %% k;
    } else {
        auto x = static_cast<double>(a);
        auto c = static_cast<uint64_t>( (x*b) / k );
        auto r = static_cast<int64_t>( (a*b) - (c*k) ) %% static_cast<int64_t>(k);
        return r < 0 ? static_cast<uint64_t>(r+k) : static_cast<uint64_t>(r);
    }
}

mathfn real_t warponce(const real_t a, const real_t b) { return a >= b ? a - b : a; }
mathfn uint64_t warponce_i(const uint64_t a, const uint64_t b) { return a >= b ? a - b : a; }

#pragma nv_diag_suppress 177 // disable the "function was declared but never referenced" warning

#include <cinttypes>
#include <cub/block/block_reduce.cuh>
#include <thrust/complex.h>

#define mathfn __device__ static inline

typedef double real_t;
typedef thrust::complex<real_t> complex_t;

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
    mathfn real_t SecDecInternalPow(real_t x, real_t n)
    { return pow(x, n); }
#endif

mathfn real_t SecDecInternalRealPart(const real_t x) { return x; }
mathfn real_t SecDecInternalRealPart(const complex_t x) { return x.real(); }
mathfn real_t SecDecInternalImagPart(const real_t x) { return 0; }
mathfn real_t SecDecInternalImagPart(const complex_t x) { return x.imag(); }
mathfn complex_t SecDecInternalI(const real_t x) { return complex_t{0, x}; }
mathfn complex_t SecDecInternalI(const complex_t x) { return complex_t{-x.imag(), x.real()}; }

mathfn real_t exp(int n) { return exp(real_t(n)); }

mathfn real_t clamp01(const real_t &a)
{ real_t b = a < 1 ? a : 1; return b > 0 ? b : 0; }

#define SecDecInternalDenominator(x) (1.0/(x))
#define i_ (complex_t{0,1})

#define likely(x) __builtin_expect((x), 1)
#define unlikely(x) __builtin_expect((x), 0)

mathfn real_t none_f(real_t x) { return x; }
mathfn real_t none_w(real_t x) { return 1; }
mathfn real_t baker_f(real_t x) { auto a = 2*x; auto b = 2-a; return (a <= b) ? a : b; }
mathfn real_t baker_w(real_t x) { return 1; }
mathfn real_t korobov1x1_f(real_t x) { return x*x*((-2)*x + 3); }
mathfn real_t korobov1x1_w(real_t x) { return (1 - x)*x*6; }
mathfn real_t korobov2x2_f(real_t x) { return x*x*x*((6*x - 15)*x + 10); }
mathfn real_t korobov2x2_w(real_t x) { auto xx = (1 - x)*x; return xx*xx*30; }
mathfn real_t korobov3x3_f(real_t x) { auto xx = x*x; return xx*xx*((((-20)*x + 70)*x - 84)*x + 35); }
mathfn real_t korobov3x3_w(real_t x) { auto xx = (1 - x)*x; return xx*xx*xx*140; }
mathfn real_t korobov4x4_f(real_t x) { auto xx = x*x; return xx*xx*x*((((70*x - 315)*x + 540)*x - 420)*x + 126); }
mathfn real_t korobov4x4_w(real_t x) { auto xx = (1 - x)*x; auto xx2 = xx*xx; return xx2*xx2*630; }
mathfn real_t korobov5x5_f(real_t x) { auto x3 = x*x*x; return x3*x3*((((((-252)*x + 1386)*x - 3080)*x + 3465)*x - 1980)*x + 462); }
mathfn real_t korobov5x5_w(real_t x) { auto xx = (1 - x)*x; auto xx2 = xx*xx; return xx2*xx2*xx*2772; }

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

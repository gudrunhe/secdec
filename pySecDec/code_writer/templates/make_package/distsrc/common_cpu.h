#include <cinttypes>
#include <complex>

typedef int64_t int_t;
typedef double real_t;
typedef std::complex<real_t> complex_t;

#if defined(__GNUC__) || defined(__clang__)
    #define GNU_VECTORS 1
    struct alignas(32) realvec_t { real_t x __attribute__((vector_size(32))); };
    struct alignas(32) complexvec_t { realvec_t re, im; };
    #define likely(x) __builtin_expect((x), 1)
    #define unlikely(x) __builtin_expect((x), 0)
    #define restrict __restrict__
#else
    #define GNU_VECTORS 0
    struct alignas(32) realvec_t { real_t x[4]; };
    struct alignas(32) complexvec_t { realvec_t re, im; };
    #define likely(x) (x)
    #define unlikely(x) (x)
    #define restrict
#endif

#define mathfn static inline

// Misc

mathfn int_t warponce_i(const int_t a, const int_t b)
{ int_t ab = a - b; return ab >= 0 ? ab : a; }

#define SecDecInternalDenominator(x) 1/(x)

mathfn real_t SecDecInternalRealPart(const real_t &a) { return a; }
mathfn real_t SecDecInternalImagPart(const real_t &a) { return 0; }
mathfn complex_t SecDecInternalI(const real_t &a) { return complex_t{0, a}; }
mathfn complex_t SecDecInternalI(const complex_t &a) { return complex_t{-a.imag(), a.real()}; }

static unsigned long mulmod(unsigned long a, unsigned long b, unsigned long k) {
    // assume 0 <= a,b <= k < 2^53
    if (k <= 3037000499) { // floor(Int, sqrt(2^63-1))
        return (a*b) %% k;
    } else {
        auto x = static_cast<double>(a);
        auto c = static_cast<unsigned long>( (x*b) / k );
        auto r = static_cast<signed long>( (a*b) - (c*k) ) %% static_cast<signed long>(k);
        return r < 0 ? static_cast<unsigned long>(r+k) : static_cast<unsigned long>(r);
    }
}

// Operators for complex (op) int

#define INT_COMPLEX(OP) \
  mathfn complex_t operator OP(const int x, const complex_t &y) { return (real_t)(x) OP y; } \
  mathfn complex_t operator OP(const complex_t &x, const int y) { return x OP (real_t)(y); }

INT_COMPLEX(+)
INT_COMPLEX(-)
INT_COMPLEX(*)
INT_COMPLEX(/)

// Real vectors

#define REALVEC_CONST(c) (realvec_t{{c,c,c,c}})
#define REALVEC_ZERO REALVEC_CONST(0)

#if GNU_VECTORS

    #define DEF_OPERATOR(ret_t, fname, arg1_t, arg2_t, op) \
        mathfn ret_t fname(const arg1_t &a, const arg2_t &b) { return ret_t{ a.x op b.x }; }

    #define DEF_SCALAR_OPERATOR(ret_t, fname, vec_t, scalar_t, op) \
        mathfn ret_t fname(const scalar_t &a, const vec_t &b) { return ret_t{ a op b.x }; } \
        mathfn ret_t fname(const vec_t &a, const scalar_t &b) { return ret_t{ a.x op b }; }

#else

    #define DEF_OPERATOR(ret_t, fname, arg1_t, arg2_t, op) \
        mathfn ret_t fname(const arg1_t &a, const arg2_t &b) \
        { return ret_t{{ a.x[0] op b.x[0], a.x[1] op b.x[1], a.x[2] op b.x[2], a.x[3] op b.x[3] }}; }

    #define DEF_SCALAR_OPERATOR(ret_t, fname, vec_t, scalar_t, op) \
        mathfn ret_t fname(const scalar_t &a, const vec_t &b) \
        { return ret_t{{ a op b.x[0], a op b.x[1], a op b.x[2], a op b.x[3] }}; } \
        mathfn ret_t fname(const vec_t &a, const scalar_t &b) \
        { return ret_t{{ a.x[0] op b, a.x[1] op b, a.x[2] op b, a.x[3] op b }}; }

#endif

#define DEF_FUNCTION(ret_t, fname, arg_t, f) \
    mathfn ret_t fname(const arg_t &a) \
    { return ret_t{{ f(a.x[0]), f(a.x[1]), f(a.x[2]), f(a.x[3]) }}; }
#define DEF_SCALAR_BOOL_OPERATOR(fname, vec_t, scalar_t, op, redop) \
    mathfn bool fname(const vec_t &a, const scalar_t &b) \
    { return (a.x[0] op b) redop (a.x[1] op b) redop (a.x[2] op b) redop (a.x[3] op b); }

DEF_OPERATOR(realvec_t, operator +, realvec_t, realvec_t, +)
DEF_OPERATOR(realvec_t, operator -, realvec_t, realvec_t, -)
DEF_OPERATOR(realvec_t, operator *, realvec_t, realvec_t, *)
DEF_OPERATOR(realvec_t, operator /, realvec_t, realvec_t, /)
DEF_SCALAR_OPERATOR(realvec_t, operator +, realvec_t, real_t, +)
DEF_SCALAR_OPERATOR(realvec_t, operator -, realvec_t, real_t, -)
DEF_SCALAR_OPERATOR(realvec_t, operator *, realvec_t, real_t, *)
DEF_SCALAR_OPERATOR(realvec_t, operator /, realvec_t, real_t, /)
DEF_SCALAR_BOOL_OPERATOR(operator >, realvec_t, real_t, >, ||)
DEF_SCALAR_BOOL_OPERATOR(operator <, realvec_t, real_t, <, ||)
DEF_FUNCTION(realvec_t, operator -, realvec_t, -)
DEF_FUNCTION(realvec_t, SecDecInternalAbs, realvec_t, std::abs)

mathfn realvec_t SecDecInternalRealPart(const realvec_t &a) { return a; }
mathfn realvec_t SecDecInternalImagPart(const realvec_t &a) { return REALVEC_ZERO; }

mathfn real_t componentsum(const realvec_t &a)
{ return a.x[0] + a.x[1] + a.x[2] + a.x[3]; }

mathfn realvec_t warponce(const realvec_t &a, const real_t b)
{ realvec_t ab = a - b;
  return realvec_t{{ ab.x[0] >= 0 ? ab.x[0] : a.x[0],
                     ab.x[1] >= 0 ? ab.x[1] : a.x[1],
                     ab.x[2] >= 0 ? ab.x[2] : a.x[2],
                     ab.x[3] >= 0 ? ab.x[3] : a.x[3] }}; }

mathfn realvec_t korobov1_f(const realvec_t &x) { return x*x*((-2)*x + 3); }
mathfn realvec_t korobov1_w(const realvec_t &x) { return (1 - x)*x*6; }
mathfn realvec_t korobov2_f(const realvec_t &x) { return x*x*x*((6*x - 15)*x + 10); }
mathfn realvec_t korobov2_w(const realvec_t &x) { auto xx = (1 - x)*x; return xx*xx*30; }
mathfn realvec_t korobov3_f(const realvec_t &x) { auto xx = x*x; return xx*xx*((((-20)*x + 70)*x - 84)*x + 35); }
mathfn realvec_t korobov3_w(const realvec_t &x) { auto xx = (1 - x)*x; return xx*xx*xx*140; }
mathfn realvec_t korobov4_f(const realvec_t &x) { auto xx = x*x; return xx*xx*x*((((70*x - 315)*x + 540)*x - 420)*x + 126); }
mathfn realvec_t korobov4_w(const realvec_t &x) { auto xx = (1 - x)*x; auto xx2 = xx*xx; return xx2*xx2*630; }
mathfn realvec_t korobov5_f(const realvec_t &x) { auto x3 = x*x*x; return x3*x3*((((((-252)*x + 1386)*x - 3080)*x + 3465)*x - 1980)*x + 462); }
mathfn realvec_t korobov5_w(const realvec_t &x) { auto xx = (1 - x)*x; auto xx2 = xx*xx; return xx2*xx2*xx*2772; }

// Complex vectors

#define COMPLEXVEC_ZERO (complexvec_t{REALVEC_ZERO, REALVEC_ZERO})

mathfn complexvec_t operator +(const complexvec_t &a, const complexvec_t &b)
{ return complexvec_t{a.re + b.re, a.im + b.im}; }
mathfn complexvec_t operator +(const complexvec_t &a, const realvec_t &b)
{ return complexvec_t{a.re + b, a.im}; }
mathfn complexvec_t operator +(const realvec_t &a, const complexvec_t &b)
{ return complexvec_t{a + b.re, b.im}; }
mathfn complexvec_t operator +(const complexvec_t &a, const real_t &b)
{ return complexvec_t{a.re + b, a.im}; }
mathfn complexvec_t operator +(const real_t &a, const complexvec_t &b)
{ return complexvec_t{a + b.re, b.im}; }
mathfn complexvec_t operator +(const complexvec_t &a, const complex_t &b)
{ return complexvec_t{a.re + b.real(), a.im + b.imag()}; }
mathfn complexvec_t operator +(const complex_t &a, const complexvec_t &b)
{ return complexvec_t{a.real() + b.re, a.imag() + b.im}; }
mathfn complexvec_t operator +(const realvec_t &a, const complex_t &b)
{ return complexvec_t{a + b.real(), REALVEC_CONST(b.imag())}; }
mathfn complexvec_t operator +(const complex_t &a, const realvec_t &b)
{ return complexvec_t{a.real() + b, REALVEC_CONST(a.imag())}; }

mathfn complexvec_t operator -(const complexvec_t &a)
{ return complexvec_t{-a.re, -a.im}; }

mathfn complexvec_t operator -(const complexvec_t &a, const complexvec_t &b)
{ return complexvec_t{a.re - b.re, a.im - b.im}; }
mathfn complexvec_t operator -(const complexvec_t &a, const realvec_t &b)
{ return complexvec_t{a.re - b, a.im}; }
mathfn complexvec_t operator -(const realvec_t &a, const complexvec_t &b)
{ return complexvec_t{a - b.re, -b.im}; }
mathfn complexvec_t operator -(const complexvec_t &a, const real_t &b)
{ return complexvec_t{a.re - b, a.im}; }
mathfn complexvec_t operator -(const real_t &a, const complexvec_t &b)
{ return complexvec_t{a - b.re, -b.im}; }
mathfn complexvec_t operator -(const complexvec_t &a, const complex_t &b)
{ return complexvec_t{a.re - b.real(), a.im - b.imag()}; }
mathfn complexvec_t operator -(const complex_t &a, const complexvec_t &b)
{ return complexvec_t{a.real() - b.re, a.imag() - b.im}; }
mathfn complexvec_t operator -(const realvec_t &a, const complex_t &b)
{ return complexvec_t{a - b.real(), REALVEC_CONST(-b.imag())}; }
mathfn complexvec_t operator -(const complex_t &a, const realvec_t &b)
{ return complexvec_t{a.real() - b, REALVEC_CONST(a.imag())}; }

mathfn complexvec_t operator *(const complexvec_t &a, const complexvec_t &b)
{ return complexvec_t{a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re}; }
mathfn complexvec_t operator *(const complexvec_t &a, const realvec_t &b)
{ return complexvec_t{a.re*b, a.im*b}; }
mathfn complexvec_t operator *(const realvec_t &a, const complexvec_t &b)
{ return complexvec_t{a*b.re, a*b.im}; }
mathfn complexvec_t operator *(const complexvec_t &a, const real_t &b)
{ return complexvec_t{a.re*b, a.im*b}; }
mathfn complexvec_t operator *(const real_t &a, const complexvec_t &b)
{ return complexvec_t{a*b.re, a*b.im}; }
mathfn complexvec_t operator *(const complexvec_t &a, const complex_t &b)
{ return complexvec_t{a.re*b.real() - a.im*b.imag(), a.re*b.imag() + a.im*b.real()}; }
mathfn complexvec_t operator *(const complex_t &a, const complexvec_t &b)
{ return complexvec_t{a.real()*b.re - a.imag()*b.im, a.real()*b.im + a.imag()*b.re}; }
mathfn complexvec_t operator *(const realvec_t &a, const complex_t &b)
{ return complexvec_t{a*b.real() , a*b.imag()}; }
mathfn complexvec_t operator *(const complex_t &a, const realvec_t &b)
{ return complexvec_t{a.real()*b, a.imag()*b}; }

mathfn complexvec_t operator /(const complexvec_t &a, const complexvec_t &b)
{ realvec_t inv_abs2_b = 1/(b.re*b.re + b.im*b.im);
  return complexvec_t{ (a.re*b.re + a.im*b.im)*inv_abs2_b,
                       (a.im*b.re - a.re*b.im)*inv_abs2_b };}
mathfn complexvec_t operator /(const realvec_t &a, const complexvec_t &b)
{ realvec_t a_over_abs2_b = a/(b.re*b.re + b.im*b.im);
  return complexvec_t{ b.re*a_over_abs2_b, -b.im*a_over_abs2_b };}
mathfn complexvec_t operator /(const complexvec_t &a, const realvec_t &b)
{ realvec_t inv_b = 1/b; return complexvec_t{ a.re*inv_b, a.im*inv_b };}
mathfn complexvec_t operator /(const real_t &a, const complexvec_t &b)
{ realvec_t a_over_abs2_b = a/(b.re*b.re + b.im*b.im);
  return complexvec_t{ b.re*a_over_abs2_b, -b.im*a_over_abs2_b };}
mathfn complexvec_t operator /(const complexvec_t &a, const real_t &b)
{ real_t inv_b = 1/b; return complexvec_t{ a.re*inv_b, a.im*inv_b };}
mathfn complexvec_t operator /(const complexvec_t &a, const complex_t &b)
{ real_t inv_abs2_b = 1/(b.real()*b.real() + b.imag()*b.imag());
  return complexvec_t{ (a.re*b.real() + a.im*b.imag())*inv_abs2_b,
                       (a.im*b.real() - a.re*b.imag())*inv_abs2_b };}
mathfn complexvec_t operator /(const complex_t &a, const complexvec_t &b)
{ realvec_t inv_abs2_b = 1/(b.re*b.re + b.im*b.im);
  return complexvec_t{ (a.real()*b.re + a.imag()*b.im)*inv_abs2_b,
                       (a.imag()*b.re - a.real()*b.im)*inv_abs2_b };}
mathfn complexvec_t operator /(const realvec_t &a, const complex_t &b)
{ real_t inv_abs2_b = 1/(b.real()*b.real() + b.imag()*b.imag());
  return complexvec_t{ a*b.real()*inv_abs2_b, a*(-b.imag())*inv_abs2_b };}
mathfn complexvec_t operator /(const complex_t &a, const realvec_t &b)
{ realvec_t inv_b = 1/b;
  return complexvec_t{ (a.real()*inv_b, a.imag()*inv_b) };}

mathfn realvec_t SecDecInternalRealPart(const complexvec_t &a) { return a.re; }
mathfn realvec_t SecDecInternalImagPart(const complexvec_t &a) { return a.im; }

#define DEF_COMPLEX_FUNCTION(fname, fn) \
    static inline complexvec_t fname(const complexvec_t &a) { \
        complex_t c0 = fn(complex_t{a.re.x[0], a.im.x[0]}); \
        complex_t c1 = fn(complex_t{a.re.x[1], a.im.x[1]}); \
        complex_t c2 = fn(complex_t{a.re.x[2], a.im.x[2]}); \
        complex_t c3 = fn(complex_t{a.re.x[3], a.im.x[3]}); \
        return complexvec_t{ \
            {{c0.real(), c1.real(), c2.real(), c3.real()}}, \
            {{c0.imag(), c1.imag(), c2.imag(), c3.imag()}} \
        }; \
    } \
    static inline complexvec_t fname(const realvec_t &a) { \
        complex_t c0 = fn(a.x[0]); \
        complex_t c1 = fn(a.x[1]); \
        complex_t c2 = fn(a.x[2]); \
        complex_t c3 = fn(a.x[3]); \
        return complexvec_t{ \
            {{c0.real(), c1.real(), c2.real(), c3.real()}}, \
            {{c0.imag(), c1.imag(), c2.imag(), c3.imag()}} \
        }; \
    }

mathfn complexvec_t SecDecInternalI(const realvec_t &a)
{ return complexvec_t{REALVEC_ZERO, a}; }
mathfn complexvec_t SecDecInternalI(const complexvec_t &a)
{ return complexvec_t{-a.im, a.re}; }

mathfn complex_t componentsum(const complexvec_t &a)
{ return complex_t{ componentsum(a.re), componentsum(a.im) }; }

// Result vectors

#if SECDEC_RESULT_IS_COMPLEX

    typedef complex_t result_t;
    typedef complexvec_t resultvec_t;
    #define RESULTVEC_ZERO COMPLEXVEC_ZERO

    static inline complex_t SecDecInternalLog(const real_t x)
    { return (x >= 0) ? std::log(x) : complex_t{std::log(-x), -M_PI}; };
    static inline complex_t SecDecInternalLog(const complex_t x)
    { return (x.imag() == 0) ? SecDecInternalLog(x.real()) : std::log(x); };

    DEF_COMPLEX_FUNCTION(SecDecInternalLog, SecDecInternalLog)

#else

    typedef real_t result_t;
    typedef realvec_t resultvec_t;
    #define RESULTVEC_ZERO REALVEC_ZERO

    static inline real_t SecDecInternalLog(real_t x)
    { return std::log(x); };

#endif

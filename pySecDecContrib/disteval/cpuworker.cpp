#define __STDC_FORMAT_MACROS
#include <assert.h>
#include <dlfcn.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <vector>

#include <ginac/ginac.h>
#include <ginac/parser.h>
#include <streambuf>
#include <istream>
#include <fstream>
#include <sstream>

#ifdef unlikely
    #undef unlikely
#endif

#if __GNUC__
    #define unlikely(x) __builtin_expect((x), 0)
#else
    #define unlikely(x) (x)
#endif

typedef double real_t;
typedef struct { double re, im; } complex_t;

#define MAXPATH 4095
#define MAXNAME 255
#define MAXDIM 32

typedef int (*IntegrateF)(
    void * presult,
    const uint64_t lattice,
    const uint64_t index1,
    const uint64_t index2,
    const uint64_t * genvec,
    const real_t * shift,
    const real_t * realp,
    const complex_t * complexp,
    const real_t * deformp
);

typedef void (*MaxdeformpF)(
    real_t * deformp,
    const uint64_t lattice,
    const uint64_t index1,
    const uint64_t index2,
    const uint64_t * genvec,
    const real_t * shift,
    const real_t * realp,
    const complex_t * complexp
);

typedef int (*FpolycheckF)(
    const uint64_t lattice,
    const uint64_t index1,
    const uint64_t index2,
    const uint64_t * genvec,
    const real_t * shift,
    const real_t * realp,
    const complex_t * complexp,
    const real_t * deformp
);

struct Family {
    uint64_t dimension;
    real_t realp[MAXDIM];
    complex_t complexp[MAXDIM];
    bool complex_result;
    void* so_handle;
    char name[MAXNAME + 1];
};

struct Kernel {
    uint64_t familyidx;
    IntegrateF fn_integrate;
    MaxdeformpF fn_maxdeformp;
    FpolycheckF fn_fpolycheck;
    char name[MAXNAME + 1];
};

// Command structures

struct StartCmd {
    char dirname[MAXPATH + 1];
};

struct FamilyCmd {
    uint64_t index;
    char name[MAXNAME + 1];
    uint64_t dimension;
    real_t realp[MAXDIM];
    complex_t complexp[MAXDIM];
    bool complex_result;
};

struct KernelCmd {
    uint64_t index;
    uint64_t familyidx;
    char name[MAXNAME + 1];
};

struct PresampleCmd {
    uint64_t kernelidx;
    uint64_t ndeformp;
    uint64_t lattice;
    uint64_t genvec[MAXDIM];
    real_t shift[MAXDIM];
};

struct IntegrateCmd {
    uint64_t kernelidx;
    uint64_t lattice;
    uint64_t i1;
    uint64_t i2;
    uint64_t genvec[MAXDIM];
    real_t shift[MAXDIM];
    real_t deformp[MAXDIM];
};

// Global data
static char workername[MAXNAME];
static std::vector<Family> families;
static std::vector<Kernel> kernels;
static char *input_line = NULL;
static char *input_p = NULL;
static size_t input_linesize = 0;

#define input_getchar() (*input_p++)
#define input_peekchar() (*input_p)

static double
timestamp()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec*1e-9;
}

static double
cmd_start(uint64_t token, StartCmd &c)
{
    int r = chdir(c.dirname);
    if (r == 0) {
        printf("@[%" PRIu64 ",\"%s\",null]\n", token, workername);
    } else {
        printf("@[%" PRIu64 ",null,\"failed to chdir '%s': %d\"]\n", token, c.dirname, r);
    }
    return 0;
}

static double
cmd_family(uint64_t token, FamilyCmd &c)
{
    assert(c.index == families.size());
    char buf[MAXNAME+8];
    snprintf(buf, sizeof(buf), "./%s.so", c.name);
    void *so_handle = dlopen(buf, RTLD_LAZY | RTLD_LOCAL);
    if (so_handle == NULL) {
        printf("@[%" PRIu64 ",null,\"failed to open '%s': %s\"]\n", token, buf, strerror(errno));
        return 0;
    }
    Family fam = {};
    fam.dimension = c.dimension;
    memcpy(fam.realp, c.realp, sizeof(fam.realp));
    memcpy(fam.complexp, c.complexp, sizeof(fam.complexp));
    fam.complex_result = c.complex_result;
    fam.so_handle = so_handle;
    memcpy(fam.name, c.name, sizeof(fam.name));
    families.push_back(fam);
    printf("@[%" PRIu64 ",null,null]\n", token);
    return 0;
}

static double
cmd_change_family_parameters(uint64_t token, FamilyCmd &c)
{
    assert(c.index < families.size());
    Family &fam = families[c.index];
    memcpy(fam.realp, c.realp, sizeof(fam.realp));
    memcpy(fam.complexp, c.complexp, sizeof(fam.complexp));
    printf("@[%" PRIu64 ",null,null]\n", token);
    return 0;
}

static double
cmd_kernel(uint64_t token, KernelCmd &c)
{
    assert(c.familyidx < families.size());
    assert(c.index == kernels.size());
    const Family &fam = families[c.familyidx];
    char buf[2*MAXNAME+18];
    Kernel ker = {};
    ker.familyidx = c.familyidx;
    snprintf(buf, sizeof(buf), "%s__%s", fam.name, c.name);
    ker.fn_integrate = (IntegrateF)dlsym(fam.so_handle, buf);
    if (ker.fn_integrate == NULL) {
        printf("@[%" PRIu64 ",null,\"function not found: %s\"]\n", token, buf);
        return 0;
    }
    snprintf(buf, sizeof(buf), "%s__%s__maxdeformp", fam.name, c.name);
    ker.fn_maxdeformp = (MaxdeformpF)dlsym(fam.so_handle, buf);
    snprintf(buf, sizeof(buf), "%s__%s__fpolycheck", fam.name, c.name);
    ker.fn_fpolycheck = (FpolycheckF)dlsym(fam.so_handle, buf);
    memcpy(ker.name, c.name, sizeof(ker.name));
    kernels.push_back(ker);
    printf("@[%" PRIu64 ",null,null]\n", token);
    return 0;
}

static double
cmd_presample(uint64_t token, PresampleCmd &c)
{
    if (unlikely(c.kernelidx >= kernels.size())) {
        printf("@[%" PRIu64 ",null,\"kernel %" PRIu64 " was not loaded\"]\n", token, c.kernelidx);
        return 0;
    }
    const Kernel &ker = kernels[c.kernelidx];
    const Family &fam = families[ker.familyidx];
    if (unlikely(c.ndeformp == 0)) {
        printf("@[%" PRIu64 ",[],null]\n", token);
        return 0;
    }
    if (unlikely(ker.fn_maxdeformp == NULL)) {
        printf("@[%" PRIu64 ",null,\"kernel %" PRIu64 " has no *__maxdefomp function\"]\n", token, c.kernelidx);
        return 0;
    }
    if (unlikely(ker.fn_fpolycheck == NULL)) {
        printf("@[%" PRIu64 ",null,\"kernel %" PRIu64 " has no *__fpolycheck function\"]\n", token, c.kernelidx);
        return 0;
    }
    double deformp[MAXDIM] = {};
    double t1 = timestamp();
    ker.fn_maxdeformp(deformp,
        c.lattice, 0, c.lattice, c.genvec, c.shift,
        fam.realp, fam.complexp);
    for (;;) {
        int r =
            ker.fn_fpolycheck(
                c.lattice, 0, c.lattice, c.genvec, c.shift,
                fam.realp, fam.complexp, deformp);
        if (r == 0) break;
        for (uint64_t i = 0; i < c.ndeformp; i++) deformp[i] *= 0.9;
    }
    double t2 = timestamp();
    printf("@[%" PRIu64 ",[", token);
    for (uint64_t i = 0; i < c.ndeformp; i++) {
        if (i != 0) putchar(',');
        printf("%.16e", deformp[i]);
    }
    printf("],null]\n");
    return t2-t1;
}

static double
cmd_integrate(uint64_t token, IntegrateCmd &c)
{
    if (unlikely(c.kernelidx >= kernels.size())) {
        printf("@[%" PRIu64 ",null,\"kernel %" PRIu64 " was not loaded\"]\n", token, c.kernelidx);
        return 0;
    }
    const Kernel &ker = kernels[c.kernelidx];
    const Family &fam = families[ker.familyidx];
    complex_t result = {};
    double t1 = timestamp();
    int r = ker.fn_integrate(&result,
        c.lattice, c.i1, c.i2, c.genvec, c.shift,
        fam.realp, fam.complexp, c.deformp);
    double t2 = timestamp();
    if (unlikely((isnan(result.re) || isnan(result.im)) ^ (r != 0))) {
        printf("@[%" PRIu64 ",[[NaN,NaN],%" PRIu64 ",%.4e],\"NaN != sign check error %d in %s.%s\"]\n", token, c.i2-c.i1, t2-t1, r, fam.name, ker.name);
    } else if (isnan(result.re) || isnan(result.im)) {
        printf("@[%" PRIu64 ",[[NaN,NaN],%" PRIu64 ",%.4e],null]\n", token, c.i2-c.i1, t2-t1);
    } else {
        printf("@[%" PRIu64 ",[[%.16e,%.16e],%" PRIu64 ",%.4e],null]\n", token, result.re, result.im, c.i2-c.i1, t2-t1);
    }
    return t2-t1;
}

static void
parse_fail()
{
    fprintf(stderr, "%s] input parsing failed:\n", workername);
    fprintf(stderr, "%s", input_line);
    for (char *p = input_line + 1; p < input_p; p++)
        putc('-', stderr);
    fprintf(stderr, "^\n");
    exit(1);
}

static void
match_c(char c)
{
    if (unlikely(c != input_getchar())) {parse_fail();}
}

static void
match_str(const char *s)
{
    for (; *s; s++) {
        int c = input_getchar();
        if (unlikely(c != *s)) {parse_fail();}
    }
}

static bool
parse_bool()
{
    int c = input_getchar();
    if (c == 't') { match_str("rue"); return true; }
    if (c == 'f') { match_str("alse"); return false; }
    parse_fail();
    exit(1);
}

static uint64_t
parse_uint()
{
    char *end = NULL;
    long long x = strtoll(input_p, &end, 10);
    if (unlikely(input_p == end)) parse_fail();
    input_p = end;
    return (uint64_t)x;
}

static int64_t
parse_int()
{
    char *end = NULL;
    long long x = strtoll(input_p, &end, 10);
    if (unlikely(input_p == end)) parse_fail();
    input_p = end;
    return (int64_t)x;
}

static real_t
parse_real()
{
    char *end = NULL;
    real_t x = strtod(input_p, &end);
    if (unlikely(input_p == end)) parse_fail();
    input_p = end;
    return x;
}

static complex_t
parse_complex()
{
    match_c('[');
    real_t re = parse_real();
    match_c(',');
    real_t im = parse_real();
    match_c(']');
    return complex_t{re, im};
}

#define define_parse_X_array(name, type, parse_X) \
    static void \
    name(type *ar, size_t maxn) \
    { \
        match_c('['); \
        int c = input_peekchar(); \
        if (c == ']') { input_getchar(); return; } \
        for (size_t i = 0;; i++) { \
            if (unlikely(i >= maxn)) parse_fail(); \
            ar[i] = parse_X(); \
            int c = input_getchar(); \
            if (c == ']') break; \
            if (unlikely(c != ',')) parse_fail(); \
        } \
    }
define_parse_X_array(parse_uint_array, uint64_t, parse_uint)
define_parse_X_array(parse_real_array, double, parse_real)
define_parse_X_array(parse_complex_array, complex_t, parse_complex)

static void
parse_str(char *str, size_t maxn)
{
    match_c('"');
    for (size_t i = 0; ; i++) {
        if (unlikely(i >= maxn)) parse_fail();
        int c = input_getchar();
        if (unlikely(c == '\\')) {
            int cc = input_getchar();
            switch (cc) {
                case '"': case '\\': case '/': str[i] = cc; break;
                case 'b': str[i] = '\b'; break;
                case 'f': str[i] = '\f'; break;
                case 'n': str[i] = '\n'; break;
                case 'r': str[i] = '\r'; break;
                case 't': str[i] = '\t'; break;
                default: parse_fail();
            }
        } else if (unlikely(c == '"')) {
            str[i] = 0;
            break;
        } else {
            str[i] = c;
        }
    }
}

// GiNaC-related code

struct FixedStreamBuf : public std::streambuf {
    FixedStreamBuf(char* s, size_t n) { setg(s, s, s + n); }
};

static GiNaC::ex
ginac_read_string(GiNaC::parser &reader, char *str, size_t size)
{
    FixedStreamBuf buf(str, size);
    std::istream i(&buf);
    return reader(i);
}

static GiNaC::ex
ginac_read_string(GiNaC::parser &reader, char *str)
{
    return ginac_read_string(reader, str, strlen(str));
}

template <typename F> void
term_iter(const GiNaC::ex &e, F yield)
{
    if (GiNaC::is_a<GiNaC::add>(e)) {
        for (const auto &t : e) {
            yield(t);
        }
    } else {
        yield(e);
    }
}

template <typename F> void
factor_iter(const GiNaC::ex &e, F yield)
{
    if (GiNaC::is_a<GiNaC::mul>(e)) {
        for (const auto &f : e) {
            if (GiNaC::is_a<GiNaC::power>(f)) {
                yield(f.op(0), GiNaC::ex_to<GiNaC::numeric>(f.op(1)).to_int());
            } else {
                yield(f, 1);
            }
        }
    } else {
        if (GiNaC::is_a<GiNaC::power>(e)) {
            yield(e.op(0), GiNaC::ex_to<GiNaC::numeric>(e.op(1)).to_int());
        } else {
            yield(e, 1);
        }
    }
}

static std::map<std::vector<int>, GiNaC::ex>
ginac_bracket(const GiNaC::ex &expr, const GiNaC::exvector &X)
{
    std::map<std::vector<int>, GiNaC::ex> result;
    std::map<GiNaC::ex, int, GiNaC::ex_is_less> x2id;
    for (unsigned i = 0; i < X.size(); i++) {
        x2id[X[i]] = i;
    }
    term_iter(expr.expand(), [&](const GiNaC::ex &term) {
        std::vector<int> stemidx(X.size());
        GiNaC::exvector coef;
        factor_iter(term, [&](const GiNaC::ex &factor, int power) {
            auto it = x2id.find(factor);
            if (it != x2id.end()) {
                stemidx[it->second] = power;
            } else {
                if (power == 1) {
                    coef.push_back(factor);
                } else {
                    coef.push_back(pow(factor, power));
                }
            }
        });
        result[stemidx] += GiNaC::mul(coef);
    });
    return result;
}

static double
parse_cmd_evalf(uint64_t token)
{
    char filename[MAXPATH];
    parse_str(filename, sizeof(filename));
    GiNaC::parser reader;
    double t1 = timestamp();
    std::ifstream inf(filename);
    if (!inf) {
        printf("@[%" PRIu64 ",null,\"failed to open '%s'\"]\n", token, filename);
        exit(1);
    }
    GiNaC::ex expr = reader(inf);
    GiNaC::exmap table;
    match_str(",{");
    char varname[MAXNAME];
    if (input_peekchar() != '}') {
        char value[MAXPATH];
        for (;;) {
            parse_str(varname, sizeof(varname));
            match_c(':');
            parse_str(value, sizeof(value));
            table[ginac_read_string(reader, varname)] = ginac_read_string(reader, value);
            if (input_peekchar() != ',') break;
            input_getchar();
        }
    }
    match_str("},[");
    expr = expr.subs(table);
    GiNaC::exvector varlist;
    for (;;) {
        match_c('[');
        parse_str(varname, sizeof(varname));
        match_c(',');
        int64_t order = parse_int();
        match_c(']');
        auto x = ginac_read_string(reader, varname);
        varlist.push_back(x);
        expr = GiNaC::series_to_poly(expr.series(x, order+1));
        if (input_peekchar() != ',') break;
        input_getchar();
    }
    match_str("]]]\n");
    auto br = ginac_bracket(expr.expand(), varlist);
    std::map<std::vector<int>, complex_t> brc;
    for (auto &&kv : br) {
        GiNaC::ex val = kv.second.evalf();
        GiNaC::ex val_re = val.real_part();
        GiNaC::ex val_im = val.imag_part();
        if (GiNaC::is_a<GiNaC::numeric>(val_re) && GiNaC::is_a<GiNaC::numeric>(val_im)) {
            double re = GiNaC::ex_to<GiNaC::numeric>(val_re).to_double();
            double im = GiNaC::ex_to<GiNaC::numeric>(val_im).to_double();
            brc[kv.first] = complex_t{re, im};
        } else {
            printf("@[%" PRIu64 ",null,\"the coefficient is not numeric after substitution\"]\n", token);
            return 0;
        }
    }
    double t2 = timestamp();
    printf("@[%" PRIu64 ",[", token);
    bool first = true;
    for (auto &&kv : brc) {
        if (first) { first = false; } else { putchar(','); }
        printf("[[");
        bool first2 = true;
        for (auto &&i : kv.first) {
            if (first2) { first2 = false; } else { putchar(','); }
            printf("%d", i);
        }
        printf("],[%.16e,%.16e]]", kv.second.re, kv.second.im);
    }
    printf("],null]\n");
    return t2 - t1;
}

// Main RPC cycle

static double
handle_one_command()
{
    match_c('[');
    uint64_t token = parse_uint();
    match_c(','); match_c('"');
    int c = input_getchar();
    if (c == 'i') {
        IntegrateCmd c = {};
        match_str("ntegrate\",[");
        c.kernelidx = parse_uint();
        match_c(',');
        c.lattice = parse_uint();
        match_c(',');
        c.i1 = parse_uint();
        match_c(',');
        c.i2 = parse_uint();
        match_c(',');
        parse_uint_array(c.genvec, MAXDIM);
        match_c(',');
        parse_real_array(c.shift, MAXDIM);
        match_c(',');
        parse_real_array(c.deformp, MAXDIM);
        match_str("]]\n");
        return cmd_integrate(token, c);
    }
    if (c == 'm') {
        PresampleCmd c = {};
        match_str("axdeformp\",[");
        c.kernelidx = parse_uint();
        match_c(',');
        c.ndeformp = parse_uint();
        match_c(',');
        c.lattice = parse_uint();
        match_c(',');
        parse_uint_array(c.genvec, MAXDIM);
        match_c(',');
        parse_real_array(c.shift, MAXDIM);
        match_str("]]\n");
        return cmd_presample(token, c);
    }
    if (c == 'p') {
        match_str("ing\",[]]\n");
        printf("@[%" PRIu64 ",null,null]\n", token);
        return 0;
    }
    if (c == 'f') {
        FamilyCmd c = {};
        match_str("amily\",[");
        c.index = parse_uint();
        match_c(',');
        parse_str(c.name, sizeof(c.name));
        match_c(',');
        c.dimension = parse_uint();
        match_c(',');
        parse_real_array(c.realp, MAXDIM);
        match_c(',');
        parse_complex_array(c.complexp, MAXDIM);
        match_c(',');
        c.complex_result = parse_bool();
        match_str("]]\n");
        return cmd_family(token, c);
    }
    if (c == 'c') {
        FamilyCmd c = {};
        match_str("hangefamily\",[");
        c.index = parse_uint();
        match_c(',');
        parse_real_array(c.realp, MAXDIM);
        match_c(',');
        parse_complex_array(c.complexp, MAXDIM);
        match_str("]]\n");
        return cmd_change_family_parameters(token, c);
    }
    if (c == 'k') {
        KernelCmd c = {};
        match_str("ernel\",[");
        c.index = parse_uint();
        match_c(',');
        c.familyidx = parse_uint();
        match_c(',');
        parse_str(c.name, sizeof(c.name));
        match_str("]]\n");
        return cmd_kernel(token, c);
    }
    if (c == 's') {
        StartCmd c = {};
        match_str("tart\",[");
        parse_str(c.dirname, sizeof(c.dirname));
        match_str("]]\n");
        return cmd_start(token, c);
    }
    if (c == 'e') {
        match_str("valf\",[");
        return parse_cmd_evalf(token);
    }
    parse_fail();
    return 0;
}

static void
fill_workername()
{
    char host[MAXNAME] = {};
    gethostname(host, sizeof(host));
    long pid = getpid();
    snprintf(workername, sizeof(workername), "%s:%ld", host, pid);
}

int main() {
    fill_workername();
    setvbuf(stdin, NULL, _IOFBF, 1024*1024);
    setvbuf(stdout, NULL, _IOLBF, 1024*1024);
    setvbuf(stderr, NULL, _IOLBF, 1024*1024);
    double readt = 0;
    double workt = 0;
    double lastt = 0;
    double t1 = timestamp();
    bool quit = false;
    while (!quit) {
        lastt = timestamp();
        if (getline(&input_line, &input_linesize, stdin) < 0) break;
        readt += timestamp() - lastt;
        input_p = input_line;
        workt += handle_one_command();
    }
    double t2 = timestamp();
    if (0) {
        fprintf(stderr, "%s] Done in %.3gs: %.3g%% useful time, %.3g%% read time; work ended %.3gs ago\n",
                workername, lastt-t1, 100*workt/(lastt-t1), 100*readt/(lastt-t1), t2-lastt);
    }
}

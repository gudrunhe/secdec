#include <assert.h>
#include <dlfcn.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <vector>

//#define unlikely(x) (x)
#define unlikely(x) __builtin_expect((x), 0)

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
        printf("@[%zu,\"%s\",null]\n", token, workername);
    } else {
        printf("@[%zu,null,\"failed to chdir %s: %d\"]\n", token, c.dirname, r);
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
        printf("@[%zu,null,\"failed to open %s: %s\"]\n", token, buf, strerror(errno));
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
    printf("@[%zu,null,null]\n", token);
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
        printf("@[%zu,null,\"function not found: %s\"]\n", token, buf);
        return 0;
    }
    snprintf(buf, sizeof(buf), "%s__%s__maxdeformp", fam.name, c.name);
    ker.fn_maxdeformp = (MaxdeformpF)dlsym(fam.so_handle, buf);
    snprintf(buf, sizeof(buf), "%s__%s__fpolycheck", fam.name, c.name);
    ker.fn_fpolycheck = (FpolycheckF)dlsym(fam.so_handle, buf);
    memcpy(ker.name, c.name, sizeof(ker.name));
    kernels.push_back(ker);
    printf("@[%zu,null,null]\n", token);
    return 0;
}

static double
cmd_presample(uint64_t token, PresampleCmd &c)
{
    if (unlikely(c.kernelidx >= kernels.size())) {
        printf("@[%zu,null,\"kernel %zu was not loaded\"]\n", token, c.kernelidx);
        return 0;
    }
    const Kernel &ker = kernels[c.kernelidx];
    const Family &fam = families[ker.familyidx];
    if (unlikely(c.ndeformp == 0)) {
        printf("@[%zu,[],null]\n", token);
        return 0;
    }
    if (unlikely(ker.fn_maxdeformp == NULL)) {
        printf("@[%zu,null,\"kernel %zu has no *__maxdefomp function\"]\n", token, c.kernelidx);
        return 0;
    }
    if (unlikely(ker.fn_fpolycheck == NULL)) {
        printf("@[%zu,null,\"kernel %zu has no *__fpolycheck function\"]\n", token, c.kernelidx);
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
        for (int i = 0; i < c.ndeformp; i++) deformp[i] *= 0.9;
    }
    double t2 = timestamp();
    printf("@[%zu,[", token);
    for (int i = 0; i < c.ndeformp; i++) {
        if (i != 0) putchar(',');
        printf("%.16e", deformp[i]);
    }
    printf("],null]\n");
    return t2-t1;
}

static void
print_real(double x)
{
    if (isnan(x)) printf("NaN");
    else printf("%.16e", x);
}

static void
print_complex(complex_t x)
{
    putchar('[');
    print_real(x.re);
    putchar(',');
    print_real(x.im);
    putchar(']');
}

static double
cmd_integrate(uint64_t token, IntegrateCmd &c)
{
    if (unlikely(c.kernelidx >= kernels.size())) {
        printf("@[%zu,null,\"kernel %zu was not loaded\"]\n", token, c.kernelidx);
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
        printf("@[%zu,[[NaN,NaN],%zu,%.4e],\"NaN != sign check error %d in %s.%s\"]", token, c.i2-c.i1, t2-t1, r, fam.name, ker.name);
    }
    printf("@[%zu,[", token); print_complex(result); printf(",%zu,%.4e],null]\n", c.i2-c.i1, t2-t1);
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
    fflush(stderr);
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
    long x = strtol(input_p, &end, 10);
    if (unlikely(input_p == end)) parse_fail();
    input_p = end;
    return (uint64_t)x;
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
            break;
        } else {
            str[i] = c;
        }
    }
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
        printf("@[%zu,null,null]\n", token);
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
    setvbuf(stdout, NULL, _IOFBF, 1024*1024);
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
        fflush(stdout);
    }
    double t2 = timestamp();
    fprintf(stderr, "%s] Done in %.3gs: %.3g%% useful time, %.3g%% read time; work ended %.3gs ago\n",
            workername, lastt-t1, 100*workt/(lastt-t1), 100*readt/(lastt-t1), t2-lastt);
    fflush(stderr);
}

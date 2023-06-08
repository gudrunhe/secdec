#define __STDC_FORMAT_MACROS
#include <assert.h>
#include <dlfcn.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <vector>

#include "minicuda.h"

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
typedef struct { real_t re, im; } complex_t;

#define MAXPATH 4095
#define MAXNAME 255
#define MAXDIM 32
#define NTHREADS 1
#define MAXQUEUE 16

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
    CUmodule cuda_module;
    char name[MAXNAME + 1];
};

struct Kernel {
    uint64_t familyidx;
    IntegrateF fn_integrate;
    MaxdeformpF fn_maxdeformp;
    FpolycheckF fn_fpolycheck;
    CUfunction cuda_fn_integrate;
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
    uint64_t token;
    uint64_t kernelidx;
    uint64_t lattice;
    uint64_t i1;
    uint64_t i2;
    uint64_t genvec[MAXDIM];
    real_t shift[MAXDIM];
    real_t deformp[MAXDIM];
};

struct CudaParameterData {
    uint64_t genvec[MAXDIM];
    real_t shift[MAXDIM];
    real_t realp[MAXDIM];
    complex_t complexp[MAXDIM];
    real_t deformp[MAXDIM];
};

// Global state

struct PerThreadState {
    pthread_t thread;
    CUstream stream;
    CUdeviceptr buffer_d;
    size_t buffer_size;
    CUdeviceptr params_d;
    CudaParameterData *params;
    complex_t *result;
    double useful_time;
};

static struct GlobalState {
    char workername[MAXNAME];
    std::vector<Family> families;
    std::vector<Kernel> kernels;
    char *input_line = NULL;
    char *input_p = NULL;
    size_t input_linesize = 0;
    double useful_time = 0;
    struct GlobalCudaState {
        CUdevice device;
        CUcontext context;
        CUmodule builtin_module;
        CUfunction fn_sum_d_b128_x1024;
        CUfunction fn_sum_c_b128_x1024;
    } cuda;
    struct IntegrateCmdQueue {
        IntegrateCmd cmd[MAXQUEUE];
        int head = 0;
        int tail = 0;
        pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
        pthread_cond_t cond_ins = PTHREAD_COND_INITIALIZER;
        pthread_cond_t cond_rem = PTHREAD_COND_INITIALIZER;
    } queue;
    PerThreadState threads[NTHREADS];
} G;

#define input_getchar() (*G.input_p++)
#define input_peekchar() (*G.input_p)

static double
timestamp()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec*1e-9;
}

// Work queue

static void
submit_integrate_cmd(const IntegrateCmd &cmd)
{
    pthread_mutex_lock(&G.queue.lock);
    while (((G.queue.head + 1) % MAXQUEUE) == G.queue.tail) {
        pthread_cond_wait(&G.queue.cond_rem, &G.queue.lock);
    }
    G.queue.cmd[G.queue.head] = cmd;
    G.queue.head = (G.queue.head + 1) % MAXQUEUE;
    pthread_cond_signal(&G.queue.cond_ins);
    pthread_mutex_unlock(&G.queue.lock);
}

static void
obtain_integrate_cmd(IntegrateCmd &cmd)
{
    pthread_mutex_lock(&G.queue.lock);
    while (G.queue.head == G.queue.tail) {
        pthread_cond_wait(&G.queue.cond_ins, &G.queue.lock);
    }
    cmd = G.queue.cmd[G.queue.tail];
    G.queue.tail = (G.queue.tail + 1) % MAXQUEUE;
    pthread_cond_signal(&G.queue.cond_rem);
    pthread_mutex_unlock(&G.queue.lock);
}

// CUDA

#define CU(fn, ...) { \
        CUresult curesult = fn(__VA_ARGS__); \
        if (unlikely(curesult != 0)) cuda_fail(#fn, curesult); \
    }

static void
cuda_fail(const char *what, CUresult code)
{
    const char *n = NULL, *s = NULL;
    cuGetErrorName(code, &n);
    cuGetErrorString(code, &s);
    fprintf(stderr, "%s] %s failed with code %d (%s): %s\n", G.workername, what, code, n, s);
    exit(1);
}

// Commands

static void
cmd_start(uint64_t token, StartCmd &c)
{
    int r = chdir(c.dirname);
    if (r != 0) {
        printf("@[%" PRIu64 ",null,\"failed to chdir '%s': %d\"]\n", token, c.dirname, r);
        return;
    }
    CU(cuModuleLoad, &G.cuda.builtin_module, "./builtin.fatbin");
    CU(cuModuleGetFunction, &G.cuda.fn_sum_d_b128_x1024, G.cuda.builtin_module, "sum_d_b128_x1024");
    CU(cuModuleGetFunction, &G.cuda.fn_sum_c_b128_x1024, G.cuda.builtin_module, "sum_c_b128_x1024");
    printf("@[%" PRIu64 ",\"%s\",null]\n", token, G.workername);
}

static void
cmd_family(uint64_t token, FamilyCmd &c)
{
    assert(c.index == G.families.size());
    Family fam = {};
    char buf[MAXNAME+16];
    snprintf(buf, sizeof(buf), "./%s.so", c.name);
    fam.so_handle = dlopen(buf, RTLD_LAZY | RTLD_LOCAL);
    if (fam.so_handle == NULL) {
        printf("@[%" PRIu64 ",null,\"failed to open '%s': %s\"]\n", token, buf, strerror(errno));
        return;
    }
    snprintf(buf, sizeof(buf), "./%s.fatbin", c.name);
    if (cuModuleLoad(&fam.cuda_module, buf) != 0) {
        printf("@[%" PRIu64 ",null,\"failed to open '%s'\"]\n", token, buf);
        return;
    }
    fam.dimension = c.dimension;
    memcpy(fam.realp, c.realp, sizeof(fam.realp));
    memcpy(fam.complexp, c.complexp, sizeof(fam.complexp));
    fam.complex_result = c.complex_result;
    memcpy(fam.name, c.name, sizeof(fam.name));
    G.families.push_back(fam);
    printf("@[%" PRIu64 ",null,null]\n", token);
}

static void
cmd_change_family_parameters(uint64_t token, FamilyCmd &c)
{
    assert(c.index < G.families.size());
    Family &fam = G.families[c.index];
    memcpy(fam.realp, c.realp, sizeof(fam.realp));
    memcpy(fam.complexp, c.complexp, sizeof(fam.complexp));
    printf("@[%" PRIu64 ",null,null]\n", token);
}

static void
cmd_kernel(uint64_t token, KernelCmd &c)
{
    assert(c.familyidx < G.families.size());
    assert(c.index == G.kernels.size());
    const Family &fam = G.families[c.familyidx];
    char buf[2*MAXNAME+18];
    Kernel ker = {};
    ker.familyidx = c.familyidx;
    snprintf(buf, sizeof(buf), "%s__%s", fam.name, c.name);
    ker.fn_integrate = (IntegrateF)dlsym(fam.so_handle, buf);
    if (ker.fn_integrate == NULL) {
        printf("@[%" PRIu64 ",null,\"function not found: %s\"]\n", token, buf);
        return;
    }
    if (cuModuleGetFunction(&ker.cuda_fn_integrate, fam.cuda_module, buf) != 0) {
        printf("@[%" PRIu64 ",null,\"CUDA function not found: %s\"]\n", token, buf);
        return;
    }
    snprintf(buf, sizeof(buf), "%s__%s__maxdeformp", fam.name, c.name);
    ker.fn_maxdeformp = (MaxdeformpF)dlsym(fam.so_handle, buf);
    snprintf(buf, sizeof(buf), "%s__%s__fpolycheck", fam.name, c.name);
    ker.fn_fpolycheck = (FpolycheckF)dlsym(fam.so_handle, buf);
    memcpy(ker.name, c.name, sizeof(ker.name));
    G.kernels.push_back(ker);
    printf("@[%" PRIu64 ",null,null]\n", token);
}

static void
cmd_presample(uint64_t token, PresampleCmd &c)
{
    if (unlikely(c.kernelidx >= G.kernels.size())) {
        printf("@[%" PRIu64 ",null,\"kernel %" PRIu64 " was not loaded\"]\n", token, c.kernelidx);
        return;
    }
    const Kernel &ker = G.kernels[c.kernelidx];
    const Family &fam = G.families[ker.familyidx];
    if (unlikely(c.ndeformp == 0)) {
        printf("@[%" PRIu64 ",[],null]\n", token);
        return;
    }
    if (unlikely(ker.fn_maxdeformp == NULL)) {
        printf("@[%" PRIu64 ",null,\"kernel %" PRIu64 " has no *__maxdefomp function\"]\n", token, c.kernelidx);
        return;
    }
    if (unlikely(ker.fn_fpolycheck == NULL)) {
        printf("@[%" PRIu64 ",null,\"kernel %" PRIu64 " has no *__fpolycheck function\"]\n", token, c.kernelidx);
        return;
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
        for (uint64_t i = 0; i < c.ndeformp; i++)
            deformp[i] *= 0.9;
    }
    double t2 = timestamp();
    printf("@[%" PRIu64 ",[", token);
    for (uint64_t i = 0; i < c.ndeformp; i++) {
        if (i != 0) putchar(',');
        printf("%.16e", deformp[i]);
    }
    printf("],null]\n");
    G.useful_time += t2-t1;
}

static void stupid_cuda_dummy(void *userdata)
{ (void)userdata; }

static void *
worker_thread(void *ps)
{
    PerThreadState &s = *(PerThreadState*)ps;
    CU(cuCtxSetCurrent, G.cuda.context);
    for (;;) {
        IntegrateCmd c;
        obtain_integrate_cmd(c);
        const Kernel &ker = G.kernels[c.kernelidx];
        const Family &fam = G.families[ker.familyidx];
        if (0) { // CPU path
            complex_t result = {};
            double t1 = timestamp();
            int r = ker.fn_integrate(&result,
                c.lattice, c.i1, c.i2, c.genvec, c.shift,
                fam.realp, fam.complexp, c.deformp);
            double t2 = timestamp();
            if (unlikely((isnan(result.re) || isnan(result.im)) ^ (r != 0))) {
                printf("@[%" PRIu64 ",[[NaN,NaN],%" PRIu64 ",%.4e],\"NaN != sign check error %d in %s.%s\"]\n", c.token, c.i2-c.i1, t2-t1, r, fam.name, ker.name);
            } else if (isnan(result.re) || isnan(result.im)) {
                printf("@[%" PRIu64 ",[[NaN,NaN],%" PRIu64 ",%.4e],null]\n", c.token, c.i2-c.i1, t2-t1);
            } else {
                printf("@[%" PRIu64 ",[[%.16e,%.16e],%" PRIu64 ",%.4e],null]\n", c.token, result.re, result.im, c.i2-c.i1, t2-t1);
            }
            s.useful_time += t2-t1;
        }
        if (1) { // CUDA path
            uint64_t threads = 128, pt_per_thread = 8;
            uint64_t blocksperbatch = fam.complex_result ? s.buffer_size/sizeof(complex_t) : s.buffer_size/sizeof(real_t);
            uint64_t ptperbatch = blocksperbatch * (threads*pt_per_thread);
            complex_t result = {0, 0};
            memcpy(s.params->genvec, c.genvec, sizeof(c.genvec));
            memcpy(s.params->shift, c.shift, sizeof(c.shift));
            memcpy(s.params->realp, fam.realp, sizeof(fam.realp));
            memcpy(s.params->complexp, fam.complexp, sizeof(fam.complexp));
            memcpy(s.params->deformp, c.deformp, sizeof(c.deformp));
            double t1 = timestamp();
            CU(cuMemcpyHtoDAsync, s.params_d, s.params, sizeof(CudaParameterData), s.stream);
            CUdeviceptr genvec_d = s.params_d + offsetof(CudaParameterData, genvec);
            CUdeviceptr shift_d = s.params_d + offsetof(CudaParameterData, shift);
            CUdeviceptr realp_d = s.params_d + offsetof(CudaParameterData, realp);
            CUdeviceptr complexp_d = s.params_d + offsetof(CudaParameterData, complexp);
            CUdeviceptr deformp_d = s.params_d + offsetof(CudaParameterData, deformp);
            for (uint64_t i1 = c.i1; i1 < c.i2; i1 += ptperbatch) {
                uint64_t i2 = i1 + ptperbatch < c.i2 ? i1 + ptperbatch : c.i2;
                uint64_t blocks = (i2 - i1 + threads*pt_per_thread - 1)/(threads*pt_per_thread);
                void *args[] = {&s.buffer_d, &c.lattice, &i1, &i2, &genvec_d, &shift_d, &realp_d, &complexp_d, &deformp_d, NULL };
                CU(cuLaunchKernel, ker.cuda_fn_integrate, blocks, 1, 1, threads, 1, 1, 0, s.stream, args, NULL);
                void *sum_args[] = {&s.buffer_d, &s.buffer_d, &blocks, NULL};
                CUfunction fn_sum = fam.complex_result ? G.cuda.fn_sum_c_b128_x1024 : G.cuda.fn_sum_d_b128_x1024;
                while (blocks > 1) {
                    uint64_t reduced = (blocks + 1024-1)/1024;
                    CU(cuLaunchKernel, fn_sum, reduced, 1, 1, 128, 1, 1, 0, s.stream, sum_args, NULL);
                    blocks = reduced;
                }
                s.result->re = 0;
                s.result->im = 0;
                CU(cuMemcpyDtoHAsync, s.result, s.buffer_d, fam.complex_result ? sizeof(complex_t) : sizeof(real_t), s.stream);
                // Without this CU_CTX_SCHED_BLOCKING_SYNC doesn't work,
                // and cuStreamSynchronize spins with 100% CPU usage.
                // With this, both CU_CTX_SCHED_BLOCKING_SYNC and
                // CU_CTX_SCHED_YIELD have the same result: 0% CPU usage
                // during cuStreamSynchronize.
                // It's not clear how cuLaunchHostFunc is related here
                // at all, and the whole thing is completely undocumented.
                CU(cuLaunchHostFunc, s.stream, stupid_cuda_dummy, NULL);
                CU(cuStreamSynchronize, s.stream);
                result.re += s.result->re;
                result.im += s.result->im;
            }
            double t2 = timestamp();
            if (isnan(result.re) || isnan(result.im)) {
                printf("@[%" PRIu64 ",[[NaN,NaN],%" PRIu64 ",%.4e],null]\n", c.token, c.i2-c.i1, t2-t1);
            } else {
                printf("@[%" PRIu64 ",[[%.16e,%.16e],%" PRIu64 ",%.4e],null]\n", c.token, result.re, result.im, c.i2-c.i1, t2-t1);
            }
            s.useful_time += t2-t1;
        }
    }
    return NULL;
}

static void
cmd_integrate(uint64_t token, IntegrateCmd &c)
{
    if (unlikely(c.kernelidx >= G.kernels.size())) {
        printf("@[%" PRIu64 ",null,\"kernel %" PRIu64 " was not loaded\"]\n", token, c.kernelidx);
        return;
    }
    submit_integrate_cmd(c);
}

// Initialization

static void
init(int devindex)
{
    CU(cuInit, 0);
    int ver = 0, ndev = 0;
    size_t memsize = 0;
    CU(cuDriverGetVersion, &ver);
    CU(cuDeviceGetCount, &ndev);
    CU(cuDeviceGet, &G.cuda.device, devindex);
    char buf[256];
    CU(cuDeviceGetName, buf, sizeof(buf), G.cuda.device);
    CU(cuDeviceGetName, buf, sizeof(buf), G.cuda.device);
    CU(cuDeviceTotalMem, &memsize, G.cuda.device);
    fprintf(stderr, "%s] CUDA v%d, %d devices, using #%d: '%s' with %zu MB of memory\n", G.workername, ver, ndev, devindex, buf, memsize/1024/1024);
    CU(cuDevicePrimaryCtxSetFlags, G.cuda.device, CU_CTX_SCHED_BLOCKING_SYNC);
    CU(cuDevicePrimaryCtxRetain, &G.cuda.context, G.cuda.device);
    CU(cuCtxSetCurrent, G.cuda.context);
    for (int thr = 0; thr < NTHREADS; thr++) {
        PerThreadState &ts = G.threads[thr];
        CU(cuStreamCreate, &ts.stream, CU_STREAM_NON_BLOCKING);
        CU(cuMemAlloc, &ts.params_d, sizeof(CudaParameterData));
        CU(cuMemAllocHost, (void**)&ts.params, sizeof(*ts.params));
        CU(cuMemAllocHost, (void**)&ts.result, sizeof(*ts.result));
        ts.buffer_size = 128*1024*1024;
        CU(cuMemAlloc, &ts.buffer_d, ts.buffer_size);
        CU(cuMemsetD8Async, ts.params_d, 0, sizeof(CudaParameterData), ts.stream);
        CU(cuMemsetD8Async, ts.buffer_d, 0, ts.buffer_size, ts.stream);
        pthread_create(&ts.thread, NULL, &worker_thread, (void*)&G.threads[thr]);
    }
}

// Parsing

static void
parse_fail()
{
    fprintf(stderr, "%s] input parsing failed:\n", G.workername);
    fprintf(stderr, "%s", G.input_line);
    for (char *p = G.input_line + 1; p < G.input_p; p++)
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
    long long x = strtoll(G.input_p, &end, 10);
    if (unlikely(G.input_p == end)) parse_fail();
    G.input_p = end;
    return (uint64_t)x;
}

static int64_t
parse_int()
{
    char *end = NULL;
    long long x = strtoll(G.input_p, &end, 10);
    if (unlikely(G.input_p == end)) parse_fail();
    G.input_p = end;
    return (int64_t)x;
}

static real_t
parse_real()
{
    char *end = NULL;
    real_t x = strtod(G.input_p, &end);
    if (unlikely(G.input_p == end)) parse_fail();
    G.input_p = end;
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

static void
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
            return;
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
    G.useful_time += t2 - t1;
}

// Main RPC cycle

static void
handle_one_command()
{
    match_c('[');
    uint64_t token = parse_uint();
    match_c(','); match_c('"');
    int c = input_getchar();
    if (c == 'i') {
        IntegrateCmd c = {token};
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
        return;
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
}

static void
fill_workername()
{
    char host[MAXNAME] = {};
    gethostname(host, sizeof(host));
    long pid = getpid();
    snprintf(G.workername, sizeof(G.workername), "%s:%ld:cuda", host, pid);
}

void
usage(const char *argv0)
{
    fprintf(stderr, "%s] usage: %s [-d cuda-device-index]\n", G.workername, argv0);
    exit(1);
}

int
main(int argc, char *argv[])
{
    fill_workername();
    int devindex = 0;
    for (int opt; (opt = getopt(argc, argv, "d:")) != -1;) {
        switch (opt) {
        case 'd': devindex = atoi(optarg); break;
        default: usage(argv[0]); break;
        }
    }
    if (optind < argc) usage(argv[0]);
    load_minicuda();
    init(devindex);
    setvbuf(stdin, NULL, _IOFBF, 1024*1024);
    setvbuf(stdout, NULL, _IOLBF, 1024*1024);
    setvbuf(stderr, NULL, _IOLBF, 1024*1024);
    double readt = 0;
    double lastt = 0;
    double t1 = timestamp();
    bool quit = false;
    while (!quit) {
        lastt = timestamp();
        if (getline(&G.input_line, &G.input_linesize, stdin) < 0) break;
        readt += timestamp() - lastt;
        G.input_p = G.input_line;
        handle_one_command();
    }
    double t2 = timestamp();
    for (int i = 0; i < NTHREADS; i++)
        G.useful_time += G.threads[i].useful_time;
    if (0) {
        fprintf(stderr, "%s] Done in %.3gs: %.3g%% useful time, %.3g%% read time; work ended %.3gs ago\n",
                G.workername, lastt-t1, 100*G.useful_time/(lastt-t1), 100*readt/(lastt-t1), t2-lastt);
    }
}

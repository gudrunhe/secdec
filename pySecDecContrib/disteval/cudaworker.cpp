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

#include "minicuda.h"

//#define unlikely(x) (x)
#define unlikely(x) __builtin_expect((x), 0)

typedef double real_t;
typedef struct { double re, im; } complex_t;

#define MAXPATH 4095
#define MAXNAME 255
#define MAXDIM 32

typedef int (*IntegrateF)(
    void * presult,
    const unsigned long lattice,
    const unsigned long index1,
    const unsigned long index2,
    const unsigned long * genvec,
    const real_t * shift,
    const real_t * realp,
    const complex_t * complexp,
    const real_t * deformp
);

typedef void (*MaxdeformpF)(
    real_t * deformp,
    const unsigned long lattice,
    const unsigned long index1,
    const unsigned long index2,
    const unsigned long * genvec,
    const real_t * shift,
    const real_t * realp,
    const complex_t * complexp
);

typedef int (*FpolycheckF)(
    const unsigned long lattice,
    const unsigned long index1,
    const unsigned long index2,
    const unsigned long * genvec,
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
    uint64_t kernelidx;
    uint64_t lattice;
    uint64_t i1;
    uint64_t i2;
    uint64_t genvec[MAXDIM];
    real_t shift[MAXDIM];
    real_t deformp[MAXDIM];
};

struct CudaParameterData {
    unsigned long genvec[MAXDIM];
    real_t shift[MAXDIM];
    real_t realp[MAXDIM];
    complex_t complexp[MAXDIM];
    real_t deformp[MAXDIM];
};

// Global state
static struct GlobalState {
    char workername[MAXNAME];
    std::vector<Family> families;
    std::vector<Kernel> kernels;
    char *input_line = NULL;
    char *input_p = NULL;
    size_t input_linesize = 0;
    struct GlobalCudaState {
        CUdevice device;
        CUcontext context;
        CUstream stream;
        CUdeviceptr buffer_d;
        size_t buffer_size;
        CUdeviceptr params_d;
        CudaParameterData *params;
        complex_t *result;
        CUmodule builtin_module;
        CUfunction fn_sum_d_b128_x1024;
        CUfunction fn_sum_c_b128_x1024;
    } cuda;
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
    fflush(stderr);
    exit(1);
}

static void
cuda_init()
{
    CU(cuInit, 0);
    int ver = 0, ndev = 0;
    size_t memsize = 0;
    CU(cuDriverGetVersion, &ver);
    CU(cuDeviceGetCount, &ndev);
    CU(cuDeviceGet, &G.cuda.device, 0);
    char buf[256];
    CU(cuDeviceGetName, buf, sizeof(buf), G.cuda.device);
    CU(cuDeviceGetName, buf, sizeof(buf), G.cuda.device);
    CU(cuDeviceTotalMem, &memsize, G.cuda.device);
    fprintf(stderr, "%s] CUDA v%d, %d devices, using '%s' with %zuMB of memory\n", G.workername, ver, ndev, buf, memsize/1024/1024);
    CU(cuDevicePrimaryCtxSetFlags, G.cuda.device, CU_CTX_SCHED_BLOCKING_SYNC);
    CU(cuDevicePrimaryCtxRetain, &G.cuda.context, G.cuda.device);
    CU(cuCtxPushCurrent, G.cuda.context);
    CU(cuStreamCreate, &G.cuda.stream, CU_STREAM_NON_BLOCKING);
    CU(cuMemAlloc, &G.cuda.params_d, sizeof(CudaParameterData));
    CU(cuMemAllocHost, (void**)&G.cuda.params, sizeof(*G.cuda.params));
    CU(cuMemAllocHost, (void**)&G.cuda.result, sizeof(*G.cuda.result));
    G.cuda.buffer_size = 128*1024*1024;
    CU(cuMemAlloc, &G.cuda.buffer_d, G.cuda.buffer_size);
    CU(cuMemsetD8Async, G.cuda.params_d, 0, sizeof(CudaParameterData), G.cuda.stream);
    CU(cuMemsetD8Async, G.cuda.buffer_d, 0, G.cuda.buffer_size, G.cuda.stream);
}

// Commands

static double
cmd_start(uint64_t token, StartCmd &c)
{
    int r = chdir(c.dirname);
    if (r != 0) {
        printf("@[%zu,null,\"failed to chdir %s: %d\"]\n", token, c.dirname, r);
        return 0;
    }
    CU(cuModuleLoad, &G.cuda.builtin_module, "./builtin.fatbin");
    CU(cuModuleGetFunction, &G.cuda.fn_sum_d_b128_x1024, G.cuda.builtin_module, "sum_d_b128_x1024");
    CU(cuModuleGetFunction, &G.cuda.fn_sum_c_b128_x1024, G.cuda.builtin_module, "sum_c_b128_x1024");
    printf("@[%zu,\"%s\",null]\n", token, G.workername);
    return 0;
}

static double
cmd_family(uint64_t token, FamilyCmd &c)
{
    assert(c.index == G.families.size());
    Family fam = {};
    char buf[MAXNAME+16];
    snprintf(buf, sizeof(buf), "./%s.so", c.name);
    fam.so_handle = dlopen(buf, RTLD_LAZY | RTLD_LOCAL);
    if (fam.so_handle == NULL) {
        printf("@[%zu,null,\"failed to open %s: %s\"]\n", token, buf, strerror(errno));
        return 0;
    }
    snprintf(buf, sizeof(buf), "./%s.fatbin", c.name);
    if (cuModuleLoad(&fam.cuda_module, buf) != 0) {
        printf("@[%zu,null,\"failed to open %s\"]\n", token, buf);
        return 0;
    }
    fam.dimension = c.dimension;
    memcpy(fam.realp, c.realp, sizeof(fam.realp));
    memcpy(fam.complexp, c.complexp, sizeof(fam.complexp));
    fam.complex_result = c.complex_result;
    memcpy(fam.name, c.name, sizeof(fam.name));
    G.families.push_back(fam);
    printf("@[%zu,null,null]\n", token);
    return 0;
}

static double
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
        printf("@[%zu,null,\"function not found: %s\"]\n", token, buf);
        return 0;
    }
    if (cuModuleGetFunction(&ker.cuda_fn_integrate, fam.cuda_module, buf) != 0) {
        printf("@[%zu,null,\"CUDA function not found: %s\"]\n", token, buf);
        return 0;
    }
    snprintf(buf, sizeof(buf), "%s__%s__maxdeformp", fam.name, c.name);
    ker.fn_maxdeformp = (MaxdeformpF)dlsym(fam.so_handle, buf);
    snprintf(buf, sizeof(buf), "%s__%s__fpolycheck", fam.name, c.name);
    ker.fn_fpolycheck = (FpolycheckF)dlsym(fam.so_handle, buf);
    memcpy(ker.name, c.name, sizeof(ker.name));
    G.kernels.push_back(ker);
    printf("@[%zu,null,null]\n", token);
    return 0;
}

static double
cmd_presample(uint64_t token, PresampleCmd &c)
{
    if (unlikely(c.kernelidx >= G.kernels.size())) {
        printf("@[%zu,null,\"kernel %zu was not loaded\"]\n", token, c.kernelidx);
        return 0;
    }
    const Kernel &ker = G.kernels[c.kernelidx];
    const Family &fam = G.families[ker.familyidx];
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

static void stupid_cuda_dummy(void *userdata)
{ (void)userdata; }

static double
cmd_integrate(uint64_t token, IntegrateCmd &c)
{
    if (unlikely(c.kernelidx >= G.kernels.size())) {
        printf("@[%zu,null,\"kernel %zu was not loaded\"]\n", token, c.kernelidx);
        return 0;
    }
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
            printf("@[%zu,[[NaN,NaN],%zu,%.4e],\"NaN != sign check error %d in %s.%s\"]", token, c.i2-c.i1, t2-t1, r, fam.name, ker.name);
        }
        printf("@[%zu,[", token); print_complex(result); printf(",%zu,%.4e],null]\n", c.i2-c.i1, t2-t1);
        return t2-t1;
    }
    if (1) { // CUDA path
        unsigned long threads = 128, pt_per_thread = 8;
        unsigned long blocks = (c.i2 - c.i1 + threads*pt_per_thread - 1)/(threads*pt_per_thread);
        unsigned long bufsize = fam.complex_result ? blocks*sizeof(complex_t) : blocks*sizeof(real_t);
        if (bufsize > G.cuda.buffer_size) {
            fprintf(stderr, "%s] realloc CUDA buffer to %zuMB\n", G.workername, bufsize/1024/1024);
            CU(cuMemFree, G.cuda.buffer_d);
            G.cuda.buffer_size = bufsize;
            CU(cuMemAlloc, &G.cuda.buffer_d, G.cuda.buffer_size);
            CU(cuMemsetD8Async, G.cuda.buffer_d, 0, G.cuda.buffer_size, G.cuda.stream);
        }
        memcpy(G.cuda.params->genvec, c.genvec, sizeof(c.genvec)); 
        memcpy(G.cuda.params->shift, c.shift, sizeof(c.shift)); 
        memcpy(G.cuda.params->realp, fam.realp, sizeof(fam.realp)); 
        memcpy(G.cuda.params->complexp, fam.complexp, sizeof(fam.complexp)); 
        memcpy(G.cuda.params->deformp, c.deformp, sizeof(c.deformp)); 
        double t1 = timestamp();
        CU(cuMemcpyHtoDAsync, G.cuda.params_d, G.cuda.params, sizeof(CudaParameterData), G.cuda.stream);
        CUdeviceptr genvec_d = G.cuda.params_d + offsetof(CudaParameterData, genvec);
        CUdeviceptr shift_d = G.cuda.params_d + offsetof(CudaParameterData, shift);
        CUdeviceptr realp_d = G.cuda.params_d + offsetof(CudaParameterData, realp);
        CUdeviceptr complexp_d = G.cuda.params_d + offsetof(CudaParameterData, complexp);
        CUdeviceptr deformp_d = G.cuda.params_d + offsetof(CudaParameterData, deformp);
        void *args[] = {&G.cuda.buffer_d, &c.lattice, &c.i1, &c.i2, &genvec_d, &shift_d, &realp_d, &complexp_d, &deformp_d, NULL };
        CU(cuLaunchKernel, ker.cuda_fn_integrate, blocks, 1, 1, threads, 1, 1, 0, G.cuda.stream, args, NULL);
        void *sum_args[] = {&G.cuda.buffer_d, &G.cuda.buffer_d, &blocks, NULL};
        CUfunction fn_sum = fam.complex_result ? G.cuda.fn_sum_c_b128_x1024 : G.cuda.fn_sum_d_b128_x1024;
        while (blocks > 1) {
            unsigned reduced = (blocks + 1024-1)/1024;
            CU(cuLaunchKernel, fn_sum, reduced, 1, 1, 128, 1, 1, 0, G.cuda.stream, sum_args, NULL);
            blocks = reduced;
        }
        CU(cuLaunchKernel, fn_sum, 1, 1, 1, 128, 1, 1, 0, G.cuda.stream, sum_args, NULL);
        G.cuda.result->re = 0;
        G.cuda.result->im = 0;
        CU(cuMemcpyDtoHAsync, G.cuda.result, G.cuda.buffer_d, fam.complex_result ? sizeof(complex_t) : sizeof(real_t), G.cuda.stream);
        // Without this CU_CTX_SCHED_BLOCKING_SYNC doesn't work,
        // and cuStreamSynchronize spins with 100% CPU usage.
        // With this, both CU_CTX_SCHED_BLOCKING_SYNC and
        // CU_CTX_SCHED_YIELD have the same result: 0% CPU usage
        // during cuStreamSynchronize.
        // How is cuLaunchHostFunc related though?
        // And how could one possibly find out about this?
        CU(cuLaunchHostFunc, G.cuda.stream, stupid_cuda_dummy, NULL);
        CU(cuStreamSynchronize, G.cuda.stream);
        double t2 = timestamp();
        printf("@[%zu,[", token); print_complex(*G.cuda.result); printf(",%zu,%.4e],null]\n", c.i2-c.i1, t2-t1);
        return t2-t1;
    }
}

static void
parse_fail()
{
    fprintf(stderr, "%s] input parsing failed:\n", G.workername);
    fprintf(stderr, "%s", G.input_line);
    for (char *p = G.input_line + 1; p < G.input_p; p++)
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
    long x = strtol(G.input_p, &end, 10);
    if (unlikely(G.input_p == end)) parse_fail();
    G.input_p = end;
    return (uint64_t)x;
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
    snprintf(G.workername, sizeof(G.workername), "%s:%ld:cuda", host, pid);
}

int main() {
    fill_workername();
    load_minicuda();
    cuda_init();
    setvbuf(stdout, NULL, _IOFBF, 1024*1024);
    double readt = 0;
    double workt = 0;
    double lastt = 0;
    double t1 = timestamp();
    bool quit = false;
    while (!quit) {
        lastt = timestamp();
        if (getline(&G.input_line, &G.input_linesize, stdin) < 0) break;
        readt += timestamp() - lastt;
        G.input_p = G.input_line;
        workt += handle_one_command();
        fflush(stdout);
    }
    double t2 = timestamp();
    fprintf(stderr, "%s] Done in %.3gs: %.3g%% useful time, %.3g%% read time; work ended %.3gs ago\n",
            G.workername, lastt-t1, 100*workt/(lastt-t1), 100*readt/(lastt-t1), t2-lastt);
    fflush(stderr);
}

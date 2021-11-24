/* This file provides a small subset of the CUDA driver API [1],
 * completely dynamically loaded, so that it could be compiled
 * on a machine without NVCC or any CUDA installation at all.
 *
 * [1] https://docs.nvidia.com/cuda/cuda-driver-api/
 */

#define CUDA_SONAME "libcuda.so.1"

#include <dlfcn.h>
#include <stdint.h>
#include <stdlib.h>

typedef enum CUresult_enum {} CUresult;
typedef int CUdevice;
typedef struct CUcontext_dummy *CUcontext;
typedef struct CUfunction_dummy *CUfunction;
typedef struct CUmodule_dummy *CUmodule;
typedef struct CUstream_dummy *CUstream;
typedef uintptr_t CUdeviceptr;
typedef void (*CUhostFn)(void *userData);

// cuStreamCreate
#define CU_STREAM_DEFAULT       0
#define CU_STREAM_NON_BLOCKING  1

// cuCtxCreate, cuDevicePrimaryCtxSetFlags
#define CU_CTX_SCHED_AUTO           0
#define CU_CTX_SCHED_SPIN           1
#define CU_CTX_SCHED_YIELD          2
#define CU_CTX_SCHED_BLOCKING_SYNC  4
#define CU_CTX_MAP_HOST             8

// cuDeviceGetAttribute
typedef enum CUdevice_attribute_enum {
    CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR = 75,
    CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR = 76
} CUdevice_attribute;

#define CUdeclare(fn, ...) \
    typedef CUresult (*fn ## _t)(__VA_ARGS__); \
    static fn ## _t fn = (fn ## _t)NULL;

CUdeclare(cuCtxCreate, CUcontext *ctx, unsigned int flags, CUdevice dev);
CUdeclare(cuCtxPushCurrent, CUcontext ctx);
CUdeclare(cuCtxSetCurrent, CUcontext ctx);
CUdeclare(cuDeviceGet, CUdevice *device, int ordinal);
CUdeclare(cuDeviceGetAttribute, int *val, CUdevice_attribute attrib, CUdevice dev);
CUdeclare(cuDeviceGetCount, int *count);
CUdeclare(cuDeviceGetName, char *name, int len, CUdevice dev);
CUdeclare(cuDevicePrimaryCtxRetain, CUcontext *ctx, CUdevice dev);
CUdeclare(cuDevicePrimaryCtxSetFlags, CUdevice dev, unsigned int flags);
CUdeclare(cuDeviceTotalMem, size_t *nbytes, CUdevice dev);
CUdeclare(cuDriverGetVersion, int *driverVersion);
CUdeclare(cuGetErrorName, CUresult error, const char **pstr);
CUdeclare(cuGetErrorString, CUresult error, const char **pstr);
CUdeclare(cuInit, unsigned int flags);
CUdeclare(cuLaunchHostFunc, CUstream stream, CUhostFn fn, void *userdata);
CUdeclare(cuLaunchKernel, CUfunction fn, unsigned int gridDimX, unsigned int gridDimY, unsigned int gridDimZ, unsigned int blockDimX, unsigned int blockDimY, unsigned int blockDimZ, unsigned int sharedMemBytes, CUstream stream, void **kernelParams, void **extra);
CUdeclare(cuMemAlloc, CUdeviceptr *dptr, size_t bytesize);
CUdeclare(cuMemAllocHost, void **pp, size_t bytesize);
CUdeclare(cuMemFree, CUdeviceptr dptr);
CUdeclare(cuMemcpyDtoHAsync, void *dstHost, CUdeviceptr srcDevice, size_t ByteCount, CUstream stream);
CUdeclare(cuMemcpyHtoDAsync, CUdeviceptr dstDevice, const void *srcHost, size_t ByteCount, CUstream stream);
CUdeclare(cuMemsetD8Async, CUdeviceptr dstDevice, unsigned char uc, size_t n, CUstream stream);
CUdeclare(cuModuleGetFunction, CUfunction *hfunc, CUmodule hmod, const char *name);
CUdeclare(cuModuleLoad, CUmodule *module, const char *fname);
CUdeclare(cuStreamCreate, CUstream *stream, unsigned int flags);
CUdeclare(cuStreamSynchronize, CUstream stream);

#define CUinit(name, symbol) \
    name = (name ## _t) dlsym(h, symbol); \
    if (name == (name ## _t)NULL) { \
        fprintf(stderr, "Failed to find %s in %s\n", symbol, CUDA_SONAME); \
        exit(1); \
    }

static void
load_minicuda()
{
    void *h = dlopen(CUDA_SONAME, RTLD_LAZY | RTLD_LOCAL);
    if (h == NULL) {
        fprintf(stderr, "Failed to load %s; is CUDA installed?\n", CUDA_SONAME);
        exit(1);
    }
    CUinit(cuCtxCreate,                 "cuCtxCreate_v2");
    CUinit(cuCtxPushCurrent,            "cuCtxPushCurrent_v2");
    CUinit(cuCtxSetCurrent,             "cuCtxSetCurrent");
    CUinit(cuDeviceGet,                 "cuDeviceGet");
    CUinit(cuDeviceGetAttribute,        "cuDeviceGetAttribute");
    CUinit(cuDeviceGetCount,            "cuDeviceGetCount");
    CUinit(cuDeviceGetName,             "cuDeviceGetName");
    CUinit(cuDevicePrimaryCtxRetain,    "cuDevicePrimaryCtxRetain");
    CUinit(cuDevicePrimaryCtxSetFlags,  "cuDevicePrimaryCtxSetFlags_v2");
    CUinit(cuDeviceTotalMem,            "cuDeviceTotalMem_v2");
    CUinit(cuDriverGetVersion,          "cuDriverGetVersion");
    CUinit(cuGetErrorName,              "cuGetErrorName");
    CUinit(cuGetErrorString,            "cuGetErrorString");
    CUinit(cuInit,                      "cuInit");
    CUinit(cuLaunchHostFunc,            "cuLaunchHostFunc");
    CUinit(cuLaunchKernel,              "cuLaunchKernel");
    CUinit(cuMemAlloc,                  "cuMemAlloc_v2");
    CUinit(cuMemAllocHost,              "cuMemAllocHost_v2");
    CUinit(cuMemFree,                   "cuMemFree_v2");
    CUinit(cuMemcpyDtoHAsync,           "cuMemcpyDtoHAsync_v2");
    CUinit(cuMemcpyHtoDAsync,           "cuMemcpyHtoDAsync_v2");
    CUinit(cuMemsetD8Async,             "cuMemsetD8Async");
    CUinit(cuModuleGetFunction,         "cuModuleGetFunction");
    CUinit(cuModuleLoad,                "cuModuleLoad");
    CUinit(cuStreamCreate,              "cuStreamCreate");
    CUinit(cuStreamSynchronize,         "cuStreamSynchronize");
}

#undef CUdeclare
#undef CUinit

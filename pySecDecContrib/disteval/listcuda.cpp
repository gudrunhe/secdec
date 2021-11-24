#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "minicuda.h"

#define CU(fn, ...) { \
        CUresult curesult = fn(__VA_ARGS__); \
        if (curesult != 0) cuda_fail(#fn, curesult); \
    }

static void
cuda_fail(const char *what, CUresult code)
{
    const char *n = NULL, *s = NULL;
    cuGetErrorName(code, &n);
    cuGetErrorString(code, &s);
    fprintf(stderr, "%s failed with code %d (%s): %s\n", what, code, n, s);
    fflush(stderr);
    exit(1);
}

int main() {
    load_minicuda();
    CU(cuInit, 0);
    int cudaver = 0;
    int ndev = 0;
    CU(cuDriverGetVersion, &cudaver);
    CU(cuDeviceGetCount, &ndev);
    for (int i = 0; i < ndev; i++) {
        CUdevice device;
        CU(cuDeviceGet, &device, i);
        char buf[256];
        CU(cuDeviceGetName, buf, sizeof(buf), device);
        size_t memsize = 0;
        CU(cuDeviceTotalMem, &memsize, device);
        int v_major = 0;
        int v_minor = 0;
        CU(cuDeviceGetAttribute, &v_major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, device);
        CU(cuDeviceGetAttribute, &v_minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, device);
        printf("{name: '%s', compute_capability: %d, memory: %zu}\n", buf, v_major*10 + v_minor, memsize/1024/1024);
    }
}

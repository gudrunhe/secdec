#!/bin/sh
tmp=${TMP:-/tmp}
cat >"$tmp/cudaarch.cu" << EOF
#include <stdio.h>
int main() {
  int ndevices = 0;
  if (cudaGetDeviceCount(&ndevices) != 0) {
    fprintf(stderr, "failed to list CUDA devices: %s\n",
        cudaGetErrorString(cudaGetLastError()));
  };
  for (int i = 0; i < ndevices; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    int v = prop.major * 10 + prop.minor;
    printf("-gencode arch=compute_%d,code=sm_%d\n", v, v);
  }
}
EOF
nvcc "$tmp/cudaarch.cu" -o "$tmp/cudaarch" || exit 1
"$tmp/cudaarch"
rm "$tmp/cudaarch.cu" "$tmp/cudaarch"

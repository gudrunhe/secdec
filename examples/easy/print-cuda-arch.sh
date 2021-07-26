#!/bin/sh
tmp=${TMP:-/tmp}
cat >"$tmp/cudaarch.cu" << EOF
#include <stdio.h>
int main() {
  int ndevices;
  cudaDeviceProp prop;
  cudaGetDeviceCount(&ndevices);
  for (int i = 0; i < ndevices; i++) {
    cudaGetDeviceProperties(&prop, i);
    int v = prop.major * 10 + prop.minor;
    printf("-gencode arch=compute_%d,code=sm_%d\n", v, v);
  }
}
EOF
nvcc "$tmp/cudaarch.cu" -o "$tmp/cudaarch" || exit 1
"$tmp/cudaarch"
rm "$tmp/cudaarch.cu" "$tmp/cudaarch"

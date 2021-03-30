#!/bin/sh

srcdir=$(dirname "$0")
img=baseimage.sif
dstdir=${TMP:-/tmp}/pysecdec-github-runner

cd "$dstdir"
singularity shell -C --overlay overlay.img -W work -B runner:/runner --pwd /runner $img

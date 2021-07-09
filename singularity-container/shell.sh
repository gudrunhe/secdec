#!/bin/sh

srcdir=$(dirname "$0")
img=baseimage.sif
dstdir=${TMP:-/tmp}/pysecdec-github-runner

cd "$dstdir"
imgpath=/homedir/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
singularity shell -C -W work -H homedir:/homedir -B runner:/runner --env=PATH=$imgpath --pwd /runner $img

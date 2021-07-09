#!/bin/sh

srcdir=$(dirname "$0")
img=baseimage.sif
dstdir=${TMP:-/tmp}/pysecdec-github-runner

if [ ! -d "$dstdir" ]; then
    echo "error: $0: directory does not exist: $dstdir" 1>&2
    exit 1
fi

cd "$dstdir"
imgpath=/homedir/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
singularity exec -C -W work -H homedir:/homedir -B runner:/runner --env=PATH=$imgpath --pwd /runner $img nice /runner/run.sh

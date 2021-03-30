#!/bin/sh

srcdir=$(dirname "$0")
img=baseimage.sif
dstdir=${TMP:-/tmp}/pysecdec-github-runner

if [ ! -d "$dstdir" ]; then
    echo "error: $0: directory does not exist: $dstdir" 1>&2
    exit 1
fi

cd "$dstdir"
singularity exec -C --overlay overlay.img -W work -B runner:/runner $img nice /runner/run.sh

#!/bin/sh

srcdir=$(dirname "$0")
img=baseimage.sif
tmpdir=$(mktemp -d)
tmpimg=${tmpdir}/$img

singularity build "$tmpimg" "$srcdir/${img%.sif}.def" || exit 1
mv "$tmpimg" "$srcdir/$img"

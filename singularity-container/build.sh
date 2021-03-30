#!/bin/sh

srcdir=$(dirname "$0")
img=baseimage.sif
tmpimg=${TMP:-/tmp}/$img

sudo singularity build "$tmpimg" "$srcdir/${img%.sif}.def" || exit 1
mv "$tmpimg" "$srcdir/$img"

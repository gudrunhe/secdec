#!/bin/sh

srcdir=$(dirname "$0")

daemon -n pysecdec-github-runner -P /tmp -o ${TMP:-/tmp}/pysecdec-github-runner/daemon.log ${srcdir}/run.sh

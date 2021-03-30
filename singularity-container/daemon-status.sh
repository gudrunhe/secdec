#!/bin/sh

if daemon -n pysecdec-github-runner -P /tmp --running; then
    echo "Status: running"
else
    echo "Status: not running"
fi
echo "=== Log output ==="
tail -20 ${TMP:-/tmp}/pysecdec-github-runner/daemon.log

#!/bin/sh

daemon -n pysecdec-github-runner -P /tmp --stop
sleep 1
pkill Runner.Listener
sleep 1
daemon -n pysecdec-github-runner -P /tmp --stop

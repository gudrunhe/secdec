#!/usr/bin/env python3
"""
Distributed pySecDec integral evaluator.
Usage: python3 -m pySecDecWorker [--cpu | --cuda]
"""

import os
import os.path
import sys

import pySecDecContrib

if __name__ == "__main__":
    worker_type = "cpu"
    for arg in sys.argv[1:]:
        if arg == "--cpu": worker_type = "cpu"
        elif arg == "--cuda":
            raise ValueError(f"CUDA worker does not work yet")
        else:
            raise ValueError(f"Unknown option: {arg}")
    binpath = os.path.join(pySecDecContrib.dirname, "bin", f"pysecdec_{worker_type}worker")
    os.execl(binpath, binpath)

#!/usr/bin/env python3
"""
Distributed pySecDec integral evaluator.
Usage: python3 -m pySecDec.dworker [--cpu | --cuda] [dirname]
"""

import json
import numpy as np
import os
import time
import sys
import socket
import re
import pickle
import base64

def log(text):
    print(f"{WORKERNAME}]", text, file=sys.stderr)
    sys.stderr.flush()

class CudaWorker:

    def __init__(self):
        self.narguments = self.nrealp = self.ncomplexp = self.ndeformp = self.nblocks = 0
        self.genvec_d = self.shift_d = self.realp_d = self.complexp_d = self.deformp_d = 0
        self.reallocate(16, 2**24, 16, 16, 16)
        self.cubin = {}
        self.kernels = {}
        self.first_t = None
        self.last_t = None
        self.sum_t = 0.0

    def start(self, dirname):
        self.load_builtin_kernels(dirname)

    def load_builtin_kernels(self, dirname):
        self.sum_cubin = cuda.module_from_file(os.path.join(dirname, "builtin_cuda.fatbin"))
        self.sum_kernels = {}
        self.sum_kernels[np.float64] = self.sum_cubin.get_function("sum_d_b128_x1024")
        self.sum_kernels[np.float64].prepare("PPQ")
        self.sum_kernels[np.complex128] = self.sum_cubin.get_function("sum_c_b128_x1024")
        self.sum_kernels[np.complex128].prepare("PPQ")
        self.gauge_kernel = self.sum_cubin.get_function("builtin__gauge")
        self.gauge_kernel.prepare("PQQQPPPPP")

    def load(self, jsonfiles):
        for jsonfile in jsonfiles:
            with open(jsonfile, "r") as f:
                info = json.load(f)
            dirname = os.path.dirname(jsonfile)
            self.cubin[info["name"]] = cuda.module_from_file(os.path.join(dirname, info["name"] + ".fatbin"))

    def reallocate(self, narguments, nblocks, nrealp, ncomplexp, ndeformp):
        if narguments > self.narguments:
            self.narguments = narguments
            self.genvec_d = cuda.mem_alloc(self.narguments*8)
            self.shift_d = cuda.mem_alloc(self.narguments*8)
        if nblocks > self.nblocks:
            self.nblocks = nblocks
            log(f"result realloc to {self.nblocks}")
            self.result_d = cuda.mem_alloc(self.nblocks*16)
        if nrealp > self.nrealp:
            self.nrealp = nrealp
            self.realp_d = cuda.mem_alloc(self.nrealp*8)
        if ncomplexp > self.ncomplexp:
            self.ncomplexp = ncomplexp
            self.complexp_d = cuda.mem_alloc(self.ncomplexp*16)
        if ndeformp > self.ndeformp:
            self.ndeformp = ndeformp
            self.deformp_d = cuda.mem_alloc(self.ndeformp*8)

    def integrate_kernel(self, kernel, narguments, lattice, i1, i2, genvec, shift, realp=None, complexp=None, deformp=None, complex_result=False):
        threads = 128
        ppthread = 8
        blocks = (i2 - i1 + threads*ppthread - 1)//(threads*ppthread)
        self.reallocate(narguments, blocks,
                len(realp) if realp is not None else 0,
                len(complexp) if complexp is not None else 0,
                len(deformp) if deformp is not None else 0)
        cuda.memcpy_htod(self.genvec_d, genvec)
        cuda.memcpy_htod(self.shift_d, shift)
        if realp is not None: cuda.memcpy_htod(self.realp_d, realp)
        if complexp is not None: cuda.memcpy_htod(self.complexp_d, complexp)
        if deformp is not None: cuda.memcpy_htod(self.deformp_d, deformp)
        kernel.prepared_call(
            (blocks, 1), (threads, 1, 1),
            self.result_d, lattice, i1, i2, self.genvec_d, self.shift_d,
            self.realp_d, self.complexp_d, self.deformp_d)
        dtype = np.complex128 if complex_result else np.float64
        while blocks > 1:
            red = (blocks + 1024-1)//1024
            self.sum_kernels[dtype].prepared_call((red, 1), (128, 1, 1), self.result_d, self.result_d, blocks)
            blocks = red
        self.sum_kernels[dtype].prepared_call((1, 1), (128, 1, 1), self.result_d, self.result_d, blocks)
        res = np.zeros(1, dtype=dtype)
        cuda.memcpy_dtoh(res, self.result_d)
        return res[0]

    def integrate(self, iname, kname, narguments, lattice, i1, i2, genvec, shift, realp=None, complexp=None, deformp=None, complex_result=False):
        try:
            if (iname, kname) not in self.kernels:
                self.kernels[iname, kname] = self.cubin[iname].get_function(iname + "__" + kname)
                self.kernels[iname, kname].prepare("PQQQPPPPP")
            t1 = time.time()
            res = self.integrate_kernel(self.kernels[iname, kname], narguments, lattice, i1, i2, genvec, shift,
                    realp=realp, complexp=complexp, deformp=deformp, complex_result=complex_result)
            t2 = time.time()
            if self.first_t is None: self.first_t = t1
            self.last_t = t2
            self.sum_t += t2-t1
            return res, i2-i1, t2-t1
        except Exception as e:
            log(f"integrate({iname}.{kname}) failed: {e}")
            raise

    def integrate_gauge(self, lattice, i1, i2, genvec, shift, realp, complexp, deformp):
        t1 = time.time()
        res = self.integrate_kernel(self.gauge_kernel, 2, lattice, i1, i2, genvec, shift, realp, complexp, deformp, True)
        t2 = time.time()
        if self.first_t is None: self.first_t = t1
        self.last_t = t2
        self.sum_t += t2-t1
        return res, i2-i1, t2-t1

    def ping(self, *args, **kwargs):
        return (args, kwargs)

class CPUWorker:

    def __init__(self):
        self.dlls = {}
        self.kernels = {}
        self.first_t = None
        self.last_t = None
        self.sum_t = 0.0

    def start(self, dirname):
        self.load_dll("builtin", os.path.join(dirname, "builtin_cpu.so"))

    def load(self, jsonfiles):
        for jsonfile in jsonfiles:
            with open(jsonfile, "r") as f:
                info = json.load(f)
            dirname = os.path.dirname(jsonfile)
            self.load_dll(info["name"], os.path.join(dirname, info["name"] + ".so"))

    def load_dll(self, key, sofile):
        self.dlls[key] = ctypes.cdll.LoadLibrary(os.path.abspath(sofile))

    def integrate(self, iname, kname, narguments, lattice, i1, i2, genvec, shift, realp, complexp, deformp, complex_result):
        if (iname, kname) not in self.kernels:
            fun = getattr(self.dlls[iname], iname + "__" + kname)
            fun.argtypes = [
                ctypes.c_void_p, # result
                ctypes.c_ulong, # lattice
                ctypes.c_ulong, # index1
                ctypes.c_ulong, # index2
                ctypes.c_void_p, # genvec
                ctypes.c_void_p, # shift
                ctypes.c_void_p, # realp
                ctypes.c_void_p, # complexp
                ctypes.c_void_p # deformp
            ]
            fun.restype = ctypes.c_int
            self.kernels[iname, kname] = fun
        t1 = time.time()
        result = np.zeros(1, dtype=np.complex128 if complex_result else np.float64)
        self.kernels[iname, kname](
            result.ctypes.data_as(ctypes.c_void_p),
            lattice, i1, i2,
            genvec.ctypes.data_as(ctypes.c_void_p),
            shift.ctypes.data_as(ctypes.c_void_p),
            realp.ctypes.data_as(ctypes.c_void_p) if realp is not None else 0,
            complexp.ctypes.data_as(ctypes.c_void_p) if complexp is not None else 0,
            deformp.ctypes.data_as(ctypes.c_void_p) if deformp is not None else 0)
        t2 = time.time()
        if self.first_t is None: self.first_t = t1
        self.last_t = t2
        self.sum_t += t2-t1
        return result[0], i2-i1, t2-t1

    def integrate_gauge(self, lattice, i1, i2, genvec, shift, realp, complexp, deformp):
        return self.integrate("builtin", "gauge", 2, lattice, i1, i2, genvec, shift, realp, complexp, deformp, True)

    def ping(self, *args, **kwargs):
        return (args, kwargs)

def encode_bin(obj):
    return base64.b64encode(pickle.dumps(obj))

def decode_bin(data):
    return pickle.loads(base64.b64decode(data))

def respond(message):
    sys.stdout.buffer.write(message)
    sys.stdout.buffer.flush()

if __name__ == "__main__":

    WORKERNAME = f"{socket.gethostname()}.{os.getpid()}"

    worker_type = "cpu"
    for arg in sys.argv[1:]:
        if arg == "--cpu": worker_type = "cpu"
        elif arg == "--cuda": worker_type = "cuda"
        else:
            raise ValueError(f"Unknown option: {arg}")

    if worker_type == "cuda":
        log("Starting CUDA worker")
        import pycuda.autoinit
        import pycuda.driver as cuda
        ii = CudaWorker()

    if worker_type == "cpu":
        log("Starting CPU worker")
        import ctypes
        ii = CPUWorker()

    read_t = 0.0
    while True:
        t0 = time.time()
        line = sys.stdin.readline()
        if line == "": break
        t1 = time.time()
        if ii.first_t is not None:
            read_t += t1 - t0
        m = re.match("@([a-zA-Z0-9_]+)(?: ([^\\n]*))?\n", line)
        if not m:
            print(b"??? " + repr(line).encode("utf-8"))
            sys.stdout.flush()
            continue
        cmd = m.group(1)
        arg = m.group(2)
        if cmd == "call":
            try:
                arg = decode_bin(arg)
                result = getattr(ii, arg["m"])(*arg["a"])
                respond(b"@return " + encode_bin({"i": arg["i"], "r": result}) + b"\n")
            except Exception as e:
                respond(b"@return " + encode_bin({"i": arg["i"], "e":
                    Exception("Worker " + str(type(e).__name__) + ": " + str(e))}) + b"\n")
        elif cmd == "call_json":
            arg = json.loads(arg)
            kwargs = {k: np.array(v) if isinstance(v, list) else v for k, v in arg["args"].items()}
            result = getattr(ii, arg["method"])(**kwargs)
            print("@return", str(result))
        elif cmd == "start":
            ii.start(decode_bin(arg) if arg else ".")
            respond(f"@started {WORKERNAME}\n".encode("utf-8"))
        elif cmd == "stop":
            break
        else:
            print("??? " + repr(line))
        sys.stdout.flush()

    dt = time.time() - ii.last_t
    allt = ii.last_t - ii.first_t
    frac_int = ii.sum_t/allt
    frac_rd = read_t/allt
    log(f"Done in {allt:.3f}s; {100*frac_int:.1f}% int time, {100*frac_rd:.1f}% read time; work ended {dt:.3f}s ago")

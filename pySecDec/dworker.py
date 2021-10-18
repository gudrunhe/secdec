#!/usr/bin/env python3
"""
Distributed pySecDec integral evaluator.
Usage: python3 -m pySecDec.dworker [--cpu | --cuda]
"""

import base64
import json
import numpy as np
import os
import pickle
import re
import socket
import sys
import time

def log(text):
    print(f"{WORKERNAME}]", text, file=sys.stderr)
    sys.stderr.flush()

class CudaWorker:

    def __init__(self):
        self.narguments = self.nrealp = self.ncomplexp = self.ndeformp = self.nblocks = 0
        self.genvec_d = self.shift_d = self.realp_d = self.complexp_d = self.deformp_d = 0
        self._reallocate(16, 2**24, 16, 16, 16)
        self.cubin = {}
        self.kernels = {}

    def start(self, dirname):
        self.sum_cubin = cuda.module_from_file(os.path.join(dirname, "builtin_cuda.fatbin"))
        self.cubin["builtin"] = self.sum_cubin
        self.kernels["builtin", "gauge"] = self.sum_cubin.get_function("builtin__gauge")
        self.kernels["builtin", "gauge"].prepare("PQQQPPPPP")
        self.sum_kernels = {}
        self.sum_kernels[np.float64] = self.sum_cubin.get_function("sum_d_b128_x1024")
        self.sum_kernels[np.float64].prepare("PPQ")
        self.sum_kernels[np.complex128] = self.sum_cubin.get_function("sum_c_b128_x1024")
        self.sum_kernels[np.complex128].prepare("PPQ")
        return WORKERNAME, 0

    def load(self, jsonfiles):
        for jsonfile in jsonfiles:
            with open(jsonfile, "r") as f:
                info = json.load(f)
            dirname = os.path.dirname(jsonfile)
            self.cubin[info["name"]] = cuda.module_from_file(os.path.join(dirname, info["name"] + ".fatbin"))
        return None, 0

    def _reallocate(self, narguments, nblocks, nrealp, ncomplexp, ndeformp):
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

    def _integrate_kernel(self, kernel, narguments, lattice, i1, i2, genvec, shift, realp, complexp, deformp, complex_result):
        threads = 128
        ppthread = 8
        blocks = (i2 - i1 + threads*ppthread - 1)//(threads*ppthread)
        self._reallocate(narguments, blocks,
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

    def integrate(self, iname, kname, narguments, lattice, i1, i2, genvec, shift, realp, complexp, deformp, complex_result):
        try:
            kernel = self.kernels.get((iname, kname), None)
            if kernel is None:
                kernel = self.cubin[iname].get_function(iname + "__" + kname)
                kernel.prepare("PQQQPPPPP")
                self.kernels[iname, kname] = kernel
            t0 = time.time()
            res = self._integrate_kernel(kernel, narguments, lattice, i1, i2, genvec, shift,
                    realp=realp, complexp=complexp, deformp=deformp, complex_result=complex_result)
            dt = time.time() - t0
            return (res, i2-i1, dt), dt
        except Exception as e:
            log(f"integrate({iname}.{kname}) failed: {e}")
            raise

    def ping(self, *args):
        return args, 0

class CPUWorker:

    def __init__(self):
        self.dlls = {}
        self.kernels = {}

    def start(self, dirname):
        self._load_dll("builtin", os.path.join(dirname, "builtin_cpu.so"))
        return WORKERNAME, 0

    def load(self, jsonfiles):
        for jsonfile in jsonfiles:
            with open(jsonfile, "r") as f:
                info = json.load(f)
            dirname = os.path.dirname(jsonfile)
            self._load_dll(info["name"], os.path.join(dirname, info["name"] + ".so"))
        return None, 0

    def _load_dll(self, key, sofile):
        self.dlls[key] = ctypes.cdll.LoadLibrary(os.path.abspath(sofile))

    def integrate(self, iname, kname, narguments, lattice, i1, i2, genvec, shift, realp, complexp, deformp, complex_result):
        fun = self.kernels.get((iname, kname), None)
        if fun is None:
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
        result = np.zeros(1, dtype=np.complex128 if complex_result else np.float64)
        t0 = time.time()
        code = fun(
            result.ctypes.data_as(ctypes.c_void_p),
            lattice, i1, i2,
            genvec.ctypes.data_as(ctypes.c_void_p),
            shift.ctypes.data_as(ctypes.c_void_p),
            realp.ctypes.data_as(ctypes.c_void_p) if realp is not None else 0,
            complexp.ctypes.data_as(ctypes.c_void_p) if complexp is not None else 0,
            deformp.ctypes.data_as(ctypes.c_void_p) if deformp is not None else 0)
        dt = time.time() - t0
        if np.isnan(result) and code == 0:
            log(f"bad nan with {iname}.{kname} deformp={deformp} shift={shift}")
        return (result[0], i2-i1, dt), dt

    def maxdeformp(self, iname, kname, ndeformp, lattice, genvec, shift, realp, complexp):
        fun_maxdeformp = getattr(self.dlls[iname], iname + "__" + kname + "__maxdeformp")
        fun_maxdeformp.argtypes = [
            ctypes.c_void_p, # result
            ctypes.c_ulong, # lattice
            ctypes.c_ulong, # index1
            ctypes.c_ulong, # index2
            ctypes.c_void_p, # genvec
            ctypes.c_void_p, # shift
            ctypes.c_void_p, # realp
            ctypes.c_void_p # complexp
        ]
        fun_fpolycheck = getattr(self.dlls[iname], iname + "__" + kname + "__fpolycheck")
        fun_fpolycheck.argtypes = [
            ctypes.c_ulong, # lattice
            ctypes.c_ulong, # index1
            ctypes.c_ulong, # index2
            ctypes.c_void_p, # genvec
            ctypes.c_void_p, # shift
            ctypes.c_void_p, # realp
            ctypes.c_void_p, # complexp
            ctypes.c_void_p # deformp
        ]
        fun_fpolycheck.restype = ctypes.c_int
        t0 = time.time()
        deformp = np.zeros(ndeformp, dtype=np.float64)
        fun_maxdeformp(
            deformp.ctypes.data_as(ctypes.c_void_p),
            lattice, 0, lattice,
            genvec.ctypes.data_as(ctypes.c_void_p),
            shift.ctypes.data_as(ctypes.c_void_p),
            realp.ctypes.data_as(ctypes.c_void_p) if realp is not None else 0,
            complexp.ctypes.data_as(ctypes.c_void_p) if complexp is not None else 0)
        while True:
            if fun_fpolycheck(lattice, 0, lattice,
                genvec.ctypes.data_as(ctypes.c_void_p),
                shift.ctypes.data_as(ctypes.c_void_p),
                realp.ctypes.data_as(ctypes.c_void_p) if realp is not None else 0,
                complexp.ctypes.data_as(ctypes.c_void_p) if complexp is not None else 0,
                deformp.ctypes.data_as(ctypes.c_void_p)) == 0:
                break
            deformp *= 0.9
        dt = time.time() - t0
        return deformp, dt

    def ping(self, *args):
        return args, 0

def decode_message(message):
    return pickle.loads(base64.b64decode(message))

def respond(ofile, response):
    data = base64.b64encode(pickle.dumps(response))
    ofile.write(b"@" + data + b"\n")

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
        wrkr = CudaWorker()

    if worker_type == "cpu":
        log("Starting CPU worker")
        import ctypes
        wrkr = CPUWorker()

    first_t = None
    read_t = 0.0
    work_t = 0.0
    rx = re.compile("@([a-zA-Z0-9_]+)(?: ([^\\n]*))?\n")
    ifile = sys.stdin.buffer
    ofile = sys.stdout.buffer
    while True:
        t0 = time.time()
        line = ifile.readline()
        t1 = time.time()
        if len(line) == 0: break
        for tag, cmd, arg in decode_message(line):
            try:
                result, dt = getattr(wrkr, cmd)(*arg)
                respond(ofile, (tag, result, None))
                if dt is not 0:
                    if first_t is None: first_t = t0
                    work_t += dt
            except Exception as e:
                respond(ofile, (tag, None, f"{WORKERNAME}: {type(e).__name__}: {e}"))
        ofile.flush()
        if first_t is not None:
            read_t += t1 - t0

    if first_t is not None:
        dt = t1 - t0
        all_t = t0 - first_t
        log(f"Done in {all_t:.3f}s; {100*work_t/all_t:.1f}% work time, {100*read_t/all_t:.1f}% read time; work ended {dt:.3f}s ago")

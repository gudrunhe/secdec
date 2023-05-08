#!/usr/bin/env python3
"""
Benchmark the evaluation rate (i.e. evaluations per second) of
disteval kernels in a given integral sum. Print the resulting
evaluation rates in a comma-separated value (CSV) format.
Usage:
    python3 -m pySecDec.disteval_benchmark integrand.json [options] <var>=value ...
Options:
    --points=X          evaluate using this many points per batch (default: 1e5)
    --repetitions=X     repeat the measurement this many times (default: 10)
    --help              show this help message
Arguments:
    <var>=X             set this integral or coefficient variable to a given value
"""

import asyncio
import getopt
import json
import math
import numpy as np
import os
import random
import re
import subprocess
import sympy as sp
import sys
import time

from .generating_vectors import generating_vector, max_lattice_size

from pySecDecContrib import dirname as contrib_dirname

log_starttime = time.time()
def log(*args):
    t = time.time() - log_starttime
    if t < 60:
        print(f"{t:.3f}]", *args, file=sys.stderr)
    else:
        m, s = divmod(t, 60)
        if m < 60:
            print(f"{m:.0f}:{s:06.3f}]", *args, file=sys.stderr)
        else:
            h, m = divmod(m, 60)
            print(f"{h:.0f}:{m:02.0f}:{s:06.3f}]", *args, file=sys.stderr)
    sys.stderr.flush()

# Generic RPC

def encode_message(data):
    return json.dumps(data,separators=(',',':')).encode("ascii") + b"\n"

def decode_message(binary):
    return json.loads(binary[1:])

class WorkerException(Exception):
    pass

class Worker:

    def __init__(self, process, name=None):
        self.name = name
        self.process = process
        self.serial = 0
        self.callbacks = {}
        self.reader_task = asyncio.get_event_loop().create_task(self._reader())

    def queue_size(self):
        return len(self.callbacks)

    def call_cb(self, method, args, callback, callback_args=()):
        token = self.serial = self.serial + 1
        self.callbacks[token] = (callback, callback_args)
        message = encode_message((token, method, args))
        self.process.stdin.write(message)
        return token
    
    def cancel_cb(self, token):
        if token in self.callbacks:
            self.callbacks[token] = (lambda a,b,c:None, ())

    def call(self, method, *args):
        fut = asyncio.futures.Future()
        def call_return(result, error, w):
            if error is None: fut.set_result(result)
            else: fut.set_exception(WorkerException(error))
        self.call_cb(method, args, call_return)
        return fut
    
    def multicall(self, calls):
        fut = asyncio.futures.Future()
        results = [None]*len(calls)
        ntodo = [len(calls)]
        def multicall_return(result, error, w, i):
            if error is None:
                results[i] = result
                ntodo[0] -= 1
                if ntodo[0] == 0:
                    fut.set_result(results)
            else:
                fut.set_exception(WorkerException(error))
        s0 = self.serial
        self.serial += len(calls)
        parts = []
        for i, (method, args) in enumerate(calls):
            parts.append(encode_message((s0 + i, method, args)))
            self.callbacks[s0 + i] = (multicall_return, (i,))
        self.process.stdin.write(b"".join(parts))
        return fut

    async def _reader(self):
        rx = re.compile(b"^@([a-zA-Z0-9_]+)(?: ([^\\n]*))?\n$")
        line = None
        try:
            while True:
                line = await self.process.stdout.readline()
                if len(line) == 0: break
                if line.startswith(b"@"):
                    i, res, err = decode_message(line)
                    callback, callback_args = self.callbacks[i]
                    del self.callbacks[i]
                    callback(res, err, self, *callback_args)
                else:
                    log(f"{self.name}: {line}")
        except Exception as e:
            log(f"{self.name} reader failed: {type(e).__name__}: {e}")
            log(f"{self.name} line was {line!r}")
        log(f"{self.name} reader exited")

async def launch_worker(command, dirname, maxtimeout=10):
    timeout = min(1, maxtimeout/10)
    while True:
        log(f"running: {command}")
        if isinstance(command, str):
            p = await asyncio.create_subprocess_shell(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        else:
            p = await asyncio.create_subprocess_exec(*command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        p.stdin.write(encode_message((0, "start", (dirname,))))
        answer = await p.stdout.readline()
        try:
            _, name, err = decode_message(answer)
            if err is not None:
                log(f"worker startup fail: {err}")
            else:
                w = Worker(p, name=name)
                log(f"worker {w.name} connected")
                return w
        except Exception as e:
            log(f"failed to start worker: {type(e).__name__}: {e}")
        try:
            p.stdin.close()
            p.kill()
        except ProcessLookupError:
            pass
        log(f"will retry after {timeout}s")
        await asyncio.sleep(timeout)
        timeout = min(timeout*2, maxtimeout)


# Main

async def dobenchmark(workercmd, datadir, intfile, valuemap, npoints, nreps):
    # Load the integrals from the requested json file
    with open(intfile, "r") as f:
        info = json.load(f)

    if info["type"] not in ("integral", "sum"):
        raise ValueError(f"unknown type: {info['type']}")

    for p in info["realp"] + info["complexp"]:
        if p not in valuemap:
            raise ValueError(f"missing integral parameter: {p}")

    sp_regulators = sp.var(info["regulators"])
    kernel2idx = {}
    if info["type"] == "integral":
        infos = {info["name"] : info}
        for k in info["kernels"]:
            kernel2idx[info["name"], k] = len(kernel2idx)
        log(f"got the total of {len(kernel2idx)} kernels")
    elif info["type"] == "sum":
        log(f"loading {len(info['integrals'])} integrals")
        infos = {}
        for i in info["integrals"]:
            with open(os.path.join(datadir, f"{i}.json"), "r") as f:
                infos[i] = json.load(f)
                assert infos[i]["name"] == i
            for k in infos[i]["kernels"]:
                kernel2idx[i, k] = len(kernel2idx)
        log(f"got the total of {len(kernel2idx)} kernels")
        if isinstance(info["sums"], list):
            info["sums"] = {f"sum{sumidx}" : sum for sumidx, sum in enumerate(info["sums"])}
    else:
        raise ValueError(f"unknown type: {info['type']}")
    
    realp = {
        i : [valuemap[p] for p in info["realp"]]
        for i, info in infos.items()
    }
    complexp = {
        i : [(np.real(valuemap[p]), np.imag(valuemap[p])) for p in info["complexp"]]
        for i, info in infos.items()
    }
    family2idx = {fam:i for i, fam in enumerate(infos.keys())}

    # Launch the worker
    w = await launch_worker(workercmd, datadir)
    await w.call("family", 0, "builtin", 2, (2.0, 0.1, 0.2, 0.3), (), True)
    await w.call("kernel", 0, 0, "gauge")
    await w.multicall([
        ("family", (i+1, fam, info["dimension"], realp[fam], complexp[fam], info["complex_result"]))
        for i, (fam, info) in enumerate(infos.items())
    ])
    await w.multicall([
        ("kernel", (i+1, family2idx[fam]+1, ker))
        for (fam, ker), i in kernel2idx.items()
    ])
    log(f"worker: {w.name}")

    # Benchmark each kernel
    results = []
    for (fam, ker), keridx in kernel2idx.items():
        dim = infos[fam]["dimension"]
        log(f"{fam}.{ker} (dim={dim}):")
        lattice, genvec = generating_vector(dim, npoints)
        shift = [random.random() for i in range(dim)]
        deformp = [1e-10 for i in range(dim)]
        dn = np.zeros(nreps)
        dt = np.zeros(nreps)
        for i in range(nreps):
            (re, im), dn[i], dt[i] = await w.call("integrate", keridx+1, lattice, 0, lattice, genvec, shift, deformp)
            if math.isnan(re) or math.isnan(im):
                log(f"- got ({re},{im}) at {lattice:.3e} points")
        rate = float(np.mean(dn/dt))
        rateerr = float(np.std(dn/dt))/math.sqrt(nreps)
        tm = float(np.mean(dt))
        tmerr = float(np.std(dt))/math.sqrt(nreps)
        log(f"- {rate:.3e} ± {rateerr:.1e} evals/s at {lattice:.3e} points in {tm:.3e} ± {tmerr:.1e}s")
        results.append((fam, ker, rate, rateerr))

    # Print the final statistics
    print("family,kernel,rate,rate_error")
    for result in results:
        print(",".join(map(str, result)))

def main():

    valuemap = {}
    npoints = 10**5
    nreps = 10
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "", ["points=", "repetitions=", "help"])
    except getopt.GetoptError as e:
        print(e, file=sys.stderr)
        print("use --help to see the usage", file=sys.stderr)
        exit(1)
    for key, value in opts:
        if key == "--points":
            npoints = int(float(value))
        elif key == "--repetitions":
            nreps = int(float(value))
        elif key == "--help":
            print(__doc__.strip())
            exit(0)
    if len(args) < 1:
        print(__doc__.strip(), file=sys.stderr)
        exit(1)
    intfile = args[0]
    dirname = os.path.dirname(intfile)
    log("Settings:")
    log(f"- file = {intfile}")
    log(f"- points = {npoints}")
    log(f"- repetitions = {nreps}")
    for arg in args[1:]:
        if "=" not in arg: raise ValueError(f"Bad argument: {arg}")
        key, value = arg.split("=", 1)
        value = complex(value)
        value = value.real if value.imag == 0 else value
        valuemap[key] = value
    log("Invariants:")
    for key, value in valuemap.items():
        log(f"- {key} = {value}")

    # Define the worker
    ncuda = 0
    if os.path.exists(os.path.join(dirname, "builtin.fatbin")):
        try:
            p = subprocess.run(
                    [os.path.join(contrib_dirname, "bin", "pysecdec_listcuda")],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.DEVNULL)
            ncuda = len(p.stdout.splitlines())
        except Exception as e:
            log(f"Can't determine GPU count: {e}")
    else:
        log(f"CUDA worker data was not built, skipping")
    if ncuda != 0:
        worker = ["nice", sys.executable, "-m", "pySecDecContrib", "pysecdec_cudaworker", "-d", "0"]
    else:
        worker = ["nice", sys.executable, "-m", "pySecDecContrib", "pysecdec_cpuworker"]

    # Run the benchmark
    loop = asyncio.get_event_loop()
    result = loop.run_until_complete(dobenchmark(worker, dirname, intfile, valuemap, npoints, nreps))

if __name__ == "__main__":
    main()

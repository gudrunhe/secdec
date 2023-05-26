#!/usr/bin/env python3
"""
Distributed and/or local evaluator for pySecDec integrals.
Usage:
    python3 -m pySecDec.deval integrand.json [options] <var>=value ...
Options:
    --epsabs=X          stop if this absolute precision is reached (default: 1e-10)
    --epsrel=X          stop if this relative precision is reached (default: 1e-4)
    --timeout=X         stop after at most this many seconds (defaul: inf)
    --points=X          begin integration with this lattice size (default: 1e4)
    --presamples=X      use this many points for presampling (default: 1e4)
    --shifts=X          use this many lattice shifts per integral (default: 32)
    --cluster=X         use this cluster.json file
    --coefficients=X    use coefficients from this directory
    --lattice_candidates=X    number of  generating vector candidates for median Qmc rule. (default: -11)
                        If lattice_candidates=0: Use predefined generating vectors
                        If lattice_candidates<0: Use predefined generating vectors by default, but construct larger lattices using median Qmc rule, if required.
    --help              show this help message
Arguments:
    <var>=X             set this integral or coefficient variable to a given value
    int-<var>=X         set this integral variable to a given value
    coeff-<var>=X       set this coefficient variable to a given value
"""

import asyncio
import base64
import collections
import getopt
import json
import math
import numpy as np
import os
import pickle
import random
import re
import subprocess
import sympy as sp
import sys
import tempfile
import time
import traceback

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

def abs2(x):
    return np.real(x)**2 + np.imag(x)**2

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

# Generic scheduling

class RandomScheduler:
    def __init__(self):
        self.workers = []
        self.wspeed = []

    def add_worker(self, worker):
        self.workers.append(worker)
        self.wspeed.append(worker.speed)

    def queue_size(self):
        return sum(w.queue_size() for w in self.workers)

    def call(self, method, *args):
        w = min(random.choices(self.workers, weights=self.wspeed, k=3),
                key=lambda w: w.queue_size())#/w.speed)
        return w.call(method, *args)

    def call_cb(self, method, args, callback, callback_args):
        w = min(random.choices(self.workers, weights=self.wspeed, k=3),
                key=lambda w: w.queue_size())#/w.speed)
        return (w, w.call_cb(method, args, callback, callback_args))

    def cancel_cb(self, token):
        w, t = token
        w.cancel_cb(t)

# Main

async def benchmark_worker(w):
    lattice, genvec = generating_vector(2, 10**3)
    shift = ([0.3, 0.8])
    deformp = ([1.0, 1.0])
    # Measure round-trip latency.
    latency = []
    await w.call("ping")
    for i in range(4):
        t0 = time.time()
        await w.call("ping")
        latency.append(time.time() -  t0)
    w.latency = latency = np.mean(latency)
    # Calculate worker's total per-message overhead by how many
    # empty jobs can it do per second.
    t0 = time.time()
    bench0 = await asyncio.gather(*[
        w.call("integrate", 0, lattice, 0, 1, genvec, shift, deformp)
        for i in range(1000)
    ])
    t1 = time.time()
    dt0 = min(dt for v, dn, dt in bench0)
    w.int_overhead = dt0
    w.overhead = max(w.int_overhead, (t1-t0-latency)/len(bench0))
    # Give the worker a big enough job so that the linear part of
    # the scaling would dominate over the integration overhead.
    # Figure out FLOPS that way.
    for k in (6,7,8,9):
        lattice, genvec = generating_vector(2, 10**k)
        v, dn, dt = await w.call("integrate", 0,
                lattice, 0, lattice, genvec, shift, deformp)
        if dt > dt0*1000:
            break
    w.speed = dn/(dt - dt0)

def bracket_mul(br1, br2, maxorders):
    """
    Multiply two multivariate polynomials represented in the form
    of {exponent: coefficient} dictionaries; skip any terms with
    the combined exponent higher than the specified maximum.
    """
    br = {}
    for k1, v1 in br1.items():
        for k2, v2 in br2.items():
            key = tuple(x + y for x, y in zip(k1, k2))
            if all(x <= y for x, y in zip(key, maxorders)):
                br[key] = br.get(key, 0) + v1*v2
    return br

def adjust_1d_n(W2, V, w, a, tau, nmin, nmax, lattice_candidates):
    assert np.all(W2 > 0)
    assert np.all(w > 0)
    assert np.all(tau > 0)
    assert np.all(nmin > 0)
    if W2 @ (w/nmin**a) <= V:
        return nmin
    n = (w/tau * W2)**(1/(a+1))
    assert not np.any(np.isnan(n))
    n *= (1/V * (W2 @ (w/n**a)))**(1/a)
    assert not np.any(np.isinf(n))
    assert not np.any(np.isnan(n))
    # if not using medianQMC: Enforce nmax, raising the rest
    if lattice_candidates != 0: return n
    mask = (n > nmax)
    while True:
        n[mask] = nmax[mask]
        maskc = np.count_nonzero(mask)
        if maskc == 0: break
        if maskc == len(n): return n
        VV = V - (W2[mask] @ (w[mask]/n[mask]**a))
        if VV < 0:
            log(f"Probably can't reach the target error now that {np.count_nonzero(mask)} integrals are capped at maximum lattice size")
            log(f"... will still try though")
            return n
        assert np.all(VV > 0)
        mask2 = (~mask)
        n[mask2] *= (1/VV * (W2[mask2] @ (w[mask2]/n[mask2]**a)))**(1/a)
        assert not np.any(np.isinf(n))
        assert not np.any(np.isnan(n))
        add = (n > nmax)
        if not np.any(add): break
        mask |= add
    mask = (n < nmin)
    while True:
        n[mask] = nmin[mask]
        maskc = np.count_nonzero(mask)
        if maskc == 0: break
        if maskc == len(n): return n
        VV = np.clip(V - W2[mask] @ (w[mask]/n[mask]**a), 0, np.inf)
        mask2 = (~mask)
        n[mask2] *= (1/VV * (W2[mask2] @ (w[mask2]/n[mask2]**a)))**(1/a)
        add = (n < nmin)
        if not np.any(add): break
        mask |= add
    return n

def adjust_n(W2, V, w, a, tau, nmin, nmax, lattice_candidates, names=[]):
    assert np.all(V>0)
    assert len(W2) == len(V)
    n = nmin.copy()
    for i in range(len(W2)-1, -1, -1):
        mask = (w != 0) & (W2[i,:] != 0)
        n[mask] = adjust_1d_n(W2[i,mask], V[i], w[mask], a, tau[mask], n[mask], nmax[mask], lattice_candidates)
        assert not np.any(np.isnan(n))
    return n

async def doeval(workers, datadir, coeffsdir, intfile, epsabs, epsrel, npresample, npoints0, nshifts, lattice_candidates, valuemap, valuemap_coeff, deadline):
    # Load the integrals from the requested json file
    t0 = time.time()

    with open(intfile, "r") as f:
        info = json.load(f)

    if info["type"] not in ("integral", "sum"):
        raise ValueError(f"unknown type: {info['type']}")

    for p in info["realp"] + info["complexp"]:
        if p not in valuemap:
            raise ValueError(f"missing integral parameter: {p}")

    sp_regulators = sp.var(info["regulators"])
    requested_orders = info["requested_orders"]
    kernel2idx = {}
    valuemap_rat = {
        k:sp.nsimplify(v, rational=True, tolerance=np.abs(v)*1e-13)
        for k, v in valuemap_coeff.items()
    }
    if info["type"] == "integral":
        infos = {info["name"] : info}
        for k in info["kernels"]:
            kernel2idx[info["name"], k] = len(kernel2idx)
        log(f"got the total of {len(kernel2idx)} kernels")
        ampcount = 1
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
        ampcount = len(info["sums"])
    else:
        raise ValueError(f"unknown type: {info['type']}")
    korders = {}
    for fam, ii in infos.items():
        for i, oo in enumerate(ii["orders"]):
            for ker in oo["kernels"]:
                korders.setdefault((fam, ker), i)
    
    log("Kernel ids:")
    for (fam, ker), i in kernel2idx.items():
        log(f"- ({fam}, {ker}) = k{i}")

    realp = {
        i : [valuemap[p] for p in info["realp"]]
        for i, info in infos.items()
    }
    complexp = {
        i : [(np.real(valuemap[p]), np.imag(valuemap[p])) for p in info["complexp"]]
        for i, info in infos.items()
    }
    family2idx = {fam:i for i, fam in enumerate(infos.keys())}

    log(f"parsing {len(infos)} integral prefactors")
    for ii in infos.values():
        ii["expanded_prefactor_value"] = {
            tuple(t["regulator_powers"]) : complex(sp.sympify(t["coefficient"]).subs(valuemap))
            for t in ii["expanded_prefactor"]
        }

    # Launch all the workers
    t1 = time.time()

    par = RandomScheduler()

    async def add_worker(cmd):
        w = await launch_worker(cmd, datadir)
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
        await benchmark_worker(w)
        par.add_worker(w)
    await asyncio.gather(*[add_worker(cmd) for cmd in workers])
    log("workers:")
    for w in par.workers:
        log(f"- {w.name}: int speed={w.speed:.2e}bps, int overhead={w.int_overhead:.2e}s, total overhead={w.overhead:.2e}s, latency={w.latency:.2e}s")

    # Load the integral coefficients
    ap2coeffs = {} # (ampid, powerlist) -> coeflist
    if info["type"] == "integral":
        br_coef = {(0,)*len(info["regulators"]): sp.sympify(1)}
        split_integral_into_orders(ap2coeffs, 0, kernel2idx, info, br_coef, valuemap, sp_regulators, requested_orders)
    elif info["type"] == "sum":
        log("loading amplitude coefficients")
        done_evalf = asyncio.Future()
        done_evalf.todo = sum(len(terms) for terms in info["sums"])
        def evalf_cb(br_coef, exception, w, a, t):
            if exception is not None:
                done_evalf.set_exception(WorkerException(exception))
                return
            log("-", t["coefficient"])
            br_coef = {tuple(k):complex(re, im) for k, (re, im) in br_coef}
            split_integral_into_orders(ap2coeffs, a, kernel2idx, infos[t["integral"]], br_coef, valuemap, sp_regulators, requested_orders)
            done_evalf.todo -= 1
            if done_evalf.todo == 0:
                done_evalf.set_result(None)
        for a, terms in enumerate(info["sums"]):
            for t in terms:
                intinfo = infos[t["integral"]]
                pref_lord = np.min([o["regulator_powers"] for o in intinfo["expanded_prefactor"]], axis=0)
                kern_lord = np.min([o["regulator_powers"] for o in intinfo["orders"]], axis=0)
                coef_ord = - kern_lord - pref_lord + requested_orders
                par.call_cb("evalf", (
                        os.path.relpath(os.path.join(coeffsdir, t["coefficient"]), datadir),
                        {k:str(v) for k,v in valuemap_rat.items()},
                        [[str(var), int(order)] for var, order in zip(sp_regulators, coef_ord)]
                    ),
                    evalf_cb,
                    (a, t)
                )
        await done_evalf

    W = np.stack([w for w in ap2coeffs.values()])
    W2 = abs2(W)
    log(f"will consider {len(ap2coeffs)} sums:")
    for a, p in sorted(ap2coeffs.keys()):
        log(f"- amp{a},", " ".join(f"{r}^{e}" for r, e in zip(sp_regulators, p)))

    epsrel = [epsrel[a] if a < len(epsrel) else epsrel[-1] for a, p in ap2coeffs.keys()]
    epsabsOrig = epsabs.copy()
    epsabs = [epsabs[a] if a < len(epsabs) else epsabs[-1] for a, p in ap2coeffs.keys()]

    # Presample all kernels
    kern_rng = [np.random.RandomState(0) for fam, ker in kernel2idx.keys()]
    t2 = time.time()
    log("distributing presampling jobs")
    results = []
    for i, (fam, ker) in enumerate(kernel2idx.keys()):
        lattice, genvec = generating_vector(infos[fam]["dimension"], npresample)
        f = par.call("maxdeformp", i+1, infos[fam]["deformp_count"],
            lattice, genvec, kern_rng[i].rand(infos[fam]["dimension"]).tolist())
        results.append(f)
    log("waiting for the presampling results")
    deformp = await asyncio.gather(*results)
    deformp = [[min(max(x, 1e-6), 1.0) for x in defp] for defp in deformp]
    for i, d in enumerate(deformp):
        log(f"maxdeformp of k{i} is {d}")

    # Integrate the weighted sum
    t3 = time.time()

    fams = [fam for fam, ker in kernel2idx.keys()]
    kers = [ker for fam, ker in kernel2idx.keys()]
    dims = [infos[fam]["dimension"] for fam in fams]
    genvecs = [None] * len(kernel2idx)
    oldlattices = np.zeros(len(kernel2idx), dtype=np.float64)
    maxlattices = np.array([max_lattice_size(d) for d in dims], dtype=np.float64)
    lattices = np.zeros(len(kernel2idx), dtype=np.float64)
    for i in range(len(kernel2idx)):
        lattices[i], genvecs[i] = generating_vector(dims[i], npoints0)
    shift_val = np.full((len(kernel2idx), nshifts), np.nan, dtype=np.complex128)
    shift_rnd = np.empty((len(kernel2idx), nshifts), dtype=object)
    shift_tag = np.full((len(kernel2idx), nshifts), None, dtype=object)
    kern_db = np.ones(len(kernel2idx))
    kern_dt = np.ones(len(kernel2idx))
    kern_di = np.ones(len(kernel2idx))
    kern_val = np.zeros(len(kernel2idx), dtype=np.complex128)
    kern_var = np.full(len(kernel2idx), np.inf, dtype=np.complex128)

    genvec_candidates = dict()

    def shift_done_cb(result, exception, w, idx, shift):
        (re, im), di, dt = result
        if math.isnan(re) or math.isnan(im):
            for s in range(nshifts):
                par.cancel_cb(shift_tag[idx, s])
            deformp[idx] = tuple(p*0.9 for p in deformp[idx])
            log(f"got NaN from k{idx}; decreasing deformp by 0.9 to {deformp[idx]}")
            schedule_kernel(idx)
        else:
            shift_val[idx, shift] = complex(re, im)
            if dt > 2*w.int_overhead:
                kern_db[idx] += (dt - w.int_overhead)*w.speed
                kern_di[idx] += di
                kern_dt[idx] += dt

    def schedule_kernel(idx):
        for s in range(nshifts):
            shift = kern_rng[idx].rand(dims[idx])
            shift_rnd[idx, s] = shift
            shift_tag[idx, s] = par.call_cb("integrate",
                (idx+1, int(lattices[idx]), 0, int(lattices[idx]), genvecs[idx],
                shift.tolist(),
                deformp[idx]),
                shift_done_cb, (idx, s))

    def shift_done_cb_median_lattice(result, exception, w, idx, shift):
        (re, im), di, dt = result
        if math.isnan(re) or math.isnan(im):
            for s in range(abs(lattice_candidates)):
                par.cancel_cb(shift_tag[idx, s])
            deformp[idx] = tuple(p*0.9 for p in deformp[idx])
            log(f"got NaN from k{idx}; decreasing deformp by 0.9 to {deformp[idx]}")
            schedule_kernel_median_lattice(idx)
        else:
            shift_val[idx, shift] = complex(re, im)
            if dt > 2*w.int_overhead:
                kern_db[idx] += (dt - w.int_overhead)*w.speed
                kern_di[idx] += di
                kern_dt[idx] += dt

    def schedule_kernel_median_lattice(idx):
        for s in range(abs(lattice_candidates)):
            shift = kern_rng[idx].rand(dims[idx])
            shift_rnd[idx, s] = shift
            def rand():
                r = 0
                while math.gcd(int(lattices[idx]),r) != 1:
                    r = kern_rng[idx].randint(1,(lattices[idx])-1)
                return r
            genvec_candidates[(idx, s)] = tuple( rand() for _ in range(dims[idx]) )
            shift_tag[idx, s] = par.call_cb("integrate",
                (idx+1, int(lattices[idx]), 0, int(lattices[idx]), genvec_candidates[(idx,s)],
                shift.tolist(),
                deformp[idx]),
                shift_done_cb_median_lattice, (idx, s))

    perkern_epsrel = 0.2
    perkern_epsabs = 1e-4

    def propose_lattices1(amp_val, amp_var):
        scaling = 2
        K = 20
        kern_maxvar = np.maximum(perkern_epsabs**2, abs2(kern_val)*perkern_epsrel**2)
        kern_absvar = np.real(kern_var) + np.imag(kern_var)
        if np.all(kern_absvar <= kern_maxvar):
            log("per-integral precision reached")
            oldlattices[:] = lattices
            return None
        n = lattices * (kern_absvar/kern_maxvar)**(1/scaling)
        n = np.clip(n, lattices, lattices*K)
        mask_toolo = (lattices < n) & (n < lattices * 2)
        n[mask_toolo] = lattices[mask_toolo]*2
        return n

    def propose_lattices2(amp_val, amp_var):
        scaling = 2
        K = 20
        amp_absval = np.sqrt(abs2(amp_val))
        amp_abserr = np.sqrt(np.real(amp_var) + np.imag(amp_var))
        log("absval =", amp_absval)
        log("abserr =", amp_abserr)
        amp_relerr = amp_abserr/amp_absval
        amp_relerr[amp_abserr == 0] = 0
        log("relerr =", amp_relerr)
        amp_maxerr = np.maximum(epsabs, amp_absval*epsrel)
        log("max abserr =", amp_maxerr)
        log("abserr K =", amp_abserr/amp_maxerr)
        if np.all(amp_abserr <= amp_maxerr):
            log("per-amplitude precision reached")
            oldlattices[:] = lattices
            return None
        tau = kern_db/kern_di
        kern_absvar = np.real(kern_var) + np.imag(kern_var)
        v0 = kern_absvar * lattices**scaling
        n = adjust_n(W2, amp_maxerr**2, v0, scaling, tau, lattices, maxlattices, lattice_candidates)
        n = np.clip(n, lattices, lattices*K)
        toobig = n >= lattices*K
        if np.any(toobig):
            n[toobig] = lattices[toobig]*K
            toosmall = n < lattices*2
            n[toosmall] = lattices[toosmall]
        return n

    early_exit = False

    async def iterate_integration(propose_lattices):
        nonlocal early_exit
        while True:
            mask_todo = lattices != oldlattices
            if np.any(mask_todo):
                if lattice_candidates != 0:
                    for i in mask_todo.nonzero()[0]:
                        if(lattice_candidates > 0 or lattices[i] > maxlattices[i]):
                            schedule_kernel_median_lattice(int(i))
                            await asyncio.sleep(0)
                    if par.queue_size() > 0:
                        log(f"distributing {par.queue_size()} jobs to construct generating vectors")
                        log("worker queue sizes:")
                        for w in par.workers:
                            log(f"- {w.name}: {w.queue_size()}")
                    while par.queue_size() > 0:
                        log(f"still todo: {par.queue_size()}")
                        early_exit = time.time() > deadline
                        if early_exit:
                            log("WARNING: timeout reached, will stop soon")
                            break
                        await asyncio.sleep(1)
                    def signedMax(x):
                        if isinstance(x,complex):
                            return x.real if abs(x.real) > abs(x.imag) else x.imag
                        else: return x
                    for i in mask_todo.nonzero()[0]:
                        if(lattice_candidates > 0 or lattices[i] > maxlattices[i]):
                            median = np.median([signedMax(x) for x in shift_val[i,:abs(lattice_candidates)]])
                            for s in range(abs(lattice_candidates)):
                                if signedMax(shift_val[i,s]) == median:
                                    genvecs[i] = list(genvec_candidates[(i,s)])
                                shift_val[i,s] = np.nan

                log(f"distributing {np.count_nonzero(mask_todo)}*{nshifts} integration jobs")
                for i in mask_todo.nonzero()[0]:
                    schedule_kernel(int(i))
                    await asyncio.sleep(0)
                log("worker queue sizes:")
                for w in par.workers:
                    log(f"- {w.name}: {w.queue_size()}")
                while par.queue_size() > 0:
                    log(f"still todo: {par.queue_size()}")
                    early_exit = time.time() > deadline
                    if early_exit:
                        log("WARNING: timeout reached, will stop soon")
                        break
                    await asyncio.sleep(1)
                mask_done = np.logical_and(mask_todo, ~np.any(np.isnan(shift_val), axis=1))
                log(f"integration done, updated {np.count_nonzero(mask_done)} kernels")
                shift_val_done = shift_val[mask_done]
                shift_val[mask_todo] = np.nan
                new_kern_val = np.mean(shift_val_done, axis=1)
                new_kern_val /= lattices[mask_done]
                new_kern_var = np.var(np.real(shift_val_done), axis=1) + (1j)*np.var(np.imag(shift_val_done), axis=1)
                new_kern_var /= lattices[mask_done]**2 * nshifts

                latticex = lattices/oldlattices
                precisionx = np.sqrt((np.real(kern_var[mask_done]) + np.imag(kern_var[mask_done])) / (np.real(new_kern_var) + np.imag(new_kern_var)))
                for i, idx in enumerate(mask_done.nonzero()[0]):
                    idx = int(idx)
                    if precisionx[i] < 1.0:
                        log(f"k{idx} @ {lattices[i]:.3e} = {new_kern_val[i]:.16e} ~ {new_kern_var[i]:.3e} ({1/precisionx[i]:.4g}x worse at {latticex[i]:.1f}x lattice)")
                    else:
                        log(f"k{idx} @ {lattices[i]:.3e} = {new_kern_val[i]:.16e} ~ {new_kern_var[i]:.3e} ({precisionx[i]:.4g}x better at {latticex[i]:.1f}x lattice)")
                submask_lucky = new_kern_var <= kern_var[mask_done]
                kern_val[mask_done] = np.where(submask_lucky, new_kern_val, kern_val[mask_done])
                kern_var[mask_done] = np.where(submask_lucky, new_kern_var, kern_var[mask_done])
                log(f"unlucky results: {np.count_nonzero(~submask_lucky)} out of {np.count_nonzero(mask_done)}")
            amp_val = W @ kern_val
            amp_var = W2 @ kern_var
            # Report results
            if np.any(mask_todo):
                ampids = sorted(set(a for a, p in ap2coeffs.keys()))
                relerrs = []
                for ampid in ampids:
                    relerr = []
                    log(f"amp{ampid}=(")
                    for (a, p), val, var in sorted(zip(ap2coeffs.keys(), amp_val, amp_var)):
                        if a != ampid: continue
                        stem = "*".join(f"{r}^{p}" for r, p in zip(sp_regulators, p))
                        err = np.sqrt(np.real(var)) + (1j)*np.sqrt(np.imag(var))
                        log(f"  +{stem}*({val:+.16e})")
                        log(f"  +{stem}*({err:+.16e})*plusminus")
                        abserr = np.abs(err)
                        relerr.append(abserr / np.abs(val) if abserr > epsabsOrig[ampid] else 0.0)
                    log(")")
                    log(f"amp{ampid} relative errors by order:", ", ".join(f"{e:.2e}" for e in relerr))
                    relerrs.append(np.max(relerr))
                log(f"largest relative error: {np.max(relerrs):.2e} (amp{np.argmax(relerrs)})")

            if early_exit:
                return amp_val, amp_var

            n = propose_lattices(amp_val, amp_var)
            if n is None:
                return amp_val, amp_var
            newgenvecs = [None] * len(kernel2idx)
            if lattice_candidates <= 0: 
                for i in range(len(kernel2idx)):
                    try: n[i], newgenvecs[i] = generating_vector(dims[i], n[i])
                    except ValueError: 
                        if lattice_candidates < 0: pass
            if not np.any(n != lattices):
                log("can't increase the lattice sizes any more; giving up")
                return amp_val, amp_var
            for i, (l1, l2) in enumerate(zip(lattices, n)):
                if l1 != l2:
                    log(f"lattice[k{i}] = {l1:.0f} -> {l2:.0f} ({l2/l1:.1f}x)")
            oldlattices[:] = lattices
            lattices[:] = n
            genvecs[:] = newgenvecs

    if not early_exit:
        if np.min(epsrel) < 0.1:
            log(f"trying to achieve epsrel={perkern_epsrel} and epsabs={perkern_epsabs} for each kernel")
            amp_val, amp_var = await iterate_integration(propose_lattices1)

    if not early_exit:
        log(f"trying to achieve epsrel={epsrel} and epsabs={epsabs} for each amplitude")
        amp_val, amp_var = await iterate_integration(propose_lattices2)

    # Report the results
    t4 = time.time()
    log("integral load time:", t1-t0)
    log("worker startup time:", t2-t1)
    log("presampling time:", t3-t2)
    log("integration time:", t4-t3)

    log(f"per-integral statistics:")
    for f in set(fams):
        mask = np.array([ff == f for ff in fams])
        dt = np.sum(kern_dt[mask])
        di = np.sum(kern_di[mask])
        slow = kern_db[mask]/kern_di[mask]
        minlattice = np.min(lattices[mask])
        maxlattice = np.max(lattices[mask])
        log(f"- {f}: {di:.4e} evals, {dt:.4e} sec, {np.min(slow):.4g} - {np.max(slow):.4g} bubbles, {minlattice:.4e} - {maxlattice:.4e} pts")

    relerrs = []
    print("[")
    for ampid in range(ampcount):
        relerr = []
        print("  (")
        for (a, p), val, var in sorted(zip(ap2coeffs.keys(), amp_val, amp_var)):
            if a != ampid: continue
            stem = "*".join(f"{r}^{p}" for r, p in zip(sp_regulators, p))
            err = np.sqrt(np.real(var)) + (1j)*np.sqrt(np.imag(var))
            print(f"    +{stem}*({val:+.16e})")
            print(f"    +{stem}*({err:+.16e})*plusminus")
            abserr = np.abs(err)
            #print(abserr, epsabs[ampid], ampid)
            relerr.append(abserr / np.abs(val) if abserr > epsabsOrig[ampid] else 0.0)
        if len(relerr) != 0:
            print("  )," if ampid < ampcount-1 else "  )")
            sys.stdout.flush()
            log(f"amp{ampid} relative errors by order:", ", ".join(f"{e:.2e}" for e in relerr))
            relerrs.append(np.max(relerr))
        else:
            print("    0\n  )," if ampid < ampcount-1 else "    0\n  )")
            sys.stdout.flush()
            log(f"amp{ampid} has no terms")
            relerrs.append(0.0)
    print("]")
    sys.stdout.flush()
    log(f"largest relative error: {np.max(relerrs):.2e} (amp{np.argmax(relerrs)})")

def load_cluster_json(jsonfile, dirname):
    try:
        with open(jsonfile, "r") as f:
            cluster_json = json.load(f)
            assert "cluster" in cluster_json
            assert isinstance(cluster_json["cluster"], list)
            log(f"Using cluster configuration from {jsonfile!r}")
            return cluster_json
    except FileNotFoundError:
        log(f"Can't find {jsonfile}; will run locally")
        pass
    ncpu = 0
    try:
        ncpu = max(1, len(os.sched_getaffinity(0)) - 1)
    except AttributeError:
        ncpu = max(1, os.cpu_count() - 1)
    ncuda = 0
    if os.path.exists(os.path.join(dirname, "builtin.fatbin")):
        try:
            p = subprocess.run([os.path.join(contrib_dirname, "bin", "pysecdec_listcuda")], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            ncuda = len(p.stdout.splitlines())
        except Exception as e:
            log(f"Can't determine GPU count: {e}")
    else:
        log(f"CUDA worker data was not built, skipping")
    if ncuda > 0: ncpu = 0
    log(f"local CPU worker count: {ncpu}, GPU worker count: {ncuda}")
    return {
        "cluster":
            [{"count": ncpu, "command": ["nice", sys.executable, "-m", "pySecDecContrib", "pysecdec_cpuworker"]}] +
            [
                {"count": 1, "command": ["nice", sys.executable, "-m", "pySecDecContrib", "pysecdec_cudaworker", "-d", str(i)]}
                for i in range(ncuda)
            ]
    }

def split_integral_into_orders(orders, ampid, kernel2idx, info, br_coef, valmap, sp_regulators, requested_orders):
    br_pref = info["expanded_prefactor_value"]
    kern_leading_orders = np.min([o["regulator_powers"] for o in info["orders"]], axis=0)
    prefactor = bracket_mul(br_pref, br_coef, -kern_leading_orders + requested_orders)
    for o in info["orders"]:
        powers = np.array(o["regulator_powers"])
        for pow, coef in prefactor.items():
            p = powers + pow
            if np.all(p <= requested_orders):
                key = (ampid, tuple(p.tolist()))
                orders.setdefault(key, np.zeros(len(kernel2idx), dtype=np.complex128))
                for k in o["kernels"]:
                    orders[key][kernel2idx[info["name"], k]] += coef

def parse_array_shorthand(text):
    """
    Parse a shorthand notation for a array:
        array ::= item [ "," item ]*
        item ::= [ repetition-count "x" ]? value
    Example:
        >>> parse_array_shorthand("2x1,3x2,3")
        [1,1,2,2,2,3]
    """
    result = []
    for item in text.split(","):
        if "x" in item:
            repeat, value = item.split("x", 1)
            value = float(value)
            for i in range(int(repeat)):
                result.append(value)
        else:
            result.append(float(item))
    return result

def parse_unit(text, units):
    for unit, magnitude in units.items():
        if text.endswith(unit):
            return float(text[:-len(unit)])*magnitude
    return float(text)

def main():

    valuemap = {}
    valuemap_coeff = {}
    valuemap_int = {}
    npoints = 10**4
    npresamples = 10**4
    epsabs = [1e-10]
    epsrel = [1e-4]
    nshifts = 32
    clusterfile = None
    coeffsdir = None
    lattice_candidates = -11
    deadline = math.inf
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "", ["cluster=", "coefficients=", "epsabs=", "epsrel=", "points=", "presamples=", "shifts=", "timeout=", "lattice_candidates=", "help"])
    except getopt.GetoptError as e:
        print(e, file=sys.stderr)
        print("use --help to see the usage", file=sys.stderr)
        exit(1)
    for key, value in opts:
        if key == "--cluster": clusterfile = value
        elif key == "--coefficients": coeffsdir = value
        elif key == "--epsabs": epsabs = parse_array_shorthand(value)
        elif key == "--epsrel": epsrel = parse_array_shorthand(value)
        elif key == "--points": npoints = int(float(value))
        elif key == "--presamples": npresamples = int(float(value))
        elif key == "--shifts": nshifts = int(float(value))
        elif key == "--timeout": deadline = time.time() + parse_unit(value, {"s": 1, "m": 60, "h": 60*60, "d": 24*60*60})
        elif key == "--lattice_candidates": lattice_candidates = int(float(value))
        elif key == "--help":
            print(__doc__.strip())
            exit(0)
    if len(args) < 1:
        print(__doc__.strip(), file=sys.stderr)
        exit(1)
    intfile = args[0]
    dirname = os.path.dirname(intfile)
    if coeffsdir is None: coeffsdir = os.path.join(dirname, "coefficients")
    clusterfile = os.path.join(dirname, "cluster.json") if clusterfile is None else clusterfile
    if lattice_candidates!=0 and lattice_candidates % 2 == 0:
        lattice_candidates += 1 if lattice_candidates > 0 else -1
    log("Settings:")
    log(f"- file = {intfile}")
    log(f"- epsabs = {epsabs}")
    log(f"- epsrel = {epsrel}")
    log(f"- points = {npoints}")
    log(f"- presamples = {npresamples}")
    log(f"- shifts = {nshifts}")
    log(f"- lattice_candidates = {lattice_candidates}")
    for arg in args[1:]:
        if "=" not in arg: raise ValueError(f"Bad argument: {arg}")
        key, value = arg.split("=", 1)
        value = complex(value)
        value = value.real if value.imag == 0 else value
        if key.startswith("coeff-"):
            valuemap_coeff[key[6:]] = value
        elif key.startswith("int-"):
            valuemap_int[key[4:]] = value
        else:
            valuemap[key] = value
    valuemap_int = {**valuemap, **valuemap_int}
    valuemap_coeff = {**valuemap, **valuemap_coeff}
    log("Invariants:")
    for key, value in valuemap_int.items():
        if valuemap_coeff.get(key, None) == value:
            log(f"- {key} = {value}")
    for key, value in valuemap_int.items():
        if valuemap_coeff.get(key, None) != value:
            log(f"- integral {key} = {value}")
    for key, value in valuemap_coeff.items():
        if valuemap_int.get(key, None) != value:
            log(f"- coefficient {key} = {value}")

    # Load worker list
    cluster = load_cluster_json(clusterfile, dirname)
    workers = []
    for w in cluster["cluster"]:
        workers.extend([w["command"]] * w.get("count", 1))
    if len(workers) == 0:
        log("No workers defined")
        exit(1)

    # Begin evaluation
    loop = asyncio.get_event_loop()
    loop.run_until_complete(doeval(workers, dirname, coeffsdir, intfile, epsabs, epsrel, npresamples, npoints, nshifts, lattice_candidates, valuemap_int, valuemap_coeff, deadline))

if __name__ == "__main__":
    main()

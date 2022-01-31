#!/usr/bin/env python3
"""
Distributed (and local) evaluator for pySecDec integrals.
Usage:
    python3 -m pySecDec.deval integrand.json [options] param=value ...
Options:
    --epsabs=X          stop if this absolute precision is achieved (default: 1e-10)
    --epsrel=X          stop if this relative precision is achieved (default: 1e-4)
    --points=X          begin integration with this lattice size (default: 1e4)
    --shifts=X          use this many lattice shifts per integral (default: 32)
    --cluster=X         use this cluster.json file
    --coefficients=X    use coefficients from this directory
    --help              show this help message
Arguments:
    param=value use this value for the given integral parameter
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
        p = await asyncio.create_subprocess_shell(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
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
        time.sleep(timeout)
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

def bracket(expr, varlist):
    """
    Collect terms in the SymPy expression by powers of the
    variables; return a dict. Example:
    >>> bracket(sp.sympify("1 + 2*a + 3*b + 4*a*b"), sp.var(["a", "b"]))
    {(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4}
    """
    result = {}
    orders = sp.collect(expr, varlist[0], evaluate=False)
    othervars = varlist[1:]
    for stem, coef in orders.items():
        power = int(stem.exp) if stem.is_Pow else 0 if stem == 1 else 1
        if othervars:
            for powers, c in bracket(coef, othervars).items():
                result[(power,) + powers] = c
        else:
            result[(power,)] = coef
    return result

def ginsh_series(ex, var, order):
    if not ex.has(var):
        return ex
    hashed = {}
    def hashfn(m):
        v = f"hash{len(hashed)}"
        hashed[v] = m.group(0)
        return v
    text = str(ex).replace("**", "^")
    text = re.sub(r"polygamma\(([^)]*)\)", hashfn, text)
    with tempfile.NamedTemporaryFile(prefix="psd_ginsh", mode="w") as f:
        f.write("START;\nseries((")
        f.write(text)
        f.write("),(")
        f.write(str(var))
        f.write("),(")
        f.write(str(order))
        f.write("));\nquit;")
        f.flush()
        subprocess.check_call(f"'{contrib_dirname}/bin/ginsh' '{f.name}' > '{f.name}.out'", shell=True)
        with open(f"{f.name}.out", "r") as f2:
            result = f2.read()
        os.unlink(f"{f.name}.out")
    result = re.sub(r".*START\n", "", result, flags=re.DOTALL)
    result = re.sub(r"[+]Order\([^)]*\)", "", result, flags=re.DOTALL)
    result = result.strip()
    assert result != ""
    return sp.sympify(result).subs(hashed)

def ginsh_symplify(expr, valmap):
    with tempfile.NamedTemporaryFile(prefix="psd_ginsh", mode="w") as f:
        f.write("START;\n");
        for k, v in valmap.items():
            f.write(f"{k} = {v}:\n")
        f.write("(")
        f.write(expr)
        f.write(");\nquit;\n")
        f.flush()
        subprocess.check_call(f"'{contrib_dirname}/bin/ginsh' '{f.name}' > '{f.name}.out'", shell=True)
        with open(f"{f.name}.out", "r") as f2:
            result = f2.read()
        os.unlink(f"{f.name}.out")
    result = re.sub(r".*START\n", "", result, flags=re.DOTALL)
    result = re.sub(r"[+]Order\([^)]*\)", "", result, flags=re.DOTALL)
    result = result.strip()
    assert result != ""
    return sp.sympify(result)

def series_bracket(expr, varlist, orderlist):
    if expr == 0:
        return {}
    result = {}
    orders = sp.collect(ginsh_series(expr, varlist[0], orderlist[0]+1), varlist[0], evaluate=False)
    othervars = varlist[1:]
    otherorders = orderlist[1:]
    for stem, coef in orders.items():
        power = int(stem.exp) if stem.is_Pow else 0 if stem == 1 else 1
        if othervars:
            for powers, c in series_bracket(coef, othervars, otherorders).items():
                result[(power,) + powers] = c
        else:
            result[(power,)] = coef
    return result

def adjust_1d_n(W2, V, w, a, tau, nmin, nmax):
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
    # Enforce nmax, raising the rest
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

def adjust_n(W2, V, w, a, tau, nmin, nmax, names=[]):
    assert np.all(V>0)
    assert len(W2) == len(V)
    n = nmin.copy()
    for i in range(len(W2)-1, -1, -1):
        mask = (w != 0) & (W2[i,:] != 0)
        n[mask] = adjust_1d_n(W2[i,mask], V[i], w[mask], a, tau[mask], n[mask], nmax[mask])
        assert not np.any(np.isnan(n))
    return n

async def doeval(workers, datadir, coeffsdir, intfile, epsabs, epsrel, npresample, npoints0, nshifts, valuemap):
    # Load the integrals from the requested json file
    t0 = time.time()

    with open(intfile, "r") as f:
        info = json.load(f)

    sp_regulators = sp.var(info["regulators"])
    requested_orders = info["requested_orders"]
    ap2coeffs = {} # (ampid, powerlist) -> coeflist
    valuemap_rat = {k:sp.nsimplify(v, rational=True, tolerance=np.abs(v)*1e-13) for k, v in valuemap.items()}
    if info["type"] == "integral":
        infos = {info["name"] : info}
        kernel2idx = {}
        for k in info["kernels"]:
            kernel2idx[info["name"], k] = len(kernel2idx)
        log(f"got the total of {len(kernel2idx)} kernels")
        split_integral_into_orders(ap2coeffs, 0, kernel2idx, info, 1, valuemap_rat, sp_regulators, requested_orders)
    elif info["type"] == "sum":
        log(f"loading {len(info['integrals'])} integrals")
        infos = {}
        kernel2idx = {}
        for i in info["integrals"]:
            with open(os.path.join(datadir, f"{i}.json"), "r") as f:
                infos[i] = json.load(f)
                assert infos[i]["name"] == i
            for k in infos[i]["kernels"]:
                kernel2idx[i, k] = len(kernel2idx)
        log(f"got the total of {len(kernel2idx)} kernels")
        log("loading amplitude coefficients")
        for a, terms in enumerate(info["sums"]):
            for t in terms:
                log("-", t["coefficient"])
                co = load_coefficient(os.path.join(coeffsdir, t["coefficient"]), valuemap_rat)
                split_integral_into_orders(ap2coeffs, a, kernel2idx, infos[t["integral"]], co, valuemap_rat, sp_regulators, requested_orders)
    else:
        raise ValueError(f"unknown type: {info['type']}")
    korders = {}
    for fam, ii in infos.items():
        for i, oo in enumerate(ii["orders"]):
            for ker in oo["kernels"]:
                korders.setdefault((fam, ker), i)

    realp = {
        i : [valuemap[p] for p in info["realp"]]
        for i, info in infos.items()
    }
    complexp = {
        i : [(np.real(valuemap[p]), np.imag(valuemap[p])) for p in info["complexp"]]
        for i, info in infos.items()
    }
    family2idx = {fam:i for i, fam in enumerate(infos.keys())}
    W = np.stack([w for w in ap2coeffs.values()])
    W2 = abs2(W)
    log(f"will consider {len(ap2coeffs)} sums:")
    for a, p in sorted(ap2coeffs.keys()):
        log(f"- amp{a},", " ".join(f"{r}^{e}" for r, e in zip(sp_regulators, p)))

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
    shift_val = np.zeros((len(kernel2idx), nshifts), dtype=np.complex128)
    shift_rnd = np.empty((len(kernel2idx), nshifts), dtype=np.object)
    shift_tag = np.full((len(kernel2idx), nshifts), None, dtype=np.object)
    kern_db = np.ones(len(kernel2idx))
    kern_dt = np.ones(len(kernel2idx))
    kern_di = np.ones(len(kernel2idx))
    kern_val = np.zeros(len(kernel2idx), dtype=np.complex128)
    kern_var = np.full(len(kernel2idx), np.inf, dtype=np.complex128)

    def shift_done_cb(result, exception, w, idx, shift):
        (re, im), di, dt = result
        if math.isnan(re) or math.isnan(im):
            for s in range(nshifts):
                par.cancel_cb(shift_tag[idx, s])
            deformp[idx] = tuple(p*0.9 for p in deformp[idx])
            log(f"got NaN from k{idx}; decreasing deformp by 0.9")
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
        n = adjust_n(W2, amp_maxerr**2, v0, scaling, tau, lattices, maxlattices)
        n = np.clip(n, lattices, lattices*K)
        toobig = n >= lattices*K
        if np.any(toobig):
            n[toobig] = lattices[toobig]*K
            toosmall = n < lattices*2
            n[toosmall] = lattices[toosmall]
        return n

    async def iterate_integration(propose_lattices):
        while True:
            mask_todo = lattices != oldlattices
            if np.any(mask_todo):
                log(f"distributing {np.count_nonzero(mask_todo)*nshifts} integration jobs")
                for i in mask_todo.nonzero()[0]:
                    schedule_kernel(int(i))
                    asyncio.sleep(0)
                log("worker queue sizes:")
                for w in par.workers:
                    log(f"- {w.name}: {w.queue_size()}")
                while par.queue_size() > 0:
                    log(f"still todo: {par.queue_size()}")
                    await asyncio.sleep(1)
                log("integration done")
                new_kern_val = np.mean(shift_val, axis=1)
                new_kern_val /= lattices
                new_kern_var = np.var(np.real(shift_val), axis=1) + (1j)*np.var(np.imag(shift_val), axis=1)
                new_kern_var /= lattices**2 * nshifts

                latticex = lattices/oldlattices
                precisionx = np.sqrt((np.real(kern_var) + np.imag(kern_var)) / (np.real(new_kern_var) + np.imag(new_kern_var)))
                for i in mask_todo.nonzero()[0]:
                    i = int(i)
                    if precisionx[i] < 1.0:
                        log(f"k{i} @ {lattices[i]:.3e} = {new_kern_val[i]:.16e} ~ {new_kern_var[i]:.3e} ({1/precisionx[i]:.4g}x worse at {latticex[i]:.1f}x lattice)")
                    else:
                        log(f"k{i} @ {lattices[i]:.3e} = {new_kern_val[i]:.16e} ~ {new_kern_var[i]:.3e} ({precisionx[i]:.4g}x better at {latticex[i]:.1f}x lattice)")
                mask_lucky = np.logical_and(new_kern_var <= kern_var, mask_todo)
                kern_val[mask_lucky] = new_kern_val[mask_lucky]
                kern_var[mask_lucky] = new_kern_var[mask_lucky]
                mask_unlucky = np.logical_and(~mask_lucky, mask_todo)
                log(f"unlucky results: {np.count_nonzero(mask_unlucky)} out of {np.count_nonzero(mask_todo)}")

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
                        relerr.append(abserr / np.abs(val) if abserr > epsabs else 0.0)
                    log(")")
                    log(f"amp{ampid} relative errors by order:", ", ".join(f"{e:.2e}" for e in relerr))
                    relerrs.append(np.max(relerr))
                log(f"largest relative error: {np.max(relerrs):.2e} (amp{np.argmax(relerrs)})")
            n = propose_lattices(amp_val, amp_var)
            if n is None:
                return amp_val, amp_var
            newlattices = np.zeros(len(kernel2idx), dtype=np.float64)
            newgenvecs = [None] * len(kernel2idx)
            for i in range(len(kernel2idx)):
                newlattices[i], newgenvecs[i] = generating_vector(dims[i], n[i])
            if not np.any(newlattices != lattices):
                log("can't increase the lattice sizes any more; giving up")
                return amp_val, amp_var
            for i, (l1, l2) in enumerate(zip(lattices, newlattices)):
                if l1 != l2:
                    log(f"lattice[k{i}] = {l1:.0f} -> {l2:.0f} ({l2/l1:.1f}x)")
            oldlattices[:] = lattices
            lattices[:] = newlattices
            genvecs[:] = newgenvecs

    if epsrel < 0.1:
        log(f"trying to achieve epsrel={perkern_epsrel} and epsabs={perkern_epsabs} for each kernel")
        await iterate_integration(propose_lattices1)

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

    ampids = sorted(set(a for a, p in ap2coeffs.keys()))
    relerrs = []
    print("[")
    for i, ampid in enumerate(ampids):
        relerr = []
        print("  (")
        for (a, p), val, var in sorted(zip(ap2coeffs.keys(), amp_val, amp_var)):
            if a != ampid: continue
            stem = "*".join(f"{r}^{p}" for r, p in zip(sp_regulators, p))
            err = np.sqrt(np.real(var)) + (1j)*np.sqrt(np.imag(var))
            print(f"    +{stem}*({val:+.16e})")
            print(f"    +{stem}*({err:+.16e})*plusminus")
            abserr = np.abs(err)
            relerr.append(abserr / np.abs(val) if abserr > epsabs else 0.0)
        print("  )," if i < len(ampids)-1 else "  )")
        sys.stdout.flush()
        log(f"amp{ampid} relative errors by order:", ", ".join(f"{e:.2e}" for e in relerr))
        relerrs.append(np.max(relerr))
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
            [{"count": ncpu, "command": f"nice python3 -m pySecDecContrib pysecdec_cpuworker"}] +
            [
                {"count": 1, "command": f"nice python3 -m pySecDecContrib pysecdec_cudaworker -d {i}"}
                for i in range(ncuda)
            ]
    }

def load_coefficient(filename, valmap):
    tr = {ord(" "): None, ord("\n"): None, ord("\\"): None}
    coeff = sp.sympify(1)
    with open(filename, "r") as f:
        text = f.read().translate(tr)
    for part in text.split(";", 3):
        part = part.strip()
        if not part: continue
        key, value = part.split("=")
        key = key.strip()
        value = ginsh_symplify(value, valmap)
        if key == "numerator": coeff *= value
        if key == "denominator": coeff /= value
        if key == "regulator_factor": coeff *= value
    return coeff

def split_integral_into_orders(orders, ampid, kernel2idx, info, coefficient, valmap, sp_regulators, requested_orders):
    highest_orders = np.min([o["regulator_powers"] for o in info["orders"]], axis=0)
    prefactor = series_bracket(coefficient*sp.sympify(info["prefactor"]).subs(valmap), sp_regulators, -highest_orders + requested_orders)
    prefactor = {p : complex(c) for p, c in prefactor.items()}
    oo = {}
    for o in info["orders"]:
        powers = np.array(o["regulator_powers"])
        for pow, coef in prefactor.items():
            p = powers + pow
            c = complex(coef)
            if np.all(p <= requested_orders):
                key = (ampid, tuple(p.tolist()))
                orders.setdefault(key, np.zeros(len(kernel2idx), dtype=np.complex128))
                oo.setdefault(key, np.zeros(len(kernel2idx), dtype=np.complex128))
                for k in o["kernels"]:
                    orders[key][kernel2idx[info["name"], k]] += coef
                    oo[key][kernel2idx[info["name"], k]] += coef

def main():

    valuemap = {}
    npoints = 10**4
    epsabs = 1e-10
    epsrel = 1e-4
    nshifts = 32
    clusterfile = None
    coeffsdir = None
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "", ["cluster=", "coefficients=", "epsabs=", "epsrel=", "points=", "shifts=", "help"])
    except getopt.GetoptError as e:
        print(e, file=sys.stderr)
        print("use --help to see the usage", file=sys.stderr)
        exit(1)
    for key, value in opts:
        if key == "--cluster": clusterfile = value
        elif key == "--coefficients": coeffsdir = value
        elif key == "--epsabs": epsabs = float(value)
        elif key == "--epsrel": epsrel = float(value)
        elif key == "--points": npoints = int(float(value))
        elif key == "--shifts": nshifts = int(float(value))
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
    log("Settings:")
    log(f"- file = {intfile}")
    log(f"- epsabs = {epsabs}")
    log(f"- epsrel = {epsrel}")
    log(f"- points = {npoints}")
    log(f"- shifts = {nshifts}")
    log("Invariants:")
    for arg in args[1:]:
        if "=" not in arg: raise ValueError(f"Bad argument: {arg}")
        key, value = arg.split("=", 1)
        value = complex(value)
        value = value.real if value.imag == 0 else value
        valuemap[key] = value
        log(f"- {key} = {value}")

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
    loop.run_until_complete(doeval(workers, dirname, coeffsdir, intfile, epsabs, epsrel, npoints, npoints, nshifts, valuemap))

if __name__ == "__main__":
    main()

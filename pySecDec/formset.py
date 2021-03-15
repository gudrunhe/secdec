#!/bin/sh
""":" .

exec python "$0" "$@"
"""

# MIT License
#
# Copyright (c) 2021 Takahiro Ueda
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from __future__ import print_function

import argparse
import contextlib
import copy
import math
import os
import re
import subprocess
import sys

try:
    from typing import TYPE_CHECKING
except ImportError:
    TYPE_CHECKING = False

if TYPE_CHECKING:
    from typing import (
        Any,
        Dict,
        Iterator,
        List,
        Optional,
        TextIO,
        Tuple,
        Union,
        overload,
    )

    from typing_extensions import Literal
else:

    def overload(f):  # noqa: D103
        return None


__doc__ = """\
Generate form.set suited for the local machine.

Example
-------
$ formset.py -o
$ tform `formset.py -f` calcdia.frm
$ minos `formset.py -m` minos.file

Python versions
---------------
2.7, 3.2, 3.3, 3.4, 3.5
"""


if "check_output" not in dir(subprocess):
    # For old systems where Python 2.6 + argparse available.
    def check_output(*popenargs, **kwargs):  # type: ignore[no-untyped-def]
        """Run a command."""
        if "stdout" in kwargs:  # pragma: no cover
            raise ValueError("stdout argument not allowed, " "it will be overridden.")
        process = subprocess.Popen(  # type: ignore[call-overload]  # noqa: E501,S603
            stdout=subprocess.PIPE, *popenargs, **kwargs
        )
        output, _ = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            # `output` keyword is not available in 2.6.
            raise subprocess.CalledProcessError(retcode, cmd)
        return output

    subprocess.check_output = check_output


@contextlib.contextmanager
def open_w_or_stdout(filename=None):
    # type: (Optional[str]) -> Iterator[TextIO]
    """Context manager for a file or stdout."""
    if filename:
        # See https://stackoverflow.com/a/2333979.
        tmpfilename = "{0}.tmp{1}".format(filename, os.getpid())
        f = open(tmpfilename, "w")
        try:
            yield f
        finally:
            f.flush()
            os.fsync(f.fileno())
            f.close()
            os.rename(tmpfilename, filename)
    else:
        yield sys.stdout


def round_down(x, n):
    # type: (int, int) -> int
    """Round down `x` to nearest `n`."""
    return x // n * n


def round_up(x, n):
    # type: (int, int) -> int
    """Round up `x` to nearest `n`."""
    return (x + (n - 1)) // n * n


def metric_prefix(s):
    # type: (str) -> int
    """Parse a metric prefix as a number."""
    s_old = s
    s = s.strip().lower()
    if s == "":
        return 1
    if s == "k":
        return 1000
    if s == "m":
        return 1000 ** 2
    if s == "g":
        return 1000 ** 3
    if s == "t":
        return 1000 ** 4
    raise ValueError("unknown metric prefix: {0}".format(s_old))


def parse_number(s):
    # type: (str) -> int
    """Parse a string as a number with a possible metric prefix."""
    scale = 1
    m = re.match(r"(.*)([kmgtKMGT])$", s)
    if m:
        s = m.group(1)
        scale = metric_prefix(m.group(2))
    # May raise ValueError for bad `s`.
    return int(float(s) * scale)


@overload
def round_human_readable(x, up, tostring):  # noqa: D103
    # type: (int, bool, Literal[True]) -> str
    pass


@overload  # noqa: F811
def round_human_readable(x, up, tostring):  # noqa: D103, F811
    # type: (int, bool, Literal[False]) -> int
    pass


def round_human_readable(x, up, tostring):  # noqa: F811
    # type: (int, bool, bool) -> Union[int, str]
    """Round off `x` within a human readable form."""
    round_off = round_up if up else round_down
    # Take 3 significant figures.
    n = 10 ** (int(math.floor(math.log10(x))) - 2)
    x = round_off(x, n)
    # Find a good suffix which doesn't change the value.
    xx = round_off(x, 1000 ** 4)
    if xx == x:
        return "{0}T".format(xx // 1000 ** 4) if tostring else xx
    xx = round_off(x, 1000 ** 3)
    if xx == x:
        return "{0}G".format(xx // 1000 ** 3) if tostring else xx
    xx = round_off(x, 1000 ** 2)
    if xx == x:
        return "{0}M".format(xx // 1000 ** 2) if tostring else xx
    xx = round_off(x, 1000)
    if xx == x:
        return "{0}K".format(xx // 1000) if tostring else xx
    return x


class classproperty(property):  # noqa: N801
    """Decorator to make a property of a class."""

    def __get__(self, cls, owner=None):
        # type: (Any, Optional[type]) -> Any
        """Getter."""
        return classmethod(self.fget).__get__(None, owner)()


class SystemInfo(object):
    """System information."""

    _cpu_info = None  # type: Optional[Dict[str, str]]
    _mem_info = None  # type: Optional[Dict[str, List[str]]]

    verbose = False

    @classproperty
    def number_of_nodes(cls):  # noqa: N805
        # type: () -> int
        """Return the number of nodes."""
        info = cls._get_cpu_info()
        if "NUMA node(s)" in info:
            return int(info["NUMA node(s)"])
        else:
            return 1

    @classproperty
    def number_of_cpus(cls):  # noqa: N805
        # type: () -> int
        """Return the number of cpus."""
        info = cls._get_cpu_info()
        return int(info["CPU(s)"])

    @classproperty
    def number_of_physical_cores(cls):  # noqa: N805
        # type: () -> int
        """Return the number of physical cores."""
        info = cls._get_cpu_info()
        return int(info["Socket(s)"]) * int(info["Core(s) per socket"])

    @classproperty
    def total_memory(cls):  # noqa: N805
        # type: () -> int
        """Return the total physical memory in bytes."""
        info = cls._get_mem_info()
        return int(info["Mem"][0])

    @classmethod
    def _get_cpu_info(cls):
        # type: () -> Dict[str, str]
        if cls._cpu_info is None:
            if cls.verbose:
                sys.stderr.write("running lscpu...\n")
            info = subprocess.check_output(  # noqa: S603,S607
                ["lscpu"], env={"LANG": "C"}
            ).decode("utf-8")
            info_list = info.strip().split("\n")
            info_list_list = [[ss.strip() for ss in s.split(":")] for s in info_list]
            info_items = [(s[0], s[1]) for s in info_list_list]
            cls._cpu_info = dict(info_items)
        return cls._cpu_info

    @classmethod
    def _get_mem_info(cls):
        # type: () -> Dict[str, List[str]]
        if cls._mem_info is None:
            if cls.verbose:
                sys.stderr.write("running free...\n")
            info = subprocess.check_output(  # noqa: S603,S607
                ["free", "-b"], env={"LANG": "C"}
            ).decode("utf-8")
            info_list = info.strip().split("\n")
            info_list_list = [[ss.strip() for ss in s.split(":")] for s in info_list]
            info_pairs = [s for s in info_list_list if len(s) == 2]
            info_items = [(s[0], s[1].split()) for s in info_pairs]
            cls._mem_info = dict(info_items)
        return cls._mem_info


class Setup(object):
    """Setup parameters."""

    def __init__(self, target):
        # type: (Tuple[int, int, int]) -> None
        """Construct a set of setup parameters."""
        self._target = target  # the target version (major, minor, patch).

        # v4.2.0
        # We take "WORDSIZE32" (64-bit) values.

        self.compresssize = 90000
        self.filepatches = 256
        self.hidesize = 0
        self.largepatches = 256
        self.largesize = 50000000
        self.maxtermsize = 40000  # 64-bit
        self.numstorecaches = 4
        self.scratchsize = 50000000
        self.sizestorecache = 32768
        self.smallextension = 20000000
        self.smallsize = 10000000
        self.sortiosize = 100000
        self.termsinsmall = 100000
        self.threadbucketsize = 500
        self.threads = -1  # form
        self.threadscratchoutsize = 2500000
        self.threadscratchsize = 100000
        self.workspace = 40000000  # 64-bit

        self.bracketindexsize = 200000
        self.constindex = 128
        self.continuationlines = 15
        self.functionlevels = 30
        self.maxnumbersize = 200
        self.maxwildcards = 100
        self.parentheses = 100
        self.processbucketsize = 1000
        self.subfilepatches = 64
        self.sublargepatches = 64
        self.sublargesize = 4000000
        self.subsmallextension = 800000
        self.subsmallsize = 500000
        self.subsortiosize = 32768
        self.subtermsinsmall = 10000

        # 64-bit
        self._ptrsize = 8
        self._possize = 8
        self._wordsize = 4

        if self._target >= (4, 2, 1):
            # v4.2.1
            # We take "WITHPTHREADS" (TFORM) values.
            self.largesize = 1500000000
            self.scratchsize = 500000000
            self.smallextension = 600000000
            self.smallsize = 300000000
            self.sortiosize = 200000
            self.termsinsmall = 3000000

    def items(self):
        # type: () -> Tuple[Tuple[str, int]]
        """Return pairs of parameters and values."""
        items = [(k, v) for (k, v) in self.__dict__.items() if k[0] != "_"]
        items.sort()
        return tuple(items)  # type: ignore[return-value]

    def __str__(self):
        # type: () -> str
        """Return the string representation."""
        mem = self.calc()
        params = ["{0}: {1}".format(k, v) for (k, v) in self.items()]
        return "<Setup: {0} bytes, {1}>".format(mem, ", ".join(params))

    def copy(self):
        # type: () -> Setup
        """Return a shallow copy."""
        return copy.copy(self)

    def calc(self):
        # type: () -> int
        """Return an estimation of memory usage."""
        self.maxtermsize = max(self.maxtermsize, 200)

        self.compresssize = max(
            self.compresssize, 2 * self.maxtermsize * self._wordsize
        )
        self.sortiosize = max(self.sortiosize, self.maxtermsize * self._wordsize)

        # The strange factor WordSize**2 is used in the FORM source...
        self.scratchsize = max(
            self.scratchsize, 4 * self.maxtermsize * self._wordsize ** 2
        )
        if self.hidesize > 0:
            self.hidesize = max(
                self.hidesize, 4 * self.maxtermsize * self._wordsize ** 2
            )

        self.threadscratchsize = max(
            self.threadscratchsize, 4 * self.maxtermsize * self._wordsize ** 2
        )
        self.threadscratchoutsize = max(
            self.threadscratchoutsize, 4 * self.maxtermsize * self._wordsize ** 2
        )

        # constraints in RecalcSetups()

        self.filepatches = max(self.filepatches, self.threads)

        self.termsinsmall = round_up(self.termsinsmall, 16)

        numberofblocksinsort = 10
        minimumnumberofterms = 10
        n = numberofblocksinsort * minimumnumberofterms
        if self.threads >= 0:
            minbufsize = self.threads * (1 + n) * self.maxtermsize * self._wordsize
            if self.largesize + self.smallextension < minbufsize:
                self.largesize = minbufsize - self.smallextension

        # constraints in AllocSort()

        self.filepatches = max(self.filepatches, 4)

        self.smallsize = max(self.smallsize, 16 * self.maxtermsize * self._wordsize)

        self.smallextension = max(self.smallextension, self.smallsize * 3 // 2)

        if self.largesize > 0:
            self.largesize = max(self.largesize, 2 * self.smallsize)

        compinc = 2
        minbufsize = self.filepatches * (
            self.sortiosize + (compinc + 2 * self.maxtermsize) * self._wordsize
        )
        if self.largesize + self.smallextension < minbufsize:
            if self.largesize == 0:
                self.smallextension = minbufsize
            else:
                self.largesize = minbufsize - self.smallextension

        iotry = (
            (
                (self.largesize + self.smallextension)
                // self.filepatches
                // self._wordsize
            )
            - 2 * self.maxtermsize
            - compinc
        )  # in words
        self.sortiosize = max(self.sortiosize, iotry)  # bytes vs. words??

        # Compute the memory usage.

        mem = 0
        mem += self.scratchsize * 2 + (
            self.hidesize if self.hidesize > 0 else self.scratchsize
        )
        mem += self.workspace * self._wordsize
        mem += (self.compresssize + 10) * self._wordsize
        mem += (
            self.largesize
            + self.smallextension
            + 3 * self.termsinsmall * self._ptrsize
            + self.sortiosize
        )

        storecachesize = self._possize * 2 * self._ptrsize + self._wordsize
        # ignore the padding
        storecachesize += self.sizestorecache
        mem += storecachesize * self.numstorecaches

        if self.threads >= 1:
            mem += (
                self.threadscratchoutsize + self.threadscratchsize * 2
            ) * self.threads
            mem += self.workspace * self._wordsize * self.threads
            mem += (self.compresssize + 10) * self._wordsize * self.threads

            mem += (
                self._thread_alloc_sort(
                    self.largesize // self.threads,
                    self.smallsize // self.threads,
                    self.smallextension // self.threads,
                    self.termsinsmall,
                    self.largepatches,
                    self.filepatches // self.threads,
                    self.sortiosize,
                )
                * self.threads
            )

            mem += storecachesize * self.numstorecaches * self.threads

            sizethreadbuckets = (
                (self.threadbucketsize + 1) * self.maxtermsize + 2
            ) * self._wordsize
            if self.threadbucketsize >= 250:
                sizethreadbuckets //= 4
            elif self.threadbucketsize >= 90:
                sizethreadbuckets //= 3
            elif self.threadbucketsize >= 40:
                sizethreadbuckets //= 2
            sizethreadbuckets //= self._wordsize
            mem += (
                (
                    2 * sizethreadbuckets * self._wordsize
                    + (self.threadbucketsize + 1) * self._possize
                )
                * 2
                * self.threads
            )
            if self.threads >= 3:
                mem += (
                    self.workspace * self._wordsize // 8
                    + 2 * self.maxtermsize * self._wordsize
                ) * (self.threads - 2)

        return mem

    def _thread_alloc_sort(
        self,
        largesize,
        smallsize,
        smallextension,
        termsinsmall,
        largepatches,
        filepatches,
        sortiosize,
    ):
        # type: (int, int, int, int, int, int, int) -> int

        filepatches = max(filepatches, 4)

        smallsize = max(smallsize, 16 * self.maxtermsize * self._wordsize)

        smallextension = max(smallextension, smallsize * 3 // 2)

        if largesize > 0:
            largesize = max(largesize, 2 * smallsize)

        compinc = 2
        minbufsize = filepatches * (
            sortiosize + (compinc + 2 * self.maxtermsize) * self._wordsize
        )
        if largesize + smallextension < minbufsize:
            if largesize == 0:
                smallextension = minbufsize
            else:
                largesize = minbufsize - smallextension

        iotry = (
            ((largesize + smallextension) // filepatches // self._wordsize)
            - 2 * self.maxtermsize
            - compinc
        )  # in words
        sortiosize = max(sortiosize, iotry)  # bytes vs. words??

        return (
            largesize + smallextension + 3 * termsinsmall * self._ptrsize + sortiosize
        )

    def scale(self, total_memory, lowest_scale=0.0, human_readable=False):
        # type: (int, float, bool) -> Setup
        """
        Scale to the given memory usage.

        Search for a scaling of the given setup parameters that
        results in the requested total memory usage, and return
        the rescaled setup object. If the requested memory usage
        is too high, return the parameters with the lowest
        possible usage (scaled to lowest_scale).
        """
        sp0 = self.copy()
        # Presumably increasing MaxTermSize requires increasing WorkSpace, too.
        sp0.workspace = max(sp0.workspace, sp0.maxtermsize * 250)

        def f(x):
            # type: (float) -> Tuple[int, Setup]
            # Hopefully monotonically increasing.
            sp = sp0.copy()
            sp.smallsize = int(sp.smallsize * x)
            sp.largesize = int(sp.largesize * x)
            sp.termsinsmall = int(sp.termsinsmall * x)
            sp.scratchsize = int(sp.scratchsize * x)
            m = sp.calc()
            if human_readable:
                m = round_human_readable(m, True, False)
            return (-(total_memory - m), sp)

        miny, minsp = f(lowest_scale)
        if miny >= 0:
            return minsp
        # Optimize the memory usage by bisection.
        max_iteration = 50
        x1 = 1.0
        x2 = None  # type: Optional[float]
        y1 = f(x1)[0]
        y2 = None  # type: Optional[int]
        for _i in range(max_iteration):
            if x2 is None:
                if y1 < 0:
                    x = x1 * 2.0
                    y = f(x)[0]
                    if y > 0:
                        x2 = x
                        y2 = y
                    else:
                        x1 = x
                        y1 = y
                else:
                    x = x1 * 0.5
                    y = f(x)[0]
                    if y < 0:
                        x2 = x1
                        y2 = y1
                        x1 = x
                        y1 = y
                    else:
                        x1 = x
                        y1 = y
            else:
                x = (x1 + x2) * 0.5
                y = f(x)[0]
                if y < 0:
                    x1 = x
                    y1 = y
                else:
                    x2 = x
                    y2 = y
            if x2 is not None:
                if not (y2 is not None):
                    raise AssertionError()
                if not (x1 < x2):
                    raise AssertionError()
                if not (y1 < y2):
                    raise AssertionError()
        return f(x1)[1]


def main():
    # type: () -> None
    """Entry point."""
    # Parse the command line arguments.
    parser = argparse.ArgumentParser(
        usage=("%(prog)s [options] [--] " "[par=val].. [par+=int].. [par*=float].."),
        epilog=(
            "On non-Linux systems, the number of physical CPUs and memory "
            "available on the machine may be not automatically detected. "
            "In such a case, one cannot use the default parameters "
            "depending on those values and needs to explicitly specify "
            "--ncpus, --total-cpus and --total-memory."
        ),
        add_help=False,
    )
    parser.add_argument(
        "-h",
        "--help",
        action="store_const",
        const=True,
        help="show this help message and exit",
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        nargs="?",
        const="form.set",
        help=("output to FILE (default: no (stdout), " "FILE=form.set)"),
        metavar="FILE",
    )
    parser.add_argument(
        "-f",
        "--form",
        action="store_const",
        const=True,
        help="print tform options (e.g., -w4) and exit",
    )
    parser.add_argument(
        "-m",
        "--minos",
        action="store_const",
        const=True,
        help="print minos options (e.g., -m2x4) and exit",
    )
    parser.add_argument(
        "-u",
        "--usage",
        action="store_const",
        const=True,
        help="print expected initial memory usage and exit",
    )
    parser.add_argument(
        "-H",
        "--human-readable",
        action="store_const",
        const=True,
        help=("adjust to human-readable numbers " "(e.g., 1K, 23M, 456G)"),
    )
    parser.add_argument(
        "-1",
        "--one",
        action="store_const",
        const=-1,
        dest="ncpus",
        help="use cpus in a node on the machine (default)",
    )
    parser.add_argument(
        "--full",
        action="store_const",
        const=-99999,
        dest="ncpus",
        help="use cpus in all nodes on the machine",
    )
    parser.add_argument(
        "-n", "--ncpus", action="store", type=int, help="use N cpus", metavar="N"
    )
    parser.add_argument(
        "-p",
        "--percentage",
        action="store",
        default=75.0,
        type=float,
        help=("percentage of initial memory usage " "(default: 75.0)"),
        metavar="N",
    )
    parser.add_argument(
        "-t",
        "--target",
        action="store",
        default="4.2.1",
        type=str,
        help="target version of FORM (default: 4.2.1)",
        metavar="VER",
    )
    parser.add_argument(
        "--total-cpus",
        action="store",
        type=int,
        help="specify the total cpus on the machine",
        metavar="N",
    )
    parser.add_argument(
        "--total-memory",
        action="store",
        help="specify the total memory on the machine",
        metavar="N",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_const", const=True, help="verbose output"
    )
    parser.add_argument("args", nargs="*", help=argparse.SUPPRESS)
    args = parser.parse_args()
    pars = {}

    # NOTE: when all of `--ncpus`, `--total-cpus` and `--total-memory` are
    # specified, we don't need to access the system information.

    if args.verbose:
        SystemInfo.verbose = True

    if args.total_cpus:
        total_cpus = args.total_cpus
    else:
        total_cpus = SystemInfo.number_of_physical_cores

    if args.total_memory:
        try:
            total_memory = parse_number(args.total_memory)
        except ValueError:
            parser.error(
                "non-integer value for total memory: {0}".format(args.total_memory)
            )
    else:
        total_memory = SystemInfo.total_memory

    # Help message.
    if args.help:
        parser.print_help()
        exit(0)

    # Number of CPUs.
    if args.ncpus is not None:
        ncpus = args.ncpus
    else:
        # Use 1 node for each job by default.
        ncpus = -1
    if ncpus < 0:
        # Use (-ncpus) nodes.
        ncpus = -ncpus * (total_cpus // SystemInfo.number_of_nodes)
    ncpus = max(ncpus, 1)
    ncpus = min(ncpus, total_cpus)

    # Target version.
    target_input = args.target.split(".")
    if len(target_input) > 3 or any(not x.isdigit() for x in target_input):
        parser.error("invalid target version given: {0}".format(args.target))
    if len(target_input) == 3:
        target = (int(target_input[0]), int(target_input[1]), int(target_input[2]))
    elif len(target_input) == 2:
        target = (int(target_input[0]), int(target_input[1]), 0)
    else:
        target = (int(target_input[0]), 0, 0)

    # Setup parameter in the arguments.
    sp = Setup(target)
    sp.threads = ncpus if ncpus >= 2 else -1

    for a in args.args:
        m = re.match(r"([a-zA-Z][a-zA-Z0-9]*)([+*]?)=(.*)", a)
        if m:
            par = m.group(1).lower()
            ope = m.group(2)
            val = m.group(3)
            if par in sp.__dict__:
                # Known parameter.
                if ope == "" or ope == "+":
                    # We have par=val or par+=int.
                    try:
                        val = parse_number(val)
                    except ValueError:
                        parser.error("non-integer value for parameter: {0}".format(a))
                    if ope == "":
                        setattr(sp, par, val)
                    else:
                        setattr(sp, par, getattr(sp, par) + val)
                    continue
                else:
                    # We have par*=float.
                    try:
                        val = float(val)
                    except ValueError:
                        parser.error("non-float value for parameter: {0}".format(a))
                    setattr(sp, par, int(getattr(sp, par) * val))
                    continue
            elif ope == "":
                # Unknown parameter given by par=val. Add it to the dictionary.
                pars[par] = val
                continue
        parser.error("unrecognized argument: {0}".format(a))

    # Our resource.
    cpus = max(sp.threads, 1)
    memory = int(total_memory * args.percentage / 100.0 * cpus / total_cpus)

    # For --form option.
    if args.form:
        print("-w{0}".format(cpus))
        exit()

    # For --minos option.
    if args.minos:
        print("-m{0}x{1}".format(total_cpus // cpus, cpus))
        exit()

    sp = sp.scale(memory, human_readable=args.human_readable)

    # Final memory usage we've found.
    memory_usage = sp.calc()

    if memory_usage > memory:
        shortage = memory_usage - memory
        parser.exit(
            -1,
            ("failed to find parameters: {0} bytes shortage\n").format(
                round_human_readable(shortage, True, True)
                if args.human_readable
                else str(shortage)
            ),
        )

    # For --usage option.
    if args.usage:
        if args.human_readable:
            memory_usage_str = round_human_readable(memory_usage, True, True)
        else:
            memory_usage_str = str(memory_usage)
        print(memory_usage_str)
        exit()

    # Output.
    with open_w_or_stdout(args.output) as fi:

        def round_memory(m):
            # type: (int) -> Union[int, str]
            return round_human_readable(m, False, True) if args.human_readable else m

        print(
            (
                "# {0}{1} (cpu: {2}, mem: {3}; "
                "total cpu: {4}, total mem: {5}; {6}x{7})"
            ).format(
                parser.prog,
                (" " if len(sys.argv) >= 2 else "") + " ".join(sys.argv[1:]),
                cpus,
                round_memory(memory),
                total_cpus,
                round_memory(total_memory),
                total_cpus // cpus,
                cpus,
            ),
            file=fi,
        )

        sp0 = Setup(target)  # default value
        dic0 = dict(sp0.items())
        for k, v in sp.items():
            if k == "threads":
                # 'threads N' doesn't work, must be given by tform option -wN.
                continue
            if v == dic0[k]:
                # Don't write when same as the default value.
                continue
            if args.human_readable:
                v_str = round_human_readable(v, False, True)
            else:
                v_str = str(v)
            print("{0} {1}".format(k, v_str), file=fi)
        for k, v in pars.items():
            print("{0} {1}".format(k, v), file=fi)


if __name__ == "__main__":
    main()

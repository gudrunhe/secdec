from pySecDec.loop_integral import loop_regions
from pySecDec.code_writer import sum_package, make_package
from shutil import rmtree
from os.path import isdir
from os import chdir
import pySecDec as psd
from pySecDec.algebra import Polynomial

import sympy as sp
import numpy as np
from pprint import pprint
from pySecDec.expansion import expand_sympy

# This example is the first one loop box example in the Go Mishima paper arXiv:1812.04373

# here we define the Feynman diagram
li = psd.loop_integral.LoopIntegralFromGraph(
internal_lines = [['mt',[3,1]],['mt',[1,2]],['mt',[2,4]],['mt',[4,3]]],
external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],
powerlist=["1+n1","1+n1/2","1+n1/3","1+n1/5"],
regulators=["eps","n1"],
Feynman_parameters=["x%i" % i for i in range(1,5)], # this renames the parameters, so we get the same polynomials as in the paper

replacement_rules = [
                        # note that in those relations all momenta are incoming
                        # below we have general relations
                        ('p1*p1', 'm1sq'),
                        ('p2*p2', 'm2sq'),
                        ('p3*p3', 'm3sq'),
                        ('p4*p4', 'm4sq'),
                        ('p1*p2', 's/2-(m1sq+m2sq)/2'),
                        ('p1*p3', 't/2-(m1sq+m3sq)/2'),
                        ('p1*p4', 'u/2-(m1sq+m4sq)/2'),
                        ('p2*p3', 'u/2-(m2sq+m3sq)/2'),
                        ('p2*p4', 't/2-(m2sq+m4sq)/2'),
                        ('p3*p4', 's/2-(m3sq+m4sq)/2'),
                        ('u', '(m1sq+m2sq+m3sq+m4sq)-s-t'),
                        # these below are for our specific case
                        ('mt**2', 'mtsq'),
                        ('m1sq',0),
                        ('m2sq',0),
                        ('m3sq','mHsq'),
                        ('m4sq','mHsq'),
                        ('mHsq', 0),
                    ]
)

def f(s,t,u, mtsq):
    """
    The resulting function of the integral from arXiv:1812.04373, to order 0 in mtsq.

    Expected Result from the paper for ['s','t','u','mtsq'] = [4.0, -2.82842712475, -2.82842712475, 0.1]:
    (-1.30718609739+1.85618421207j)
    """
    v = -t/s
    lsm = np.log(complex(s/mtsq))
    return 1/(s**2*v)*(np.pi**2-2*(lsm-1j*np.pi)*(lsm+np.log(complex(v))))

def fourvec(m,p):
    # creates a four-vector from the mass and the momentum three-vector
    return np.array([np.linalg.norm([m]+p)]+p)

def get_stu(p1,p2,p3):
    def fourproduct(v1,v2):
        return v1[0]*v2[0]-np.dot(v1[1:],v2[1:])
    _s = p1+p2
    _t = p2-p3
    _u = p1-p3
    return fourproduct(_s,_s), fourproduct(_t,_t), fourproduct(_u,_u)

# here we define the fourvectors and then calculate s,t,u from them,
# to get an example numeric value with which to compare the result of the
# integration
p1 = fourvec(0,[1,0,0])
p2 = fourvec(0,[-1,0,0])
p3 = fourvec(0,[0,1,1])

s,t,u = get_stu(p1,p2,p3)
mtsq=0.1

print("When using s, t, u, mtsq = {:.5f}, {:.5f}, {:.5f}, {:.5f}".format(s,t,u,mtsq))
print("The expected result is: {:.5f}".format(f(s,t,u,mtsq)))
print()

# find the regions
generators_args = loop_regions(
    name = "box1L_expansion_by_regions",
    loop_integral=li,
    smallness_parameter = "mtsq",
    expansion_by_regions_order=0)

# write the code to sum up the regions
sum_package("box1L_expansion_by_regions", [make_package]*len(generators_args), generators_args, li.regulators,
            requested_orders = [0,0],
            real_parameters = ['s','t','u','mtsq'], complex_parameters = [])

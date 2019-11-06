# This example is the first two loop planar box example in the Go Mishima paper arXiv:1812.04373

from pySecDec.loop_integral.loop_regions import loop_regions
from shutil import rmtree
from os.path import isdir
from os import chdir
import os
import pySecDec as psd
from pySecDec.algebra import Polynomial

import sympy as sp
import numpy as np
from pprint import pprint
from pySecDec.expansion import expand_sympy

chdir(os.path.expanduser("~/loputoo/pySecDec-1.4.2/nodist_examples/box2L_expansion_by_regions"))

li = psd.loop_integral.LoopIntegralFromGraph(
internal_lines = [['mt',[1,5]],['mt',[1,2]],['mt',[2,6]],['mt',[6,4]],['mt',[4,3]],['mt',[3,5]],[0,[5,6]]],
external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],

# powerlist=["1+n1","1+2*n1","1+3*n1","1+4*n1","1+5*n1","1+6*n1","1+7*n1"],
regulators=["eps"],
Feynman_parameters=["x%i" % i for i in range(1,8)], # this renames the parameters, so we get the same polynomials as in the paper

replacement_rules = [
                        # note that in those relations all momenta are incoming, so with a negative sign for outgoing particles
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
                        # ('u', '(m1sq+m2sq+m3sq+m4sq)-s-t'),
                        # these below are for our specific case
                        ('mt**2', 'mtsq'),
                        ('m1sq',0),
                        ('m2sq',0),
                        ('m3sq',0),
                        ('m4sq',0),
                        ('mtsq', 'z*mtsq'),
                    ]
)

def f(s,t,u, mtsq):
    """
    The resulting function of the integral from 1812.04373.

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

# print("When using s, t, u, mtsq = {:.5f}, {:.5f}, {:.5f}, {:.5f}".format(s,t,u,mtsq))
# print("The expected result is: {:.5f}".format(f(s,t,u,mtsq)))
# print()

# print(f"U:\n{li.U}")
uterm = sp.sympify("(x1+x2+x3)*(x4+x5+x6)+(x1+x2+x3+x4+x5+x6)*x7")
assert (sp.sympify(li.U)-uterm).simplify()==0
print("u poly", (sp.sympify(li.U)-uterm).simplify()==0)

# print(f"F:\n{li.F}")

f_new = Polynomial.from_expression(str(li.F), ["mtsq","s","t","z"])

print(f"F2:\n{f_new}")

mtsqterm = f_new.coeffs[0]
mtsqterm2 = sp.sympify("((x1+x2+x3)*(x4+x5+x6)+(x1+x2+x3+x4+x5+x6)*x7)*(x1+x2+x3+x4+x5+x6)")
assert (mtsqterm-mtsqterm2).simplify()==0
sterm = f_new.coeffs[1]
sterm2 = -sp.sympify("x1*(x4*(x6+x7)+x3*(x4+x5+x6+x7))+x6*((x2+x3)*x4+(x3+x4)*x7)")
assert (sterm-sterm2).simplify()==0
tterm = f_new.coeffs[2]
tterm2 = -sp.sympify("x2*x5*x7")
assert (tterm-tterm2).simplify()==0
print("mtsqterm eq",(mtsqterm-mtsqterm2).simplify()==0)
print("sterm eq",(sterm-sterm2).simplify()==0)
print("tterm eq",(tterm-tterm2).simplify()==0)

if 1:
# if 0:
    if isdir("box2L_big_loop_massive_expansion_by_regions"):
        rmtree("box2L_big_loop_massive_expansion_by_regions")

    generators_args = psd.loop_integral.loop_regions(
        name = "box2L_big_loop_massive_expansion_by_regions",
        loop_integral=li,
        smallness_parameter = "z",
        expansion_by_regions_order=0,
        add_monomial_regulator_power="n1")
    for x in generators_args:
        print(x["name"])
    # generators_args = generators_args[:1]
    if __name__ == "__main__":
        psd.code_writer.sum_package("box2L_big_loop_massive_expansion_by_regions", [psd.make_package]*len(generators_args), generators_args, li.regulators,
                    requested_orders = [0,0],
                    real_parameters = ['s','t','u','mtsq','z'], complex_parameters = [])

    import configure

#!/usr/bin/env python3
# This example is the second two loop planar box example in the Go Mishima paper arXiv:1812.04373

from pySecDec.algebra import Polynomial
from pySecDec.loop_integral.loop_regions import loop_regions
import numpy as np
import pySecDec as psd
import re
import sympy as sp

li = psd.loop_integral.LoopIntegralFromGraph(
internal_lines = [['mt',[1,5]],['mt',[1,2]],['mt',[2,6]],[0,[6,4]],[0,[4,3]],[0,[3,5]],['mt',[5,6]]],
external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],

powerlist=["1+n1","1+n1/2","1+n1/3","1+n1/5","1+n1/7","1+n1/11","1+n1/13"],
regulators=["n1","eps"],
Feynman_parameters=["x%i" % i for i in range(1,8)], # this renames the parameters, so we get the same polynomials as in the paper

replacement_rules = [
                        # in the figure 4 (right), to get the F in eq (5.6), actually the left small loop has to be massless,
                        # not like shown in the figure
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
                        # ('mtsq', 'z*mtsq'),
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

# print(f"F2:\n{f_new}")

mtsqterm = sp.sympify(f_new.coeffs[0])
mtsqterm2 = sp.sympify("((x1+x2+x3)*(x4+x5+x6)+(x1+x2+x3+x4+x5+x6)*x7)*(x1+x2+x3+x7)")
assert (mtsqterm-mtsqterm2).simplify()==0
sterm = f_new.coeffs[1]
sterm2 = -sp.sympify("x1*(x4*(x6+x7)+x3*(x4+x5+x6+x7))+x6*((x2+x3)*x4+(x3+x4)*x7)")
assert (sterm-sterm2).simplify()==0
tterm = f_new.coeffs[2]
tterm2 = -sp.sympify("x2*x5*x7")
assert (tterm-tterm2).simplify()==0

mtsqterm2 = sp.sympify("((x1+x2+x3)*(x4+x5+x6)+(x1+x2+x3+x4+x5+x6)*x7)*(x1+x2+x3+x7)") * sp.sympify("mtsq")
sterm2 = -sp.sympify("x1*(x4*(x6+x7)+x3*(x4+x5+x6+x7))+x6*((x2+x3)*x4+(x3+x4)*x7)") * sp.sympify("s")
tterm2 = -sp.sympify("x2*x5*x7") * sp.sympify("t")

F = mtsqterm2+sterm2+tterm2

print((F-sp.sympify(li.F)).expand().simplify())

print("mtsqterm eq",(mtsqterm-mtsqterm2).simplify()==0)
print("sterm eq",(sterm-sterm2).simplify()==0)
print("tterm eq",(tterm-tterm2).simplify()==0)

mathematica_vectors = "{0, -1, -1, -1, 0, 0, -1}, {0, -1, -1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, -1}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 1, 1, 0, 0}, {0, 0, 1, 1, 1, 1, 1}, {0, 1, 1, 1, 0, 0, 0}"
mathematica_vectors = re.findall(r"{(.*?)}",mathematica_vectors)
mathematica_vectors = np.array( [ np.array([int(el) for el in vec.split(",")]) for vec in mathematica_vectors])
mathematica_vectors = np.array([ vec+max(0,-np.min(vec)) for vec in mathematica_vectors])
mathematica_vectors = np.array(sorted(mathematica_vectors, key=lambda x: "".join(str(y) for y in x)))
print(f"Mathematica vectors:\n{mathematica_vectors}")

paper_vectors = "(0, 0, 0, 0, 0, 0, 0), (0, 0, 1, 1, 1, 0, 0), (0, 1, 1, 1, 0, 0, 0), (1, 0, 0, 0, 1, 1, 0), (1, 1, 0, 0, 0, 1, 0), (0, 0, 1, 1, 1, 1, 1), (1, 0, 0, 1, 1, 1, 1), (1, 1, 1, 1, 1, 1, 0)"
paper_vectors = re.findall(r"\((.*?)\)",paper_vectors)
paper_vectors = np.array( [ np.array([int(el) for el in vec.split(",")]) for vec in paper_vectors])
paper_vectors = np.array([ vec+max(0,-np.min(vec)) for vec in paper_vectors])
paper_vectors = np.array(sorted(paper_vectors, key=lambda x: "".join(str(y) for y in x)))
print(f"Paper vectors:\n{paper_vectors}")

# note that one output folder for one region is not generated, because it has no 0th order term
if 1:
    generators_args = psd.loop_integral.loop_regions(
        name = "box2L_small_loop_massive_ebr",
        loop_integral=li,
        smallness_parameter = "mtsq",
        expansion_by_regions_order=0,processes=1)

    psd.code_writer.sum_package("box2L_small_loop_massive_ebr", generators_args, li.regulators,
                requested_orders = [0,0],
                real_parameters = ['s','t','u','mtsq'], complex_parameters = [])

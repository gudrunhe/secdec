# This example is the nonplanar two loop nonplanar box example corresponding to eq 6.1 in the Go Mishima paper arXiv:1812.04373

import pySecDec as psd
from pySecDec.algebra import Polynomial

import sympy as sp

li = psd.loop_integral.LoopIntegralFromGraph(
internal_lines = [[0,[2,5]],[0,[2,6]],['mt',[1,6]],['mt',[1,5]],['mt',[5,3]],['mt',[3,4]],['mt',[4,6]]],
external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],

powerlist=["1+n1","1+n1/2","1+n1/3","1+n1/5","1+n1/7","1+n1/11","1+n1/13"],
regulators=["n1","eps"],
Feynman_parameters=["x%i" % i for i in range(1,8)], # this renames the parameters, so we get the same polynomials as in the paper

replacement_rules = [
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
                    ]
)


print(f"U:\n{li.U}")
uterm = sp.sympify("(x1+x2)*(x3+x4+x5+x6+x7)+(x3+x4)*(x5+x6+x7)")
assert (sp.sympify(li.U)-uterm).simplify()==0
print("u poly", (sp.sympify(li.U)-uterm).simplify()==0)

# print(f"F:\n{li.F}")

f_new = Polynomial.from_expression(str(li.F), ["mtsq","s","t","u"])
print(f"F2:\n{f_new}")

mtsqterm = sp.sympify(f_new.coeffs[0])
mtsqterm2 = sp.sympify("(x3+x4+x5+x6+x7)*((x1+x2)*(x3+x4+x5+x6+x7)+(x3+x4)*(x5+x6+x7))")
print("mtsqterm eq",(mtsqterm-mtsqterm2).simplify()==0)
assert (mtsqterm-mtsqterm2).simplify()==0
sterm = f_new.coeffs[1]
sterm2 = -sp.sympify("(x1*x7*(x4+x5)+x2*x5*(x3+x7)+x5*x7*(x3+x4))")
print("sterm eq",(sterm-sterm2).simplify()==0)
assert (sterm-sterm2).simplify()==0
tterm = f_new.coeffs[2]
tterm2 = -sp.sympify("x1*x3*x6")
print("tterm eq",(tterm-tterm2).simplify()==0)
assert (tterm-tterm2).simplify()==0
uterm = f_new.coeffs[3]
uterm2 = -sp.sympify("x2*x4*x6")
print("uterm eq",(uterm-uterm2).simplify()==0)
assert (uterm-uterm2).simplify()==0

F = mtsqterm2* sp.sympify("mtsq")+sterm2*sp.sympify("s")+tterm2*sp.sympify("t")+uterm2*sp.sympify("u")
print("F eq", (F-sp.sympify(li.F)).expand().simplify()==0)
assert (F-sp.sympify(li.F)).simplify()==0

name = "box2L_nonplanar_expansion_by_regions_eq6p1"

generators_args = psd.loop_integral.loop_regions(
    name = name,
    loop_integral=li,
    smallness_parameter = "mtsq",
    expansion_by_regions_order=0,
    form_work_space="5000M",
    enforce_complex=True,
    processes=None)

psd.code_writer.sum_package(name, [psd.make_package]*len(generators_args), generators_args, li.regulators,
            requested_orders = [0,0],
            real_parameters = ['s','t','u','mtsq'], complex_parameters = [])

import configure

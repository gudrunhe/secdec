from __future__ import division

from argparse import ArgumentParser
from fractions import Fraction
import pySecDec as psd
import itertools as it
import re
import sympy as sp
import os
import shutil

parser = ArgumentParser()

parser.add_argument("-i", "--intfile", dest="intfile",
                    action="store", type=str, required=True,
                    help='the Reduze file to process',
                    metavar="INTFILE")

args = parser.parse_args()


kinematic_point = {}
kinematic_point["s"] = '(405/100)'
kinematic_point["t"] = '(-2025/1000)'
kinematic_point["msq"] = '1'

mandelstam_symbols = ['s','t']
mass_symbols = ['msq']
additional_prefactor = 'eps**4*(s/msq)**2'

loop_momenta = ['p1','p2']
external_momenta = ['k1','k2','k3']
replacement_rules = [
                        ('k1*k1', '0'),
                        ('k2*k2', '0'),
                        ('k3*k3', '0'),
                        ('k1*k2', 's/2'),
                        ('k1*k3', '-s/2-t/2'),
                        ('k3*k2', 't/2'),
                    ]

families = {}
families["AA_F1"] = ['(p1)**2-msq','(p2)**2-msq','(p1-p2)**2','(p1+k1)**2-msq','(p2+k1)**2-msq','(p1-k2)**2-msq','(p2-k2)**2-msq','(p1-k2-k3)**2-msq','(p2-k2-k3)**2-msq']
families["AA_F1x123"] = ['(p1)**2-msq','(p2)**2-msq','(p1-p2)**2','(p1+k2)**2-msq','(p2+k2)**2-msq','(p1-k3)**2-msq','(p2-k3)**2-msq','(p1-k3-k1)**2-msq','(p2-k3-k1)**2-msq']
families["AA_F3"] = ['(p1)**2','(p1-p2)**2-msq','(p1+k1)**2','(p2+k1)**2-msq','(p1-k2)**2','(p2-k2)**2-msq','(p2-k2-k3)**2-msq','(p1-k2-k3)**2','(p2)**2-msq']

def parse_integral(line):

    integral_pattern = re.compile(r"(\S*)\s*[0-9]*\s*[0-9]*\s*[0-9]*\s*[0-9]*\s*(.*)", re.IGNORECASE)
    dim_pattern = re.compile('dim(inc|dec)(\d+)', re.IGNORECASE)

    # Get integral name and power list
    integral_match = integral_pattern.findall(line)
    name, powerlist = integral_match[0][0], map(int, integral_match[0][1].split())

    # Parse dimension shift
    dim = 4
    dim_match =  dim_pattern.findall(name)
    name = dim_pattern.sub('',name)
    if dim_match:
        assert len(dim_match) == 1
        if dim_match[0][0] == 'inc':
            dim += int(dim_match[0][1])
        elif dim_match[0][0] == 'dec':
            dim -= int(dim_match[0][1])
        else:
            raise "Could not parse dimension shift for integral %s" % name

    print(line, name, dim, powerlist)

    return name, dim, powerlist

def parse_coefficient(line, dimension_symbol):

    # TODO - should be independent of kinematics (automate this)
    # Insert d -> 4-2*eps using Sympy and split numerator/denominator (Slow)
    d = sp.symbols(dimension_symbol)
    eps = sp.symbols("eps")
    s = sp.symbols("s")
    t = sp.symbols("t")
    coeff_n, coeff_d = sp.fraction(sp.cancel(sp.together(sp.sympify(line).subs(d,4-2*eps).subs(s,kinematic_point["s"]).subs(t,kinematic_point["t"])))) #TODO - mass
    coeff_n = psd.algebra.Expression(coeff_n,['eps']).simplify()
    coeff_d_poly = psd.algebra.Expression(coeff_d,['eps']).simplify()
    coeff_d = psd.algebra.ExponentiatedPolynomial(coeff_d_poly.expolist,coeff_d_poly.coeffs,exponent=-1,polysymbols=coeff_d_poly.polysymbols)

    # Compute Laurent expansion of coefficient
    coeff = psd.algebra.Product(coeff_n,coeff_d)
    coeff = psd.expansion.expand_singular(coeff,0,4)
 
    return coeff

integrals = []
coefficients = []
with open(args.intfile) as f:

    for line in f.readlines():

        # Count spaces, indicates the type of line to be parsed
        spaces = sum(1 for _ in it.takewhile(lambda c: c == ' ', line))
        line = line.strip()

        if spaces == 0:
            assert line.strip() == '' or line.strip() == ';'
        elif spaces == 2:
            # LHS integral
            print('LHS',line)
        elif spaces == 3:
            assert line.strip() == '1'
        elif spaces == 4:
            # RHS integral
            integral = {}
            integral["name"], integral["dim"], integral["powerlist"] = parse_integral(line)
            integrals.append(integral)

        elif spaces == 5:
            # Coefficient
            coefficients.append(parse_coefficient(line,'d'))

assert len(integrals) == len(coefficients)

# Sector decompose integrals and build library
built_integrals = []
for integral,coefficient in zip(integrals,coefficients):

    assert integral["name"] in families.keys()

    # Prepare for creation of integral directory
    name = integral["name"] + "_dim" + str(integral["dim"]) + "_" + ''.join(map(str,integral["powerlist"])).replace('-', 'm')
    shutil.rmtree(name, ignore_errors=True)

    # Determine lowest order of eps in coefficient
    coefficient_min_order = int(min(coefficient.expolist))

    li = psd.loop_integral.LoopIntegralFromPropagators(
        loop_momenta=loop_momenta,
        external_momenta=external_momenta,
        propagators=families[integral["name"]],
        replacement_rules=replacement_rules,
        powerlist=integral["powerlist"],
        dimensionality=str(integral["dim"]) + "-2*eps"
    )

    evaluated_coefficient = ''
    for coeff,exponent in zip(coefficient.coeffs,coefficient.expolist):
        numeric_coefficient = str(coeff)
        evaluated_coefficient += "+(" + numeric_coefficient + ")*eps**(" + str(exponent[0]) + ")"

    print(name)
    print(evaluated_coefficient)

    psd.loop_integral.loop_package(
        name=name,
        loop_integral=li,
        real_parameters=mandelstam_symbols+mass_symbols,
        additional_prefactor='(' + additional_prefactor + ')' + '*(' + evaluated_coefficient + ')',
        requested_order=4, # TODO - should be requested_order=0 here
        form_optimization_level=2,
        form_work_space='1G',
        decomposition_method='iterative',
        contour_deformation=True
    )

    built_integral = {}
    built_integral["name"] = name

    built_integrals.append(built_integral)

    # Build example
    for i in range(5):
        retcode = os.system('CXX=nvcc CPATH=/home/pcl335/sjones/qmc/src make -j20 -C' + name) # TODO - j should be parameter
        if retcode == 0:
            break
        else:
            print("retrying build...")

with open('results.txt','w') as resultsfile:

    for built_integral in built_integrals:
        library = psd.integral_interface.IntegralLibrary(built_integral["name"] + "/" + built_integral["name"] + "_pylink.so")
        library.use_Qmc(devices=[3],epsrel=1e-5, maxeval=10**8)
        _, _, str_integral_with_prefactor = library([float(Fraction(kinematic_point["s"].lstrip('(').rstrip(')'))), float(Fraction(kinematic_point["t"].lstrip('(').rstrip(')'))), float(Fraction(kinematic_point["msq"].lstrip('(').rstrip(')')))])

        resultsfile.write(str_integral_with_prefactor)
        print(str_integral_with_prefactor)


#!/usr/bin/env python3

from pySecDec.loop_integral import LoopIntegralFromPropagators, LoopPackage
from pySecDec.code_writer.sum_package import sum_package, Coefficient
from itertools import repeat
import re, sympy as sp, numpy as np

# arXiv:1809.06240

Mandelstam_symbols = ['v1', 'v2','v3', 'v4', 'v5', 'sqrtDelta']

def add_loop_integral(loop_package_args, powerlist):
    name = 'I73_I_%s' % '_'.join(s.replace('-','m') for s in map(str,powerlist))
    geometric_banlist = ['I73_I_1_1_1_1_1_1_0_0_0_0_0']

    loop_package_args.append(
        LoopPackage(

            name = name,

            loop_integral = LoopIntegralFromPropagators(

            loop_momenta = ['k1','k2'],
            external_momenta = ['p1','p2','p3','p4','p5'],

            propagators = ['k1**2','(k1-p1)**2','(k1-p1-p2)**2','(k1-p1-p2-p3)**2','k2**2','(k2-p1-p2-p3-p4)**2','(k1-k2)**2','(k1-k2+p4)**2','(k2-p1)**2','(k2-p1-p2)**2','(k2-p1-p2-p3)**2'],
            powerlist = powerlist,

            replacement_rules = [
                                    ('p1*p1', 0),
                                    ('p2*p2', 0),
                                    ('p3*p3', 0),
                                    ('p4*p4', 0),
                                    ('p1*p2', 'v1/2'),
                                    ('p2*p3', 'v2/2'),
                                    ('p1*p3', '(v4-v1-v2)/2'),
                                    ('p1*p4', '(v2-v5-v4)/2'),
                                    ('p2*p4', '(-v2-v3+v5)/2'),
                                    ('p3*p4', 'v3/2')
                                ]
            ),

            additional_prefactor = 'exp(2*EulerGamma*eps)',

            # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
            form_optimization_level = 4,

            # the WorkSpace parameter for FORM
            form_work_space = '2G',

            # the method to be used for the sector decomposition
            # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
            decomposition_method = 'iterative' if name in geometric_banlist else 'geometric',
            # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
            # $PATH, you can set the path to the 'normaliz' command-line
            # executable here
            #normaliz_executable='/path/to/normaliz',

            # whether or not to produce code to perform the contour deformation
            # contour deformation is not required if we only want to compute euclidean points (all Mandelstam invariants negative)
            contour_deformation = False,

            # complex return value is required for sum_package
            enforce_complex = True,

        )

    )

# the integral "I73" in terms of its master integrals
I73 = '''
           (-(v2*v3)/2 + (v2*v3*(v1*v2 - v2*v3 + v3*v4 - v1*v5 - v4*v5))/
              (2*Sqrt[Delta]))*F[-1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0] +
           ((v3*v4)/2 - (v3*v4*(v1*v2 - v2*v3 + v3*v4 + v1*v5 - v4*v5))/
              (2*Sqrt[Delta]))*F[0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0] +
           ((v4*v5)/2 - (v4*v5*(v1*v2 + v2*v3 - v3*v4 - v1*v5 + v4*v5))/
              (2*Sqrt[Delta]))*F[0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0] +
           ((-(v2*v3) - v1*v5)/2 + (v1*v2^2*v3 - v2^2*v3^2 + v2*v3^2*v4 +
               v1^2*v2*v5 - 2*v1*v2*v3*v5 - v1*v3*v4*v5 - v2*v3*v4*v5 - v1^2*v5^2 +
               v1*v4*v5^2)/(2*Sqrt[Delta]))*F[0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0] +
           ((v2*(v1 + v3))/2 - (v2*(v1^2*v2 - v2*v3^2 + v1*v3*v4 + v3^2*v4 -
                v1^2*v5 - v1*v3*v5 + v1*v4*v5 - v3*v4*v5))/(2*Sqrt[Delta]))*
            (F[0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0] + F[0, 1, 1, 1, 1, 0, 1, 1, 0, 0,
              0] + F[0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0] - F[0, 1, 1, 1, 1, 1, 1, 0, 0,
              0, 0] - F[0, 1, 1, 1, 1, 1, 1, 1, 0, 0, -1]) -
           (-(v2*v3*v4)/2 + (v2*v3*v4*(v1*v2 - v2*v3 + v3*v4 - 3*v1*v5 - v4*v5))/
              (2*Sqrt[Delta]))*F[0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0] +
           ((v3*v4)/2 - (v3*v4*(v1*v2 - v2*v3 + v3*v4 + v1*v5 - v4*v5))/
              (2*Sqrt[Delta]))*F[1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0] +
           (-(v3*v4)/2 + (v3*v4*(v1*v2 - v2*v3 + v3*v4 + v1*v5 - v4*v5))/
              (2*Sqrt[Delta]))*(F[1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0] +
             F[1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0] + F[1, 0, 1, 1, 1, 1, 0, 1, 0, 0,
              0] - F[1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0] - F[1, 0, 1, 1, 1, 1, 1, 1, 0,
              0, -1]) - ((v3*v4^2)/2 - (v3*v4^2*(v1*v2 - v2*v3 + v3*v4 + v1*v5 -
                v4*v5))/(2*Sqrt[Delta]))*F[1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0] +
           ((v4*v5)/2 - (v4*v5*(v1*v2 + v2*v3 - v3*v4 - v1*v5 + v4*v5))/
              (2*Sqrt[Delta]))*F[1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0] +
           (-(v4*v5)/2 + (v4*v5*(v1*v2 + v2*v3 - v3*v4 - v1*v5 + v4*v5))/
              (2*Sqrt[Delta]))*(F[1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0] +
             F[1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0] + F[1, 1, 0, 1, 1, 1, 0, 1, 0, 0,
              0] - F[1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0] - F[1, 1, 0, 1, 1, 1, 1, 1, 0,
              0, -1]) - ((v4^2*v5)/2 - (v4^2*v5*(v1*v2 + v2*v3 - v3*v4 - v1*v5 +
                v4*v5))/(2*Sqrt[Delta]))*F[1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0] +
           (-(v1*v5)/2 + (v1*v5*(v1*v2 - v2*v3 - v3*v4 - v1*v5 + v4*v5))/
              (2*Sqrt[Delta]))*F[1, 1, 1, -1, 1, 1, 1, 1, 0, 0, 0] +
           ((v1*(v2 + v5))/2 - (v1*(v1*v2^2 - v2^2*v3 + v2*v3*v4 - v2*v3*v5 +
                v2*v4*v5 - v3*v4*v5 - v1*v5^2 + v4*v5^2))/(2*Sqrt[Delta]))*
            (F[1, 1, 1, -1, 1, 1, 1, 1, 0, 0, 0] + F[1, 1, 1, 0, 1, 0, 1, 1, 0, 0,
              0] + F[1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0] - F[1, 1, 1, 0, 1, 1, 1, 0, 0,
              0, 0] - F[1, 1, 1, 0, 1, 1, 1, 1, 0, 0, -1]) -
           (-(v1*v4*v5)/2 + (v1*v4*v5*(v1*v2 - 3*v2*v3 - v3*v4 - v1*v5 + v4*v5))/
              (2*Sqrt[Delta]))*F[1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0] +
           (-(v1*v2)/2 + (v1*v2*(v1*v2 - v2*v3 + v3*v4 - v1*v5 + v4*v5))/
              (2*Sqrt[Delta]))*(F[1, 1, 1, -1, 1, 1, 1, 1, 0, 0, 0] +
             2*F[1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0] +
             2*F[1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0] -
             2*F[1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0] -
             2*F[1, 1, 1, 0, 1, 1, 1, 1, 0, 0, -1] + F[1, 1, 1, 1, 1, -1, 1, 1, 0, 0,
              0] + 2*F[1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0] -
             2*F[1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0] -
             2*F[1, 1, 1, 1, 1, 0, 1, 1, 0, 0, -1] + F[1, 1, 1, 1, 1, 1, -1, 1, 0, 0,
              0] - 2*F[1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0] -
             2*F[1, 1, 1, 1, 1, 1, 0, 1, 0, 0, -1] + F[1, 1, 1, 1, 1, 1, 1, -1, 0, 0,
              0] + 2*F[1, 1, 1, 1, 1, 1, 1, 0, 0, 0, -1] + F[1, 1, 1, 1, 1, 1, 1, 1,
              0, 0, -2]) + ((v1*v2*v4)/2 -
             (v1*v2*v4*(v1*v2 - v2*v3 + v3*v4 - v1*v5 - 2*v3*v5 + v4*v5))/
              (2*Sqrt[Delta]))*(-F[1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0] -
             F[1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0] - F[1, 1, 1, 1, 1, 1, 0, 1, 0, 0,
              0] + F[1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0] + F[1, 1, 1, 1, 1, 1, 1, 1, 0,
              0, -1]) - (v1*v2*v3*v4^2*v5*F[1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0])/
            Sqrt[Delta]
    '''.replace(
        'Sqrt[Delta]','sqrtDelta' # sqrtDelta = sqrt(v1^2 * (v2-v5)^2 + (v2*v3 + v4*(v5-v3))^2 + 2*v1 * (-v2^2*v3 + v4*(v3-v5)*v5 + v2*(v3*v4 + (v3 + v4) * v5)))
      ).replace(' ', '').replace('\n','')

# get the powerlists of the contributing integrals
integrals = sorted(  set( re.findall('F\[[^\]]+\]', I73) )  )
powerlists = []
for integral in integrals:
    assert integral.startswith('F[') and integral.endswith(']')
    powerlist = map(int,integral[len('F['):-len(']')].split(','))
    powerlists.append( list(powerlist) )

# get the coefficients of the integrals
I73_expr = sp.sympify(  I73.replace('[','(').replace(']',')')  ).expand()
coefficients = []
for integral in integrals:
    coefficient = I73_expr.coeff(  integral.replace('[','(').replace(']',')')  )
    numerator, denominator = sp.fraction( coefficient )
    assert (numerator / denominator - coefficient).simplify() == 0
    coefficients.append( Coefficient(['sqrtDelta*('+str(numerator)+')'],['('+str(denominator)+')*sqrtDelta'],['eps'],Mandelstam_symbols) )

# check the extracted coefficients
reconstructed_I73 = sum(
    np.prod(list(map(sp.sympify,coeff.numerators))) / np.prod(list(map(sp.sympify,coeff.denominators))) * sp.sympify('F(' + ','.join(map(str,powerlist)) + ')')
    for coeff,powerlist in zip(coefficients,powerlists)
)
assert (reconstructed_I73 - I73_expr).simplify() == 0

# prepare arguments for calls to loop_package
assert len(coefficients) == len(powerlists)
loop_package_args = []
zero_free_coefficients = []
for coeff,powerlist in zip(coefficients,powerlists):
    # only add integrals with a nonzero coefficient
    for numerator in coeff.numerators:
        if sp.sympify(numerator) == 0:
            print('I73_I_%s' % '_'.join(s.replace('-','m') for s in map(str,powerlist)), 'has zero coefficient')
            break
    else:
        zero_free_coefficients.append(coeff)
        add_loop_integral(loop_package_args, powerlist)

# call sum_package
sum_package(
    name = 'I73',
    package_generators = loop_package_args,
    regulators = ['eps'],
    requested_orders = [0],
    real_parameters = Mandelstam_symbols,
    coefficients = [zero_free_coefficients]
)

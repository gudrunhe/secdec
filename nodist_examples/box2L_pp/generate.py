#!/usr/bin/env python3

from pySecDec import sum_package, LoopPackage, Coefficient

from fractions import Fraction
import pySecDec as psd
import itertools as it
import re
import sympy as sp
import os
import shutil

def parse_integral(line, families, args):
    ''' Extracts the name, integral family, dimensions and powerlist from \
        a line that was identified as an integral. Returns a dict with \
        the arguments for loop_package '''
    integral_pattern = re.compile(r"(\S*)\s*[0-9]*\s*[0-9]*\s*[0-9]*\s*[0-9]*\s*(.*)", re.IGNORECASE)
    dim_pattern = re.compile('dim(inc|dec)(\d+)', re.IGNORECASE)

    # Get integral name and power list
    integral_match = integral_pattern.findall(line)
    name, powerlist = integral_match[0][0], map(int, integral_match[0][1].split())
    # Parse dimension shift
    dim = 4
    dim_match = dim_pattern.findall(name)

    name_family = dim_pattern.sub('',name)
    name = name + '_' + integral_match[0][1].replace(' ','_').replace('-','m')

    if dim_match:
        assert len(dim_match) == 1
        if dim_match[0][0] == 'inc':
            dim += int(dim_match[0][1])
        elif dim_match[0][0] == 'dec':
            dim -= int(dim_match[0][1])
        else:
            raise "Could not parse dimension shift for integral %s" % name

    li = psd.loop_integral.LoopIntegralFromPropagators(
            loop_momenta=args['loop_momenta'],
            external_momenta=args['external_momenta'],
            propagators=families[name_family],
            replacement_rules=args['replacement_rules'],
            powerlist=powerlist,
            dimensionality=str(dim) + "-2*eps"
        )

    #Define the package_args dict with all arguments for loop_package
    package_args = LoopPackage(
        name=name,
        loop_integral=li,
        #The requested order can not *really* be calculated yet as it depends on the order of the coefficient
        requested_order=args['requested_order'],
        real_parameters=args['mandelstam_symbols']+args['mass_symbols'],
        contour_deformation=args['contour_deformation'],
        additional_prefactor=args['additional_prefactor']
    )

    return package_args

def parse_coefficient(line, dimension_symbol, args):
    '''Decomposes the coefficent into numerator and denominator and \
    calculates the order of epsilon of the coefficient.'''
    d = sp.symbols(dimension_symbol)
    eps = sp.symbols("eps")
    coeff_n, coeff_d = sp.fraction(sp.cancel(sp.together(sp.sympify(line).subs(d,4-2*eps))))
    coefficient = Coefficient([str(coeff_n)],[str(coeff_d)],[eps],args['mandelstam_symbols']+args['mass_symbols'])
    return coefficient

def parse_reduze_file(reduze_file, families, args):

    list_package_args = []
    coefficients = []

    with open(reduze_file) as f:

        for line in f.readlines():

            # Count spaces, indicates the type of line to be parsed
            spaces = sum(1 for _ in it.takewhile(lambda c: c == ' ', line))
            line = line.strip()

            if spaces == 0:
                assert line.strip() == '' or line.strip() == ';'
            elif spaces == 2:
                #Get LHS integral name
                integral_pattern = re.compile(r"(\S*)\s*[0-9]*\s*[0-9]*\s*[0-9]*\s*[0-9]*\s*(.*)", re.IGNORECASE)
                integral_match = integral_pattern.findall(line)
                integral_name = integral_match[0][0]+ '_' + integral_match[0][1].replace(' ','_').replace('-','m')
            elif spaces == 3:
                assert line.strip() == '1'
            elif spaces == 4:
                # RHS integral
                list_package_args.append(parse_integral(line, families, args))

            elif spaces == 5:
                # Coefficient
                coefficients.append(parse_coefficient(line, 'd', args))

    assert len(list_package_args) == len(coefficients)

    # get real_parameters
    real_parameters=list_package_args[0].real_parameters
    for package_args in list_package_args:
        assert package_args.real_parameters == real_parameters

    sum_package(integral_name, list_package_args, ['eps'], [args['requested_order']],real_parameters=real_parameters,coefficients=[coefficients])

if __name__ == "__main__":

    #User Input
    name_reduze_file = 'F3_121111100.red'

    args = {}
    args['mandelstam_symbols'] = ['s','t']
    args['mass_symbols'] = ['msq']
    args['additional_prefactor'] = 'eps**4*(s/msq)**2'
    args['contour_deformation'] = True
    args['requested_order'] = 3

    args['loop_momenta'] = ['p1','p2']
    args['external_momenta'] = ['k1','k2','k3']
    args['replacement_rules'] = [
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

    parse_reduze_file(name_reduze_file, families, args)

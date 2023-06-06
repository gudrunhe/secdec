#!/usr/bin/env python3

import sympy as sp
import pySecDec as psd
from pySecDec.algebra import Polynomial
from pySecDec.loop_integral.loop_regions import  extra_regulator_constraints, suggested_extra_regulator_exponent
from pySecDec.make_regions import find_regions
from pySecDec.misc import sympify_expression

if __name__ == "__main__":

    # scalar integral
    name = 'box'
    smallness_parameter = sp.sympify('m')
    internal_lines = [['m', [1, 2]], ['m', [2, 3]], ['m', [3, 4]], ['m', [4, 1]]] 
    external_lines = [['p1', 1], ['p2', 2], ['p3', 3], ['p4', 4]]
    loop_integral = psd.LoopIntegralFromGraph(
        internal_lines = internal_lines,
        external_lines = external_lines,
        replacement_rules = [
            ('p4', '-p1-p2-p3'),
            ('p1*p1', '0'),
            ('p2*p2', '0'),
            ('p3*p3', '0'),
            ('p1*p2', 's12/2'),
            ('p1*p3', 's13/2'),
            ('p2*p3', '-s12/2-s13/2'),
        ],
    )
    # (Optional) Draw diagrams
    #psd.loop_integral.draw.plot_diagram(internal_lines,external_lines,name)

    print('-- suggested_extra_regulator_exponent --')
    regulator_exponents = suggested_extra_regulator_exponent(loop_integral, smallness_parameter)

    # pretty print
    print('  regulator_exponents: ', regulator_exponents)
    print('  I[' + ','.join([f'v{i}' for i in range(len(regulator_exponents))]) + '] -> ', 
          'I[' + ','.join([str(sympify_expression(f'v{i}+({e})*n1')) for i,e in enumerate(regulator_exponents)]) + ']')
    print()

    print('-- extra_regulator_constraints --')
    # Construct Lee-Pomeransky (U+F) polytope with variables [x0,...,xn,m]
    poly = Polynomial.from_expression( str(loop_integral.F + loop_integral.U), loop_integral.preliminary_F.symbols + [smallness_parameter])
    idx = poly.symbols.index(smallness_parameter)
    all_constraints = extra_regulator_constraints(idx, poly, 0, loop_integral.regulators, loop_integral.powerlist, loop_integral.dimensionality)

    # pretty print
    for region, constraints in all_constraints.items():
        print('  region: ', region)
        print('    constraints: ', constraints.tolist())
        for constraint in constraints:
            expr = '+'.join([ f'({e})*v{i}' for i,e in enumerate(constraint) if e != 0])
            print(f'     {constraint} => ({sympify_expression(expr)}) != 0')
    print()

    print('-- find_regions --')
    print('  regions: ', find_regions(idx, poly).tolist())

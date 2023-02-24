#!/usr/bin/env python3
import pySecDec as psd
from pySecDec import sum_package, loop_regions
import os
import json
import sympy as sp

if __name__ == "__main__":
    li = psd.LoopIntegralFromGraph(
        internal_lines = [['mt',[1,4]],['mw',[4,2]],[0,[2,3]],[0,[4,5]], [0, [1,5]], ['mz', [5,3]]],
        external_lines = [['p1',1],['p2',2],['p3',3]],
        regulators=['eps'],
        replacement_rules = [
            ('p1*p1', 's'),
            ('p2*p2', 0),
            ('p3*p3', 0),
            ('p1*p2', 's/2'),
            ('p1*p3', 's/2'),
            ('p2*p3', 's/2'),
            ('mw**2', 'mwsq'),
            ('mz**2', 'mzsq'),
            ('mt**2', 'mtsq')
        ]
    )

    name = 'muon_decay2L_mass_expansion'
    smallness_parameter = 'mtsq'

    terms = loop_regions(
        name = name,
        loop_integral=li,
        smallness_parameter = smallness_parameter,
        decomposition_method = 'geometric',
        form_optimization_level = 2,
        expansion_by_regions_order=1)

    term_by_prefactor_exponent = {}
    for term in terms:
        prefactor_coefficient, prefactor_smallvar_exponent = sp.sympify(str(term.prefactor)).as_coeff_exponent(sp.sympify(smallness_parameter))
        term = term._replace(prefactor = prefactor_coefficient) #remove mass prefactor from full prefactor
        term_by_prefactor_exponent.setdefault(str(prefactor_smallvar_exponent), [])
        term_by_prefactor_exponent[str(prefactor_smallvar_exponent)].append(term)

    if not os.path.exists(name):
        os.mkdir(name)
    os.chdir(name)

    prefactor_exponent_by_name = {}
    for prefactor_number, prefactor_exponent in enumerate(sorted(term_by_prefactor_exponent.keys())):
        prefactor_number += 1
        prefactor_exponent_by_name['prefactor_' + str(prefactor_number)] = prefactor_exponent
        sum_package(
            'prefactor_' + str(prefactor_number),
            term_by_prefactor_exponent[prefactor_exponent],
            regulators = ['eps'],
            requested_orders = [0],
            real_parameters = ['s', 'mwsq', 'mzsq'],
            complex_parameters = [],
            processes = 30)
    with open('prefactor_exponent_by_name.json', 'w') as f:
        json.dump(prefactor_exponent_by_name, f)
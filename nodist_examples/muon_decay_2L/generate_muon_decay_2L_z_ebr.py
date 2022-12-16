#!/usr/bin/env python3
import pySecDec as psd
from pySecDec import sum_package, loop_regions

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
            ('mt**2', 'z*mtsq')
        ]
    )

    # find the regions
    generators_args = loop_regions(
        name = "twoloop_EBR_z",
        loop_integral=li,
        smallness_parameter = 'z',
        decomposition_method = 'geometric',
        form_optimization_level = 4,
        expansion_by_regions_order=1)

    zeroth_order_contributions = []
    first_order_contributions = []

    # Separate the zeroth and first order contributions to integrate separately later
    for packages in generators_args:
        print(packages.name[-1])
        order = packages.name[-1]
        print(packages)
        if order == '0':
            print('appending zero order')
            zeroth_order_contributions.append(packages)
        if order == '1':
            print('appending first order')
            first_order_contributions.append(packages)

    # write the code to sum up the regions
    if zeroth_order_contributions:
        sum_package("twoloop_EBR_zeroth_order_z",
                zeroth_order_contributions,
                regulators = ['eps'],
                requested_orders = [0],
                real_parameters = ['z', 's', 'mwsq', 'mzsq', 'mtsq'],
                complex_parameters = [],
                processes = 30)

    if first_order_contributions:
        sum_package("twoloop_EBR_first_order_z",
                first_order_contributions,
                regulators = ['eps'],
                requested_orders = [0],
                real_parameters = ['z', 's', 'mwsq', 'mzsq', 'mtsq'],
                complex_parameters = [],
                processes = 30)

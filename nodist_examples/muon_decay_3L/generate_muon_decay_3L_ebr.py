#!/usr/bin/env python3
import pySecDec as psd
from pySecDec import sum_package, loop_regions

if __name__ == "__main__":

    li = psd.LoopIntegralFromGraph(
        internal_lines = [['mt',[1,4]],['mw',[4,2]],[0,[2,3]],[0,[4,5]], [0, [1,6]], [0, [6,7]], [0, [5,7]], ['mz', [7,3]],['mz', [6,5]]],
        external_lines = [['p1',1],['p2',2],['p3',3]],
        replacement_rules = [
            ('p1*p1', 'mwsq'),
            ('p2*p2', 0),
            ('p3*p3', 0),
            ('p1*p2', 'mwsq/2'),
            ('p1*p3', 'mwsq/2'),
            ('p2*p3', 'mwsq/2'),
            ('mw**2', 'mwsq'),
            ('mz**2', 'mzsq'),
            ('mt**2', 'mtsq')
        ]
    )

    # find the regions
    generators_args = loop_regions(
        name = name,
        loop_integral=li,
        smallness_parameter = "mtsq",
        decomposition_method = 'geometric',
        form_optimization_level = 4,
        expansion_by_regions_order=1)

    zeroth_order_contributions = []
    first_order_contributions = []

    # Separate the zeroth and first order contributions to integrate separately later
    for packages in generators_args:
        order = packages.name[-1]
        if order == '0':
            zeroth_order_contributions.append(packages)
        if order == '1':
            first_order_contributions.append(packages)

    if zeroth_order_contributions:
        sum_package("threeloop_EBR_zeroth_order",
                zeroth_order_contributions,
                li.regulators,
                requested_orders = [0],
                real_parameters = ['mwsq', 'mzsq', 'mtsq'],
                complex_parameters = [],
                processes = 30)
    if first_order_contributions:
        sum_package("threeloop_EBR_first_order",
                first_order_contributions,
                li.regulators,
                requested_orders = [0],
                real_parameters = ['mwsq', 'mzsq', 'mtsq'],
                complex_parameters = [],
                processes = 30)

#!/usr/bin/env python3

from pySecDec import sum_package, make_package, make_regions

if __name__ == "__main__":

    # Expected Result:
    # -Log[t] + O[t]

    regions_generators = make_regions(
        name = 'make_regions_ebr',
        integration_variables = ['x'],
        regulators = ['delta'],
        requested_orders = [0],
        smallness_parameter = 't',
        polynomials_to_decompose = ['(x)**(delta)','(t + x + x**2)**(-1)'],
        expansion_by_regions_order = 0,
        real_parameters = ['t'],
        complex_parameters = [],
        decomposition_method = 'geometric_infinity_no_primary',
        polytope_from_sum_of=[1]
    )

    sum_package(
        'make_regions_ebr',
        regions_generators,
        regulators = ['delta'],
        requested_orders = [0],
        real_parameters = ['t']
    )

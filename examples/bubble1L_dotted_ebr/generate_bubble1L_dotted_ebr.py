#!/usr/bin/env python3

# in this example we expand the dotted loop using the mass of the particle
# in the loop as the small parameter
# here we do it both ways simultaneously: expanding in the msq directly, 
# or replacing it by z*msq, expanding in z and then using z=1

from pySecDec.loop_integral import LoopIntegralFromPropagators
from pySecDec.loop_integral import loop_regions
from pySecDec.code_writer import sum_package

if __name__ == "__main__":

    # define the loop integral for the case where we expand in m
    li_m = LoopIntegralFromPropagators(
        propagators=("((k+p)**2)", "(k**2-msq_)"),
        loop_momenta=["k"],
        powerlist=[1,2],
        regulators=['eps'],
        replacement_rules = [('p*p', 'psq'), ('msq_', 'msq')]
    )

    # define the loop integeral for the case where we expand in z
    li_z = LoopIntegralFromPropagators(
        propagators=li_m.propagators,
        loop_momenta=li_m.loop_momenta,
        regulators=li_m.regulators,
        powerlist=li_m.powerlist,
        replacement_rules = [('p*p', 'psq'), ('msq_', 'msq*z')]
    )

    for name, real_parameters, smallness_parameter, li in (
        ("bubble1L_dotted_z", ['psq','msq','z'], "z", li_z),
        ("bubble1L_dotted_m", ['psq','msq'], "msq", li_m)
    ):
        # find the regions and expand the integrals using expansion by regions
        generators_args = loop_regions(
                name = name,
                loop_integral=li,
                smallness_parameter = smallness_parameter,
                expansion_by_regions_order=1)

        # generate code that will calculate the sum of all regions and all orders in
        # the smallness parameter
        sum_package(name, generators_args, li.regulators,
                        requested_orders = [0],
                        real_parameters = real_parameters,
                        complex_parameters = [])

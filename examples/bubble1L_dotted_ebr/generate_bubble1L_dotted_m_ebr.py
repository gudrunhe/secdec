#!/usr/bin/env python3

# In this example we expand the 1loop bubble with one massive dotted propagator,
# using the mass of the propagator as the small parameter.
# Here we demonstrate the method where we expand directly in msq,

from pySecDec.loop_integral import LoopIntegralFromPropagators
from pySecDec.loop_integral import loop_regions
from pySecDec.code_writer import sum_package

if __name__ == "__main__":

    # define the loop integral for the case where we expand in msq
    li_m = LoopIntegralFromPropagators(
        propagators=("((k+p)**2)", "(k**2-msq_)"),
        loop_momenta=["k"],
        powerlist=[1,2],
        regulators=['eps'],
        replacement_rules = [('p*p', 'psq'), ('msq_', 'msq')]
    )

    # find the regions and expand the integrals using expansion by regions
    regions_generator_args = loop_regions(
        name = "bubble1L_dotted_m",
        loop_integral=li_m,
        smallness_parameter = "msq",
        expansion_by_regions_order=1)

    # generate code that will calculate the sum of all regions and the requested
    # orders in the smallness parameter
    sum_package("bubble1L_dotted_m", regions_generator_args, li_m.regulators,
                        requested_orders = [0],
                        real_parameters = ['psq','msq'],
                        complex_parameters = [])

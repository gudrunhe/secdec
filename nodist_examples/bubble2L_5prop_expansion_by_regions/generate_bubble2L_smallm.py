#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd

# example is 2-loop bubble with five propagators, two of them massive


li = psd.loop_integral.LoopIntegralFromPropagators(
propagators = ['k1**2-msq','(k1+k2)**2-msq','(k1+p1)**2','k2**2','(k1+k2+p1)**2'],
loop_momenta = ['k1','k2'],
external_momenta = ['p1'],
replacement_rules = [ ('p1*p1', 'psq'), ('msq', 'z*msq') ]
)

# find the regions and expand the integrals using expansion by regions
regions_generator_args = psd.loop_integral.loop_regions(
    name = 'bubble2L_smallm',
    loop_integral = li,
    smallness_parameter = 'z',
    expansion_by_regions_order = 2,
)

# generate code that will calculate the sum of all regions and all orders in
# the smallness parameter
psd.code_writer.sum_package('bubble2L_smallm',
    [psd.make_package]*len(regions_generator_args),
    regions_generator_args, regulators = ['eps'],requested_orders = [0],
    real_parameters = ['psq','msq','z'],
    complex_parameters = [],)

# the following Python script will set the integrator, see the file for more
import configure

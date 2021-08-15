#! /usr/bin/env python
from pySecDec.loop_integral import loop_regions
from pySecDec.code_writer import sum_package, make_package
from pySecDec.algebra import Polynomial
import pySecDec as psd

# define triangle2L (case 8)
li = psd.loop_integral.LoopIntegralFromPropagators(
propagators = ['(k1+p1)**2','(k1-p2)**2','(k2+p1)**2','(k2-p2)**2','k2**2','(k1-k2)**2-msq'],
loop_momenta = ['k1','k2'],
external_momenta = ['p1','p2'],
replacement_rules = [('p1*p1', 0), ('p2*p2', 0), ('p1*p2', 's/2') ])

# define real parameters
Mandelstam_symbols = ['s']
mass_symbols = ['msq']
realp = Mandelstam_symbols + mass_symbols

# find the regions and expand the integrals using expansion by regions
regions_generator_args = loop_regions(
    name = "triangle2L_ebr",
    loop_integral = li,
    smallness_parameter = 's',
    decomposition_method = 'geometric',
    additional_prefactor='1/gamma(2+2*eps)',    
    expansion_by_regions_order = -1)

# call to sum_package
sum_package("triangle2L_ebr",
    regions_generator_args,
    regulators = ['eps'],
    requested_orders = [0],
    real_parameters = realp,
    complex_parameters = [])

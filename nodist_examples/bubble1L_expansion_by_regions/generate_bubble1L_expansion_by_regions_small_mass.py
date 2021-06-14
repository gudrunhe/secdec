from pySecDec.code_writer import sum_package, make_package
from pySecDec.loop_integral import loop_regions
import pySecDec as psd

# Expected Result for s=1, msq=0.01
# + ((0,0) +/- (0,0))*eps^-2 + ((1,-1.38219e-09) +/- (1.93775e-07,2.76537e-07))*eps^-1 + ((1.53572,3.07812) +/- (7.28103e-07,9.42614e-07)) + O(eps)

# here we define the Feynman graph we are using
li = psd.loop_integral.LoopIntegralFromGraph(
        internal_lines = [['m',[1,2]],['m',[2,1]]],
        external_lines = [['p',1],['p',2]],
        replacement_rules = [('p*p', 's'), ('m*m', 'z*msq')])

# find the regions and expand the integrals using expansion by regions
regions_generator_args = loop_regions(
    name = 'bubble1L_expansion_by_regions_small_mass',
    loop_integral = li,
    smallness_parameter = 'z',
    expansion_by_regions_order = 2)

# generate code that will calculate the sum of all regions and all orders in
# the smallness parameter
sum_package('bubble1L_expansion_by_regions_small_mass',
    [make_package]*len(regions_generator_args),
    regions_generator_args, regulators = ['eps'],requested_orders = [0],
    real_parameters = ['s','msq','z'],
    complex_parameters = [])

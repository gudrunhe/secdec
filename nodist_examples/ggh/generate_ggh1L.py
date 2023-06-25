import pySecDec as psd
from pySecDec.code_writer.sum_package import sum_package, Coefficient
from pySecDec.loop_integral import loop_package

from parser_tools import *

if __name__ == '__main__':
    config_dict = parse_config('config_1loop')
    integrals, coefficients = parse_amplitude('amp_ggh1L.sum',config_dict)

    # Undo simplifications made in coefficients file
    coefficients = [[x.replace('d','(4-2*eps)').replace('x','(s12/mt2)') for x in coefficients[0]]]

    sum_package('ggh1L', integrals, coefficients = coefficients, 
                real_parameters=config_dict['real_parameters'], 
                regulators=['eps'], requested_orders=[2]
               )

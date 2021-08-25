#! /usr/bin/env python
#import pySecDec as psd
from pySecDec.loop_integral import LoopPackage
from pySecDec.code_writer import sum_package
#from pySecDec.code_writer.sum_package import Coefficient

# g g > g h  amplitude M++- 
# M++- = Prefactor*( P1*li1 + P2*li2 + ... + P14*li14 ) 

# import integrals and coefficients
import integralsF212
import coefficientsF212

# define input to loop_package
# args = []
mynames = ['g1_mtsq', 'g2_s12_mtsq', 'g2_s13_mtsq', 'g3_s23_mtsq','g4_mhsq_mtsq','g5_s12_mtsq','g5_s13_mtsq','g6_s23_mtsq','g7_s12_mhsq_mtsq',
                                   'g7_s13_mhsq_mtsq','g8_s23_mhsq_mtsq','g9_s12_s23_mhsq_mtsq','g9_s12_s13_mhsq_mtsq','g9_s23_s13_mhsq_mtsq']

if __name__ == '__main__':

    # define input to loop_package
    package_generators = []
    names = mynames
    for name, integral in zip(names, integralsF212.I):
        package_generators.append(
            LoopPackage(
                name = name,
                loop_integral = integral,
                requested_orders = [0]
                )
        )

    # generate code sum of (int * coeff)
    sum_package(
        'F212',
        package_generators,
        regulators = ['eps'],
        requested_orders = [0],
        real_parameters = ['s12','s13','s23','mtsq','mhsq'],
        coefficients = coefficientsF212.coeff,
        complex_parameters = []
    )

# import configure 
# import configure

# , complex_parameters = []
























    

#! /usr/bin/env python
import pySecDec as psd
from pySecDec.code_writer.sum_package import Coefficient

# 4-photon amplitude M++-- 
# Amp4y_ppmm = -8*( 1 + (t**2+u**2)/s *Box6dim(t,u) +  (t-u)/s*(  BubbleD(u)-BubbleD(t) ) )

# import integrals and coefficients
import integrals
import coefficients 

# define input to loop_package
args = []
names = ['bubble_u', 'bubble_t', 'box_6', 'box_8']
for i in range(len(integrals.I)):
    args.append({
        'name'              : names[i],
        'loop_integral'     : integrals.I[i],
        'requested_orders'  : [0]
        })

# generate code that will calculate the sum of all integrals times the corresponding coefficients 
psd.code_writer.sum_package('yyyy1L_amp',
    [psd.loop_integral.loop_package]*len(args), 
    args, regulators = ['eps'], requested_orders = [0],
    real_parameters = ['t','u'],
    coefficients = coefficients.coeff,
    complex_parameters = [])

# import configure 
import configure


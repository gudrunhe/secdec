#! /usr/bin/env python
import pySecDec as psd
from pySecDec.code_writer.sum_package import Coefficient

# 4-photon amplitude M++-- 
# Amp4y_ppmm = -8*( 1 + (t**2+u**2)/s *Box6dim(t,u) +  (t-u)/s*(  BubbleD(u)-BubbleD(t) ) )

# import integrals and coefficients
import integrals
import coefficients 

# define input to loop_package for one loop bubble (u)
args_bu = {
    'name' : 'bu',
    'loop_integral' : integrals.li_bu,
    'requested_orders' : [0]
    }

# define input to loop_package for one loop bubble (t)
args_bt = {
    'name' : 'bt',
    'loop_integral' : integrals.li_bt,
    'requested_orders' : [0]
    }

# define input to loop_package for one loop box
args_box = {
    'name' : 'box',
    'loop_integral' : integrals.li_box,
    'requested_orders' : [0]
    }

# collect args_integral in args as list
args = [args_bu, args_bt, args_box]

# generate code that will calculate the sum of all integrals times the corresponding coefficients 
psd.code_writer.sum_package('yyyy1L_amp',
    [psd.loop_integral.loop_package]*len(args),
    args, regulators = ['eps'], requested_orders = [0],
    real_parameters = ['t','u'],
    coefficients = coefficients.coeff,
    complex_parameters = [])

# import configure 
import configure


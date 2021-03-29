#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

	#
	# Example provided by Go Mishima 16.10.2019
	# - Demonstrates VEGAS with a low nstart giving wrong results (with an underestimated error)
	#

	li = psd.loop_integral.LoopIntegralFromPropagators(
	propagators =[
			'-m2+(l1+q1)**2',
			'-m2+l1**2',
			'-m2+(l1-q2)**2',
			'-m2+(l2-q2)**2',
			'-m2+(l2+q1+q3)**2',
			'-m2+(l2+q1)**2',
			'(l1-l2)**2',
		],
	powerlist = [1,1,1,1,1,1,1],
	loop_momenta = ['l1','l2'],
	replacement_rules = [
			('q1*q1',0),
			('q2*q2',0),
			('q3*q3',0),
			('q1*q2','s/2'),
			('q1*q3','t/2'),
			('q2*q3','-s/2-t/2')
		]
	)
	Mandelstam_symbols = ['s','t']
	mass_symbols = ['m2']
	loop_package(
	name = 'box2L',
	loop_integral = li,
	real_parameters = Mandelstam_symbols + mass_symbols,
	additional_prefactor = 'm2*s**2',
	requested_order = 0,
	# the optimization level to use in FORM (can be 0, 1, 2, 3)
	form_optimization_level = 2,
	decomposition_method = 'iterative',
	contour_deformation = True
	)

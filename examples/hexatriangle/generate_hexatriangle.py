#!/usr/bin/env python3
import os
import subprocess
import tempfile
import multiprocessing
import pySecDec as psd

if __name__ == '__main__':
    psd.loop_package(
        name = 'hexatriangle',
        loop_integral = psd.loop_integral.LoopIntegralFromPropagators(
            loop_momenta = ['l1','l2'],
            external_momenta = ['q1','q2','p1','p2'],
            regulator = 'eps',
            propagators = [
                '(l1)^2',
                '(l1 - p1 - p2)^2',
                '(l1 - q1)^2',
                '(l1 - q1 - q2)^2',
                '(l2)^2-mt2',
                '(l1 - p1)^2-mt2',
                '(l1 + l2 - p1 - p2)^2-mt2',
                '(l1 + l2 - q1 - q2)^2-mt2',
            ],
            powerlist = [1,1,1,1,2,1,2,1],
            dimensionality = '8-2*eps',
            replacement_rules = [
                ('p1*p1', 'mt2'),
                ('p1*p2', 'mh2/2 - mt2 + x12/2 - x35/2 - x54/2'),
                ('p1*q1', 'x12/2 + x23/2 - x54/2'),
                ('p1*q2', '-x23/2'),
                ('p2*p2', 'mt2'),
                ('p2*q1', '-x41/2'),
                ('p2*q2', 'x12/2 - x35/2 + x41/2'),
                ('q1*q1', '0'),
                ('q1*q2', 'x12/2'),
                ('q2*q2', '0')
            ]
        ),
        additional_prefactor = '-(2^(-8 + 4*eps)*pi^(-4 + 2*eps))',
        requested_orders = [0],
        contour_deformation = True,
        real_parameters = ['mh2', 'mt2', 'x12', 'x23', 'x35', 'x41', 'x54'],
        decomposition_method = 'geometric_ku',
        form_optimization_level = 4
    )

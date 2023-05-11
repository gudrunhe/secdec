#!/usr/bin/env python3
from pySecDec.integral_interface import DistevalLibrary

library = DistevalLibrary("pentabox_offshell/disteval/pentabox_offshell.json")
result = library(parameters=dict(mm=0.5, s01=2.2, s02=2.3, s03=2.4, s12=2.5, s13=2.6, s23=2.7), epsrel=1e-2)
print(result)

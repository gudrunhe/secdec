from __future__ import division

from argparse import ArgumentParser
from fractions import Fraction
import pySecDec as psd
import itertools as it
import re
import sympy as sp
import os
import shutil

kinematic_point = {}
kinematic_point["s"] = '405/100'
kinematic_point["t"] = '-2025/1000'
kinematic_point["msq"] = '1'

built_integrals = []
built_integrals.append({"name": "AA_F1_dim4_001011000"})
built_integrals.append({"name": "AA_F1_dim6_021101000"})
built_integrals.append({"name": "AA_F1_dim8_023101000"})
built_integrals.append({"name": "AA_F3_dim4_001111000"})
built_integrals.append({"name": "AA_F3_dim6_011210100"})
built_integrals.append({"name": "AA_F3_dim8_111111100"})
built_integrals.append({"name": "AA_F3_dim8_211111100"})
built_integrals.append({"name": "AA_F1_dim4_010101000"})
built_integrals.append({"name": "AA_F1_dim6_021101010"})
built_integrals.append({"name": "AA_F1_dim8_031101010"})
built_integrals.append({"name": "AA_F3_dim4_011010000"})
built_integrals.append({"name": "AA_F3_dim6_021010200"})
built_integrals.append({"name": "AA_F3_dim8_111111200"})
built_integrals.append({"name": "AA_F3_dim8_220110100"})
built_integrals.append({"name": "AA_F1_dim4_110000000"})
built_integrals.append({"name": "AA_F1_dim6_022101000"})
built_integrals.append({"name": "AA_F1x123_dim4_001011000"})
built_integrals.append({"name": "AA_F3_dim4_011010100"})
built_integrals.append({"name": "AA_F3_dim6_021111000"})
built_integrals.append({"name": "AA_F3_dim8_121111100"})
built_integrals.append({"name": "AA_F1_dim4_110101000"})
built_integrals.append({"name": "AA_F1_dim6_102012000"})
built_integrals.append({"name": "AA_F1x123_dim6_002022000"})
built_integrals.append({"name": "AA_F3_dim6_001111100"})
built_integrals.append({"name": "AA_F3_dim6_111010100"})
built_integrals.append({"name": "AA_F3_dim8_121210100"})
built_integrals.append({"name": "AA_F1_dim6_002022000"})
built_integrals.append({"name": "AA_F1_dim8_022101010"})
built_integrals.append({"name": "AA_F1x123_dim6_102012000"})
built_integrals.append({"name": "AA_F3_dim6_011111100"})
built_integrals.append({"name": "AA_F3_dim6_120110100"})
built_integrals.append({"name": "AA_F3_dim8_131010200"})

with open('results.txt','w') as resultsfile:

    for built_integral in built_integrals:
        library = psd.integral_interface.IntegralLibrary(built_integral["name"] + "/" + built_integral["name"] + "_pylink.so")
        library.use_Qmc(devices=[3],epsrel=1e-7, maxeval=10**9)
        _, _, str_integral_with_prefactor = library([float(Fraction(kinematic_point["s"])), float(Fraction(kinematic_point["t"])), float(Fraction(kinematic_point["msq"]))])

        resultsfile.write(str_integral_with_prefactor)
        print(str_integral_with_prefactor)

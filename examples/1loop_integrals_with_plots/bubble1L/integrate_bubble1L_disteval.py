from pySecDec.integral_interface import DistevalLibrary

# load each integrator
bubble1L_m0 = DistevalLibrary('bubble1L_m0/disteval/bubble1L_m0.json', verbose=False)
bubble1L_m1 = DistevalLibrary('bubble1L_m1/disteval/bubble1L_m1.json', verbose=False)
bubble1L_m2 = DistevalLibrary('bubble1L_m2/disteval/bubble1L_m2.json', verbose=False)

#integrate each configuration
result = bubble1L_m0(parameters={"s": 5.0, "m1sq": 1.0, "m2sq": 1.0, "m3sq": 1.0},
               epsrel=1e-6, epsabs=1e-12, format="json", verbose=True)
values_m0 = result["sums"]["bubble1L_m0"]
result = bubble1L_m1(parameters={"s": 5.0, "m1sq": 1.0, "m2sq": 1.0, "m3sq": 1.0},
               epsrel=1e-6, epsabs=1e-12, format="json", verbose=True)
values_m1 = result["sums"]["bubble1L_m1"]
result = bubble1L_m2(parameters={"s": 5.0, "m1sq": 1.0, "m2sq": 1.0, "m3sq": 1.0},
               epsrel=1e-6, epsabs=1e-12, format="json", verbose=True)
values_m2 = result["sums"]["bubble1L_m2"]

# Print and comapare to analytic results
# Analytic results at s=5, msq=1
print('\n2-mass:')
print('Numerical Result:')
print('eps^-1', values_m2[(-1,)][0], '+/- (', values_m2[(-1,)][1], ')')
print('eps^0', values_m2[(0,)][0], '+/- (', values_m2[(0,)][1], ')')
print('Analytic Result:')
print('eps^-1 : 1')
print('eps^0 : 0.99237539413446318 + 1.4049629462081437j')
print("-------------------------")

print('1-mass:')
print('Numerical Result:')
print('eps^-1', values_m1[(-1,)][0], '+/- (', values_m1[(-1,)][1], ')')
print('eps^0', values_m1[(0,)][0], '+/- (', values_m1[(0,)][1], ')')
print('Analytic Result:')
print('eps^-1 : 1')
print('eps^0 : 0.31374884620255472 + 2.5132741228718323j')
print("-------------------------")

print('0-mass:')
print('Numerical Result:')
print('eps^-1', values_m0[(-1,)][0], '+/- (', values_m0[(-1,)][1], ')')
print('eps^0', values_m0[(0,)][0], '+/- (', values_m0[(0,)][1], ')')
print('Analytic Result:')
print('eps^-1 : 1')
print('eps^0 : -0.18665357733563315 + 3.1415926535897900j')
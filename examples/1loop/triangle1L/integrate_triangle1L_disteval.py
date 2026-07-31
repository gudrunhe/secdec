from pySecDec.integral_interface import DistevalLibrary

# load each integrator
triangle1L_m0 = DistevalLibrary('triangle1L_m0/disteval/triangle1L_m0.json', verbose=False)
triangle1L_m1 = DistevalLibrary('triangle1L_m1/disteval/triangle1L_m1.json', verbose=False)
triangle1L_m2 = DistevalLibrary('triangle1L_m2/disteval/triangle1L_m2.json', verbose=False)
triangle1L_m3 = DistevalLibrary('triangle1L_m3/disteval/triangle1L_m3.json', verbose=False)

#integrate each configuration
result = triangle1L_m0(parameters={"s": 5.0, "m1sq": 1.0, "m2sq": 1.0, "m3sq": 1.0},
               epsrel=1e-6, epsabs=1e-12, format="json", verbose=True)
values_m0 = result["sums"]["triangle1L_m0"]
result = triangle1L_m1(parameters={"s": 5.0, "m1sq": 1.0, "m2sq": 1.0, "m3sq": 1.0},
               epsrel=1e-6, epsabs=1e-12, format="json", verbose=True)
values_m1 = result["sums"]["triangle1L_m1"]
result = triangle1L_m2(parameters={"s": 5.0, "m1sq": 1.0, "m2sq": 1.0, "m3sq": 1.0},
               epsrel=1e-6, epsabs=1e-12, format="json", verbose=True)
values_m2 = result["sums"]["triangle1L_m2"]
result = triangle1L_m3(parameters={"s": 5.0, "m1sq": 1.0, "m2sq": 1.0, "m3sq": 1.0},
               epsrel=1e-6, epsabs=1e-12, format="json", verbose=True)
values_m3 = result["sums"]["triangle1L_m3"]

# Print and comapare to analytic results
# Analytic results at s=5, msq=1
print('\n3-mass:')
print('Numerical Result:')
print('eps^0', values_m3[(0,)][0], '+/- (', values_m3[(0,)][1], ')')
print('Analytic Result:')
print('eps^0 : -0.89433451187805812 - 0.60470861377111462j')
print("-------------------------")

print('2-mass:')
print('Numerical Result:')
print('eps^0', values_m2[(0,)][0], '+/- (', values_m2[(0,)][1], ')')
print('Analytic Result:')
print('eps^0 : -0.35674383225332618 - 1.0112396644223722j')
print("-------------------------")

print('1-mass:')
print('Numerical Result:')
print('eps^0', values_m1[(0,)][0], '+/- (', values_m1[(0,)][1], ')')
print('Analytic Result:')
print('eps^0 : 0.0268892987355248345 - 1.1257956771053608j')
print("-------------------------")

print('0-mass:')
print('Numerical Result:')
print('eps^-2', values_m0[(-2,)][0], '+/- (', values_m0[(-2,)][1], ')')
print('eps^-1', values_m0[(-1,)][0], '+/- (', values_m0[(-1,)][1], ')')
print('eps^0', values_m0[(0,)][0], '+/- (', values_m0[(0,)][1], ')')
print('Analytic Result:')
print('eps^-2 : 0.2')
print('eps^-1 : -0.43733071546712660 + 0.62831853071795796j')
print('eps^0 : -0.67330846006628431 - 1.3739149629006919j')
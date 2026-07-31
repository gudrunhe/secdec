from pySecDec.integral_interface import DistevalLibrary

# load the libraries
box1L_p3sq_m4 = DistevalLibrary("box1L_p3sq_m4/disteval/box1L_p3sq_m4.json", verbose=False)
box1L_p3sq_m3 = DistevalLibrary("box1L_p3sq_m3/disteval/box1L_p3sq_m3.json", verbose=False)
box1L_p3sq_m2 = DistevalLibrary("box1L_p3sq_m2/disteval/box1L_p3sq_m2.json", verbose=False)
box1L_p3sq_m1 = DistevalLibrary("box1L_p3sq_m1/disteval/box1L_p3sq_m1.json", verbose=False)
box1L_p3sq_m0 = DistevalLibrary("box1L_p3sq_m0/disteval/box1L_p3sq_m0.json", verbose=False)

s=5
t = -s/2
p3sq = 12/23
m1sq = m2sq = m3sq = m4sq = 1.0

# integrate all configurations and store the results
result = box1L_p3sq_m4(parameters={"s": s, 't':t,"m1sq": m1sq, "m2sq": m2sq, "m3sq":m3sq, "m4sq":m4sq,  "p3sq": p3sq},
               epsrel=1e-4, epsabs=1e-10, format="json", verbose=True)
values_m4 = result["sums"]["box1L_p3sq_m4"]

result = box1L_p3sq_m3(parameters={"s": s, 't':t, "m1sq" : m1sq, "m2sq" : m2sq, "m3sq" : m3sq, "p3sq": p3sq},
               epsrel=1e-4, epsabs=1e-10, format="json", verbose=True)
values_m3 = result["sums"]["box1L_p3sq_m3"]

result = box1L_p3sq_m2(parameters={"s": s, 't':t,"m1sq": m1sq, "m2sq": m2sq, "p3sq": p3sq},
               epsrel=1e-4, epsabs=1e-10, format="json", verbose=True)
values_m2 = result["sums"]["box1L_p3sq_m2"]

result = box1L_p3sq_m1(parameters={"s": s, 't':t,"m1sq": m1sq, "p3sq": p3sq},
               epsrel=1e-4, epsabs=1e-10, format="json", verbose=True)
values_m1 = result["sums"]["box1L_p3sq_m1"]

result = box1L_p3sq_m0(parameters={"s" : s, 't' : t, "p3sq" : p3sq},
               epsrel=1e-4, epsabs=1e-10, format="json", verbose=True)
values_m0 = result["sums"]["box1L_p3sq_m0"]


# Print and comapare to analytic results
# Analytic results at s=5, t=-2.5, p3sq=12/23, msq=1
print('\n4-mass:')
print('Numerical Result:')
print('eps^0 :', values_m4[( 0,)][0], '+/- (', values_m4[( 0,)][1], ')')
print('Analytic Result:')
print('eps^0 : 0.25997725292886648+0.28745168167449897j')
print("-------------------------")

print('3-mass:')
print('Numerical Result:')
print('eps^0 :', values_m3[( 0,)][0], '+/- (', values_m3[( 0,)][1], ')')
print('Analytic Result:')
print('eps^0 : -0.06671868679726986051+0.41479836747508014j')
print("-------------------------")

print('2-mass:')
print('Numerical Result:')
print('eps^-1 :', values_m2[( -1,)][0], '+/- (', values_m2[( -1,)][1], ')')
print('eps^0 :', values_m2[( 0,)][0], '+/- (', values_m2[( 0,)][1], ')')
print('Analytic Result:')
print('eps^-1 : 0.23322250532480140-0.21698687997767321j')
print('eps^0 : 0.28798696685957281+0.53764044956477297j')
print("-------------------------")

print('1-mass:')
print('Numerical Result:')
print('eps^-1 :', values_m1[( -1,)][0], '+/- (', values_m1[( -1,)][1], ')')
print('eps^0 :', values_m1[( 0,)][0], '+/- (', values_m1[( 0,)][1], ')')
print('Analytic Result:')
print('eps^-1 : 0.20689919150479952')
print('eps^0 : -0.40122122266033350+0.34416368817352427j')
print("-------------------------")

print('0-mass:')
print('Numerical Result:')
print('eps^-2 :', values_m0[( -2,)][0], '+/- (', values_m0[( -2,)][1], ')')
print('eps^-1 :', values_m0[( -1,)][0], '+/- (', values_m0[( -1,)][1], ')')
print('eps^0 :', values_m0[( 0,)][0], '+/- (', values_m0[( 0,)][1], ')')
print('Analytic Result:')
print('eps^-2 : -0.16')
print('eps^-1 : 0.60056510005615005')
print('eps^0 : -0.36818992552308344+0.88287303864315181j')
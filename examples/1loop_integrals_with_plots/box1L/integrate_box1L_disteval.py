from pySecDec.integral_interface import DistevalLibrary

# load the library
name = "box1L_m4"
box1L = DistevalLibrary(f'{name}/disteval/{name}.json', verbose=False)

s=5
t = -s/2
p3sq = 12/23
m1sq = m2sq = m3sq = m4sq = 1.0

# integrate
result = box1L(parameters={"s": s, 't':t,"m1sq": m1sq, "m2sq": m2sq, "m3sq":m3sq, "m4sq":m4sq,  "p3sq": p3sq},
               epsrel=1e-4, epsabs=1e-10, format="json", verbose=True)
values = result["sums"][name]


print('Numerical Result:')
print('eps^0 :', values[( 0,)][0], '+/- (', values[( 0,)][1], ')')

if name == "box1L_m4":
    print('Analytic Result (s=5, t=-2.5, p3sq, msq=1):')
    print('eps^0 : 0.25997725292886648+0.28745168167449897j')
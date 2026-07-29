from pySecDec.integral_interface import DistevalLibrary

name = 'triangle1L_m3'
triangle1L = DistevalLibrary(f'{name}/disteval/{name}.json', verbose=False)

result = triangle1L(parameters={"s": 5.0, "m1sq": 1.0, "m2sq": 1.0, "m3sq": 1.0},
               epsrel=1e-6, epsabs=1e-12, format="json", verbose=True)
values = result["sums"][name]

print('Numerical Result:')
print('eps^-0', values[(0,)][0], '+/- (', values[(0,)][1], ')')

if name == 'triangle1L_m3':
    print('Analytic Result (s=5, m1sq=m2sq=m3sq=1):')
    print('eps^0 : -0.89433451187805812 - 0.60470861377111462j')

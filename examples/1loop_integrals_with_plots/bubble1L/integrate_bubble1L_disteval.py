from pySecDec.integral_interface import DistevalLibrary

name = 'bubble1L_m2'
bubble1L = DistevalLibrary(f'{name}/disteval/{name}.json', verbose=False)

result = bubble1L(parameters={"s": 5.0, "m1sq": 1.0, "m2sq": 1.0},
               epsrel=1e-7, epsabs=1e-10, format="json", verbose=True)
values = result["sums"][name]

print('Numerical result:')
print('eps^0', values[(0,)][0], '+/- (', values[(0,)][1], ')')
print('eps^-1', values[(-1,)][0], '+/- (', values[(-1,)][1], ')')

if name == 'bubble1L_m2':
    print('Analytic Result (s=5, m1sq=m2sq=1)')
    print('eps^0 : 0.99237539413446318 + 1.4049629462081437j')
    print('eps^-1: 1+0j')
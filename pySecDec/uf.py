"""The U, F routines"""

import sympy as sp
    
def uf(str_ls,str_props):
    r'''
    Construct the 1st (U) and 2nd (F) Symanzik Polynomials from a list of loop momenta and propagators.
    Outputs a tuple containing (U,F).
    
    :param str_ls:
       list of strings;
       The loop momenta
       
    :param str_props:
       list of strings:
       The propagators
    
    Adapted from A.V. Smirnov's Mathematica UF routine: http://science.sander.su/Tools-UF.htm
    '''
    ls = sp.sympify(str_ls)
    props = sp.sympify(str_props)

    # Feynman parameters
    x = sp.symbols('x')
    x = [sp.symbols('x%i' % i) for i,_ in enumerate(props)]
    
    u = 1
    f = -sum(prop*x[i] for i,prop in enumerate(props))

    for l in ls:
        t0, t1, t2 = reversed(sp.Poly(f,l).all_coeffs())
        u *= t2
        f = sp.together((t0*t2-t1**2/4)/t2)
    
    f = sp.ratsimp(-u*f) # todo: need minus sign?
    u = sp.ratsimp(u)

    return (u,f)

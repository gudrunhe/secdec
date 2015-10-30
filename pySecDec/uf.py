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
    But with conventions of Eq(8) arXiv:0803.4177
    '''
    ls = sp.sympify(str_ls)
    props = sp.sympify(str_props)

    # Feynman parameters
    #x = sp.symbols('x')
    x = [sp.symbols('x%i' % i) for i,_ in enumerate(props)]
    
    u = 1
    f = sum(prop*x[i] for i,prop in enumerate(props))

    for l in ls:
        t0, t1, t2 = reversed(sp.Poly(f,l).all_coeffs())
        u *= t2
        f = sp.together((t0*t2-t1**2/4)/t2)
    
    f = sp.ratsimp(-u*f) # todo: need minus sign?
    u = sp.ratsimp(u)

    return (u,f)


def uf2(str_ls,str_props):
    r'''
    Alternative method for construct the 1st (U) and 2nd (F) Symanzik Polynomials from a list of loop momenta and propagators.
    Outputs a tuple containing (U,F).
        
        :param str_ls:
        list of strings;
        The loop momenta
        
        :param str_props:
        list of strings:
        The propagators
    
    Adapted from Eq(8) arXiv:0803.4177
    '''
    ls = sp.sympify(str_ls)
    props = sp.sympify(str_props)

    # Feynman parameters
    x = [sp.symbols('x%i' % i) for i,_ in enumerate(props)]

    propsum = sum(prop*x[i] for i,prop in enumerate(props)).expand()
    m = sp.Matrix( len(ls), len(ls),
                  lambda i,j: propsum.coeff(ls[i]*ls[j]) if i==j
                  else propsum.coeff(ls[i]*ls[j])/sp.sympify(2) )
    q = sp.Matrix( len(ls), 1,
                  lambda i,j: (propsum.coeff(ls[i])/sp.sympify(-2)).subs([(l,0) for l in ls]))
    j = propsum.subs([(l,0) for l in ls])
    am = m.adjugate()
    u = m.det()
    f = ((q.transpose()*am*q)[0,0]-u*j).expand().together()

    return (u,f)

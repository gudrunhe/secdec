from pySecDec import Coefficient

### Coefficients definitions ###

coeff = [
    # M++--
    [
        # bubble (u) coefficient
        Coefficient(['-8*(t-u)'],['-u-t'],['t','u']),
        
        # bubble (t) coefficient
        Coefficient(['-8*(u-t)'],['-u-t'],['t','u']),

        # box6 coefficient
        Coefficient(['-8*(t**2+u**2)'],['-u-t'],['t','u']),
        
        # box8 coefficient
        Coefficient(['-8*3*(2*eps)'],['1'],['t','u'])
    ]
]

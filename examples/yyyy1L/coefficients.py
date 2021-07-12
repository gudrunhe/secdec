from pySecDec.code_writer.sum_package import Coefficient

### Coefficients definitions ###

coeff = [
    # M++--
    [
        # bubble (u) coefficient
        Coefficient(['-8*(t-u)'],['-u-t'],['eps'],['t','u']),
        
        # bubble (t) coefficient
        Coefficient(['-8*(u-t)'],['-u-t'],['eps'],['t','u']),

        # box6 coefficient
        Coefficient(['-8*(t**2+u**2)'],['-u-t'],['eps'],['t','u']),
        
        # box8 coefficient
        Coefficient(['-8*3*(2*eps)'],['1'],['eps'],['t','u'])
    ]
]

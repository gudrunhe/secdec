from pySecDec.code_writer.sum_package import Coefficient

### Coefficients definitions ###

coeff = [ 
    # bubble (u) coefficient
    [Coefficient(['t-u'],['-u-t'],['eps'],['t','u']),
    
    # bubble (t) coefficient
    Coefficient(['u-t'],['-u-t'],['eps'],['t','u']),

    # box coefficient
    Coefficient(['t**2+u**2'],['-u-t'],['eps'],['t','u'])]
]

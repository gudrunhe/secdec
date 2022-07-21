from pySecDec import Coefficient

### Coefficients definitions ###

coeff = [
    # M++--
    [
        # bubble (u) coefficient
        '-8*(t-u)/(-u-t)',
        
        # bubble (t) coefficient
        '-8*(u-t)/(-u-t)',

        # box6 coefficient
        '-8*(t^2+u^2)/(-u-t)',
        
        # box8 coefficient
        '-8*3*(2*eps)'
    ]
]

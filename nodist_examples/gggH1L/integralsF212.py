### My Integrals ####

import pySecDec as psd

### Integral definitions ###

### Sets for Different Replacement Rules  ###

### Set1  =  If we want the Integrals to be function of s12, s23, mHsq and mtsq. We Leave out -> s13  subject to mHsq = s12 + s13 + s23 

set1 = [                    ('p1*p1', '0'),
                            ('p2*p2', '0'),
                            ('p3*p3', '0'),
                            ('p1*p2', 's12/2'),
                            ('p1*p3', '(mhsq - s12 - s23)/2'),
                            ('p2*p3', 's23/2'),
                            ('mt**2','mtsq')
                        ]

### Set2  =  If we want the Integrals to be function of s13, s23, mHsq and mtsq. We Leave out -> s12  subject to mHsq = s12 + s13 + s23

set2 = [                    ('p1*p1', '0'),
                            ('p2*p2', '0'),
                            ('p3*p3', '0'),
                            ('p1*p2', '(mhsq - s13 - s23)/2'),
                            ('p1*p3', 's13/2'),
                            ('p2*p3', 's23/2'),
                            ('mt**2','mtsq')
                        ]

### Set3  =  If we want the Integrals to be function of s12, s13, mHsq and mtsq. We Leave out -> s23  subject to mHsq = s12 + s13 + s23

set3 = [                    ('p1*p1', '0'),
                            ('p2*p2', '0'),
                            ('p3*p3', '0'),
                            ('p1*p2', 's12/2'),
                            ('p1*p3', 's13/2'),
                            ('p2*p3', '(mhsq - s12 - s13)/2'),
                            ('mt**2','mtsq')
                        ]

### Set4  =  If we want the Integrals to be function of s12, s13, s23 and mtsq. We Leave out -> mHsq  ( We just don't define it then ?) 

set4 = [                    ('p1*p1', '0'),
                            ('p2*p2', '0'),
                            ('p3*p3', '0'),
                            ('p1*p2', 's12/2'),
                            ('p1*p3', 's13/2'),
                            ('p2*p3', 's23/2'),
                            ('mt**2','mtsq')
                        ]

### Set  = Common Set of Replacement Rules for the First 9 Integrals 


common_set = set1    # As of Now, use set1 ( Leave out s13 for the first 9 Integrals )


###  Integrals definitions :


# g9(s12,s23,mhsq,mtsq)
li12 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p1)**2 - (mt)**2', '(k1-p1-p2)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [1, 1, 1, 1],

    replacement_rules = common_set,

    dimensionality= '4-2*eps'

    )

# g8(s23,mhsq,mtsq)
li11 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p1)**2 - (mt)**2', '(k1-p1-p2)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [1, 1, 0, 1],

    replacement_rules = common_set,

    dimensionality= '4-2*eps'

    )

# g7(s12,mhsq,mtsq)
li9 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p1)**2 - (mt)**2', '(k1-p1-p2)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [1, 0, 1, 1],

    replacement_rules = common_set,

    dimensionality= '4-2*eps'

    )

# g6(s23,mtsq)
li8 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p1)**2 - (mt)**2', '(k1-p1-p2)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [0, 1, 1, 1],

    replacement_rules = common_set,

    dimensionality= '4-2*eps'

    )

# g5(s12,mtsq)
li6 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p1)**2 - (mt)**2', '(k1-p1-p2)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [1, 1, 1, 0],

    replacement_rules = common_set,

    dimensionality= '4-2*eps'

    )

# g4(mhsq,mtsq)
li5 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p1)**2 - (mt)**2', '(k1-p1-p2)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [2, 0, 0, 1],

    replacement_rules = common_set,

    dimensionality= '4-2*eps'

    )

# g3(s23,mtsq)
li4 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p1)**2 - (mt)**2', '(k1-p1-p2)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [0, 2, 0, 1],

    replacement_rules = common_set,

    dimensionality= '4-2*eps'

    )

# g2(s12,mtsq)
li2 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p1)**2 - (mt)**2', '(k1-p1-p2)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [2, 0, 1, 0],

    replacement_rules = common_set,

    dimensionality= '4-2*eps'

    )

# g1(mtsq)
li1 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p1)**2 - (mt)**2', '(k1-p1-p2)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [2, 0, 0, 0],

    replacement_rules = common_set,

    dimensionality= '4-2*eps'
    )

### Crossed Integrals ###

# g9(s12,s13,mhsq,mtsq)      # Exchange p1 and p2 from g9(s12,s23,mhsq,mtsq)
li13 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p2)**2 - (mt)**2', '(k1-p1-p2)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [1, 1, 1, 1],

    replacement_rules = set3,

    dimensionality= '4-2*eps'

    )

# g9(s23,s13,mhsq,mtsq)      # Exchange p2 and p3 from g9(s12,s23,mhsq,mtsq)
li14 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p1)**2 - (mt)**2', '(k1-p1-p3)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [1, 1, 1, 1],

    replacement_rules = set2,

    dimensionality= '4-2*eps'

    )

# g2(s13,mtsq)        # Exchange p2 and p3 from g2(s12,mtsq)
li3 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p1)**2 - (mt)**2', '(k1-p1-p3)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [2, 0, 1, 0],

    replacement_rules = set2,

    dimensionality= '4-2*eps'
    )

# g5(s13,mtsq)       # Exchange p2 and p3 from g5(s12,mtsq)
li7 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p1)**2 - (mt)**2', '(k1-p1-p3)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [1, 1, 1, 0],

    replacement_rules = set2,

    dimensionality= '4-2*eps'

    )

# g7(s13,mhsq,mtsq)   # Exchange p2 and p3 from g7(s12,mhsq,mtsq)
li10 =   psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3'],

    propagators = ['(k1)**2 - (mt)**2', '(k1-p1)**2 - (mt)**2', '(k1-p1-p3)**2 - (mt)**2 ', '(k1-p1-p2-p3)**2 - (mt)**2 '],
    powerlist = [1, 0, 1, 1],

    replacement_rules = set2,

    dimensionality= '4-2*eps'

    )


### FINAL INTEGRALS LIST ###

I = [ li1, li2, li3, li4, li5, li6, li7, li8, li9, li10, li11, li12, li13, li14 ]

#############################






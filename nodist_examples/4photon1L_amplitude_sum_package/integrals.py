import pySecDec as psd

### Integral definitions ###

I = [

    # one loop bubble (u)
    psd.loop_integral.LoopIntegralFromGraph(
    internal_lines = [[0,[1,2]],[0,[2,1]]],
    external_lines = [['p1',1],['p2',2]],
    replacement_rules = [('p1*p1', 'u'),('p2*p2', 'u'),('p1*p2', 'u')]),

    # one loop bubble (t)
    psd.loop_integral.LoopIntegralFromGraph(
    internal_lines = [[0,[1,2]],[0,[2,1]]],
    external_lines = [['p1',1],['p2',2]],
    replacement_rules = [('p1*p1', 't'),('p2*p2', 't'),('p1*p2', 't')]),

    # one loop box (in 6 dimensions)
    psd.loop_integral.LoopIntegralFromGraph(
    internal_lines = [['0',[1,2]],[0,[2,3]],[0,[3,4]],[0,[4,1]]],
    external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],
    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 0),
                            ('p4*p4', 0),
                            ('p3*p2', 'u/2'),
                            ('p1*p2', 't/2'),
                            ('p1*p4', 'u/2'),
                            ('p1*p3', '-u/2-t/2'),
                            ('p2*p4', '-u/2-t/2'),
                            ('p3*p4', 't/2')
                        ],
    dimensionality= '6-2*eps'
    ),

    # one loop box (in 8 dimensions)
    psd.loop_integral.LoopIntegralFromGraph(
    internal_lines = [['0',[1,2]],[0,[2,3]],[0,[3,4]],[0,[4,1]]],
    external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],
    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 0),
                            ('p4*p4', 0),
                            ('p3*p2', 'u/2'),
                            ('p1*p2', 't/2'),
                            ('p1*p4', 'u/2'),
                            ('p1*p3', '-u/2-t/2'),
                            ('p2*p4', '-u/2-t/2'),
                            ('p3*p4', 't/2')
                        ],
    dimensionality= '8-2*eps'
    )

]

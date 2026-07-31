import pySecDec as psd

if __name__ == "__main__":

    li = psd.LoopIntegralFromGraph(
         internal_lines = [
            ['m1',[1,2]], 
            ['m2',[2,3]], 
            ['m3',[3,1]], 
        ],

        external_lines = [
            ['p1',1], 
            ['p2',2], 
            ['p3',3], 
            ],

        replacement_rules= [
        ('p3', '-p2-p1'),
        ('p1*p1', 0),
        ('p2*p2', 0),
        ('p1*p2', 's/2'),
        ('m1**2', 'm1sq'),
        ('m2**2', 'm2sq'),
        ('m3**2', 'm3sq'),
        ],
    )

    psd.loop_package(
        
        name = 'triangle1L_m3',
        loop_integral= li,
        real_parameters= ('s', 'm1sq', 'm2sq', 'm3sq'),
        requested_orders=[0],

        )
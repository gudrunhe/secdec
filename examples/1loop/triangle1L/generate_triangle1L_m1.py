import pySecDec as psd

if __name__ == "__main__":

    li = psd.LoopIntegralFromGraph(
         internal_lines = [
            ['m1',[1,2]], 
            [0, [2,3]], 
            [0, [3,1]], 
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
        ('m1**2', 'm1sq')
        ],
    )

    psd.loop_package(
        
        name = 'triangle1L_m1',
        loop_integral= li,
        real_parameters= ('s', 'm1sq'),
        requested_orders=[0],

        )
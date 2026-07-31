import pySecDec as psd

if __name__ == "__main__":

    li=psd.LoopIntegralFromGraph(

        internal_lines= [
            [0, [1,2]],
            [0, [2,1]],
        ],

        external_lines=[
            ['p1',1],
            ['p2',2],
        ],

        replacement_rules=[
            ['p2','-p1'],
            ['p1*p1', 's']
        ],
    )

    psd.loop_package(
        
        name = 'bubble1L_m0',
        loop_integral= li,
        real_parameters= ['s'],
        requested_orders=[0],

        )
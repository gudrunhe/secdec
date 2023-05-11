#!/usr/bin/env python3
import pySecDec as psd

if __name__ == "__main__":

    # Section 6.4 from arXiv:2302.08955
    li = psd.loop_integral.LoopIntegralFromGraph(
        internal_lines = [ 
            ("0", (0,1)),
            ("m", (1,2)),
            ("0", (2,6)),
            ("m", (6,3)),
            ("0", (3,4)),
            ("m", (4,5)),
            ("0", (5,0)),
            ("m", (5,6))
        ],
        external_lines = [("p0",0), ("p1",1), ("p2",2), ("p3",3), ("p4",4)],
        powerlist = [1,1,1,1,1,1,1,1],
        dimensionality = "6-2*eps",
        replacement_rules = [
            ("p4", "-p0-p1-p2-p3"),
            ("m*m", "mm"),
            ("p0*p0", "0"),
            ("p1*p1", "mm"),
            ("p2*p2", "mm"),
            ("p3*p3", "mm"),
            ("p0*p1", "(s01-p0*p0-p1*p1)/2"),
            ("p0*p2", "(s02-p0*p0-p2*p2)/2"),
            ("p0*p3", "(s03-p0*p0-p3*p3)/2"),
            ("p1*p2", "(s12-p1*p1-p2*p2)/2"),
            ("p1*p3", "(s13-p1*p1-p3*p3)/2"),
            ("p2*p3", "(s23-p2*p2-p3*p3)/2")
        ]
    )

    psd.loop_package(
        name = "pentabox_offshell",
        loop_integral = li,
        # Note: the natural prefactor is exactly cancelled here.
        additional_prefactor = "1/gamma(2*eps+2)",
        real_parameters = ["mm", "s01", "s02", "s03", "s12", "s13", "s23"],
        requested_orders = [4],
        decomposition_method = "geometric_ku",
        contour_deformation = True,
        processes = 4
    )

import pySecDec as psd
from .leading_order import leading_order

def li(powerlist, dim="4-2*eps"):
    return psd.loop_integral.LoopIntegralFromPropagators(
        loop_momenta=["l1", "l2"],
        external_momenta=["q1", "q2", "p1", "p2"],
        regulator="eps",
        propagators=[
            "(l1)^2",
            "(l2)^2",
            "(l1 - q1)^2",
            "(l2 - q2)^2",
            "(l1 - p1)^2-mt2",
            "(l1 + l2 - p1)^2-mt2",
            "(l1 + p2 - q1)^2-mt2",
            "(l1 + l2 + p2 - q1 - q2)^2-mt2",
            "(l1 + l2)^2",
            "(l1 + q2)^2",
            "(l2 + q1)^2",
        ],
        powerlist=powerlist,
        dimensionality=dim,
        replacement_rules=[
            ("p1*p1", "mt2"),
            ("p1*p2", "mh2/2 - mt2 + x12/2 - x35/2 - x54/2"),
            ("p1*q1", "x12/2 + x23/2 - x54/2"),
            ("p1*q2", "-1/2*x23"),
            ("p2*p2", "mt2"),
            ("p2*q1", "-1/2*x41"),
            ("p2*q2", "x12/2 - x35/2 + x41/2"),
            ("q1*q1", "0"),
            ("q1*q2", "x12/2"),
            ("q2*q2", "0"),
        ],
    )

def test_basic():
    assert leading_order(li([0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0]), decomposition_method="geometric") <= -2
    assert leading_order(li([0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0]), decomposition_method="geometric_ku") <= -2
    assert leading_order(li([0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0]), decomposition_method="iterative") <= -2

def test_basic_dots():
    assert leading_order(li([0, 0, 1, 0, 0, 2, 1, 1, 0, 0, 0]), decomposition_method="geometric") <= -1
    assert leading_order(li([0, 0, 1, 0, 0, 2, 1, 1, 0, 0, 0]), decomposition_method="geometric_ku") <= -1
    assert leading_order(li([0, 0, 1, 0, 0, 2, 1, 1, 0, 0, 0]), decomposition_method="iterative") <= -1

def test_basic_dimensionality():
    assert leading_order(li([0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0], dim="6-2*eps"), decomposition_method="geometric") <= -2
    assert leading_order(li([0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0], dim="6-2*eps"), decomposition_method="geometric_ku") <= -2
    assert leading_order(li([0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0], dim="6-2*eps"), decomposition_method="iterative") <= -2

def test_basic_numerators():
    assert leading_order(li([0, -1, 1, 0, 0, 2, 1, 1, 0, 0, 0]), decomposition_method="geometric") <= -2
    assert leading_order(li([0, -1, 1, 0, 0, 2, 1, 1, 0, 0, 0]), decomposition_method="geometric_ku") <= -2
    assert leading_order(li([0, -1, 1, 0, 0, 2, 1, 1, 0, 0, 0]), decomposition_method="iterative") <= -2

box2L_hfkk
==========

The following examples were provided by Hjalte Frellesvig and Kirill Kudashkin in July 2018. The integrals have internal propagators with mass equal to that of an external leg, due to this, some of the integrals require the use of the `Split=True` option in order for the numerical integration to be well behaved.

box2L_hfkk_a
------------

This integral has an internal mass equal to an external mass but in practice nevertheless does not need `Split=True`.

box2L_hfkk_b
------------

This integral is quasi-finite and therefore does not need `Split=True`.

box2L_hfkk_c
------------

This integral requires `Split=True` in order for the numerical integral to converge well. The use of `Split=True` substantially increases the time required for the numerical integration.
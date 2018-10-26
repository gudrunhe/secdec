Nbox2L_split
============

The following examples were provided by Hjalte Frellesvig and Kirill Kudashkin in July 2018. The integrals have internal propagators with mass equal to that of an external leg, due to this, some of the integrals require the use of the `split=True` option in order for the numerical integration to be well behaved.

Nbox2L_split_a
--------------

This integral has an internal mass equal to an external mass but in practice nevertheless does not need `split=True`.

Nbox2L_split_b
--------------

This integral requires `split=True` in order for the numerical integral to converge well. The use of `split=True` substantially increases the time required for the numerical integration.

Nbox2L_split_c
--------------

This integral is quasi-finite and therefore does not need `split=True`.

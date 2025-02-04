import numpy as np
import sympy

from .. import algebra
from .. import decomposition
from .common import LoopIntegral
from .loop_package import LoopPackage


def decompose_sectors(
    initial_sector: decomposition.Sector, nvars: int, decomposition_method: str
):
    idx1 = list(range(nvars))
    idx2 = list(range(nvars - 1))
    if decomposition_method == "geometric":
        prim_s = decomposition.geometric.Cheng_Wu(initial_sector, nvars - 1)
        for s in decomposition.geometric.geometric_decomposition(prim_s, idx2):
            yield prim_s, s
    elif decomposition_method == "geometric_ku":
        prim_sectors = decomposition.iterative.primary_decomposition(
            initial_sector, idx1
        )
        for prim_s in prim_sectors:
            for s in decomposition.geometric.geometric_decomposition_ku(prim_s, idx2):
                yield prim_s, s
    elif decomposition_method == "iterative":
        prim_sectors = decomposition.iterative.primary_decomposition(
            initial_sector, idx1
        )
        for prim_s in prim_sectors:
            for s in decomposition.iterative.iterative_decomposition(prim_s, idx2):
                yield prim_s, s
    else:
        raise ValueError(f"Unknown decomposition method: {decomposition_method!r}")


def leading_order(
    li: LoopIntegral,
    decomposition_method: str = "geometric",
    faster: bool = False,
) -> int:
    """
    Return the estimated leading expansion order of a loop
    integral in its dimensional regulator.

    Note that this is an upper bound, in the sense that all
    orders of the given integral lower than the returned one
    are zero, but it is not guaranteed that the returned one is
    itself non-zero.

    :param li:
        Loop integral.
    :param decomposition_method:
        Decomposition method, same as in :func:`.loop_package`.
    :param faster:
        If `True`, save some time by ignoring possible term
        cancellations in the numerator.
    """
    nvars = len(li.integration_variables)
    nregs = len(li.regulators)
    reg_eq_0 = {r: 0 for r in li.regulators}
    pkg = LoopPackage(name="psd_tmplo", loop_integral=li)
    n_poly_names = len(pkg.polynomial_names)
    # all_symbols = pkg.integration_variables + pkg.regulators + pkg.polynomial_names
    # pkg polynomials have symbols: [x1, ..., xn, reg1, .., regn, U, F]
    polys = [p.copy() for p in pkg.polynomials_to_decompose]
    other = [p.copy() for p in pkg.other_polynomials]
    orig_sector = decomposition.Sector(polys, other)
    pole_exp = 0
    for prim_sector, sector in decompose_sectors(
        orig_sector, nvars, decomposition_method
    ):
        assert prim_sector.Jacobian.symbols == sector.Jacobian.symbols
        assert len(sector.Jacobian.expolist) == 1
        # U and F, i.e. sector.cast
        # sector ~ prim.Jacobian * sector.Jacobian * Prod_i(sector.poly[i].monomial)
        expolist = prim_sector.Jacobian.expolist[0] + sector.Jacobian.expolist[0]
        for prod in sector.cast:
            assert type(prod) == algebra.Product
            mono, rest = prod.factors
            assert type(mono) == algebra.ExponentiatedPolynomial
            assert len(mono.expolist) == 1
            assert type(rest) == algebra.ExponentiatedPolynomial
            assert rest.has_constant_term()
            exp0 = mono.exponent.subs(reg_eq_0)
            expolist = expolist + mono.expolist[0] * exp0
        assert (expolist[-n_poly_names:] == 0).all()
        expolist = expolist[:-n_poly_names]
        # Numerators, i.e. sector.other
        assert len(sector.cast) >= n_poly_names
        replacements = {}
        for poly_name, prod in zip(pkg.polynomial_names, sector.cast):
            mono, rest = prod.factors
            if faster:
                replacements[poly_name] = sympy.sympify(
                    str(
                        algebra.Polynomial(
                            mono.expolist, mono.coeffs, polysymbols=mono.symbols
                        )
                    )
                )
            else:
                replacements[poly_name] = sympy.sympify(
                    str(
                        algebra.Polynomial(
                            rest.expolist + mono.expolist[0],
                            rest.coeffs,
                            polysymbols=rest.symbols,
                        )
                    )
                )
        for i, poly in enumerate(sector.other):
            assert type(poly) == algebra.Polynomial
            poly2 = sympy.sympify(str(poly)).subs(replacements).expand()
            prod = algebra.Product(
                algebra.Polynomial(
                    np.zeros([1, len(poly.symbols) - n_poly_names], dtype=int),
                    np.array([1]),
                    polysymbols=poly.symbols[:-n_poly_names],
                ),
                algebra.Polynomial.from_expression(
                    str(poly2), poly.symbols[:-n_poly_names]
                ),
            )
            decomposition.refactorize(prod)
            mono, rest = prod.factors
            expolist = expolist + mono.expolist[0]
        ndiv = sum(1 for e in expolist[:nvars] if e < 0)
        ndiv += sum(expolist[nvars:])
        pole_exp = min(pole_exp, -ndiv)
    gamma_exp = 0
    arg = li.Gamma_argument.subs(reg_eq_0)
    if arg.is_integer and int(arg) <= 0:
        pole_exp -= 1
    return pole_exp

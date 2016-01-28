"""The Sector class"""

from ..algebra import Polynomial, ExponentiatedPolynomial, Product
import numpy as np

class Sector(object):
    '''
    Container class for sectors that arise during the
    sector decomposition.

    :param cast:
        iterable of :class:`.algebra.Product` or
        of :class:`.algebra.Polynomial`;
        The polynomials to be cast to the form
        `<monomial> * (const + ...)`

    :param other:
        iterable of :class:`.algebra.Polynomial`, optional;
        All variable transformations are applied to these
        polynomials but it is not attempted to achieve the
        form `<monomial> * (const + ...)`

    :param Jacobian:
        :class:`.algebra.Polynomial` with one term, optional;
        The Jacobian determinant of this sector. If not provided,
        the according unit monomial (1*x0^0*x1^0...) is assumed.

    '''
    def __init__(self, cast, other=[], Jacobian=None):
        cast  = list(cast)
        other = list(other)

        assert cast, "Must have at least one input polynomial"

        # The elements of `cast` may be of type `Polynomial` or `Product`
        try:
            assert len(cast[0].factors) == 2, "Every `Product` must have exactly two factors" # raises AttributeError if type is `Polynomial`
            poly = cast[0].factors[1]
        except AttributeError:
            poly = cast[0]
        self.number_of_variables = poly.expolist.shape[1]

        initial_monomial_factor = Polynomial([ [0]*self.number_of_variables ], [1], poly.polysymbols)

        if Jacobian is not None:
            assert len(Jacobian.coeffs) == 1, "`Jacobian` must be a monomial"
        else:
            Jacobian = initial_monomial_factor


        for item in cast + other + [Jacobian]:
            # explanation see above
            try:
                assert len(item.factors) == 2, "Every `Product` must have exactly two factors" # raises AttributeError if type is `Polynomial`
                # if control reaches this point, assume that `item` has type `Product`
                poly = item.factors[1]
                assert len(item.factors[0].coeffs) == 1, 'The first factor of every `Product` must be a monomial'
            except AttributeError:
                poly = item
            # number of variables must match for all input polynomials
            assert self.number_of_variables == poly.expolist.shape[1], "The number of variables must be equal for all input polynomials"

        self.Jacobian = Jacobian.copy()
        self.other = [poly.copy() for poly in other]

        self.cast = []
        for item in cast:
            if hasattr(item, 'factors'): # expect type `Product`
                self.cast.append(item.copy())
            else: # expect type `Polynomial`
                if hasattr(item, 'exponent'): # `ExponentiatedPolynomial` --> same exponent for the factorizable monomial
                    monomial = ExponentiatedPolynomial(initial_monomial_factor.expolist,
                                                       initial_monomial_factor.coeffs,
                                                       polysymbols=initial_monomial_factor.polysymbols,
                                                       exponent=item.exponent)
                    item = Product(monomial, item)
                else:
                    item = Product(initial_monomial_factor, item)

                # try to factorize only if constructed from `Polynomial`
                refactorize(item)
                self.cast.append(item)

    def __repr__(self):
        return 'Sector(Jacobian=%s, cast=%s, other=%s)' % (self.Jacobian, self.cast, self.other)

    def __str__(self):
        return 'Sector:\nJacobian=%s\ncast=%s\nother=%s' % (self.Jacobian, self.cast, self.other)

    def copy(self):
        "Return a copy of a :class:`.Sector`."
        return Sector(self.cast, self.other, self.Jacobian)

def refactorize(polyprod, parameter=None):
    '''
    In a :class:`.algebra.Product` of
    the form `<monomial> * <polynomial>`, check if
    a parameter in `<polynomial>` can be shifted to
    the `<monomial>`.
    If possible, modify `polyprod` accordingly.

    :param polyprod:
        :class:`.algebra.Product` of the
        form <monomial> * <polynomial>`;
        The product to refactorize.

    :param parameter:
        integer, optional;
        Check only the parameter with this index.
        If not provided, all parameters are checked.

    '''
    expolist_mono = polyprod.factors[0].expolist
    expolist_poly = polyprod.factors[1].expolist

    if parameter is None:
        factorizable_powers = expolist_poly.min(axis=0)
        expolist_mono[:] += factorizable_powers
        expolist_poly[:] -= factorizable_powers
    else:
        factorizable_power = expolist_poly[:,parameter].min()
        expolist_mono[:,parameter] += factorizable_power
        expolist_poly[:,parameter] -= factorizable_power

# -------------------- `hide` and `unhide` --------------------

def hide(polynomial, count):
    '''
    Hide the last `count` variables of a
    polynomial.
    This function is meant to be used
    before instantiating a :class:`.Sector`.
    It splits the ``expolist`` and the
    ``polysymbols`` at the index ``count``.

    .. seealso::
        :func:`.unhide`

    .. warning::
        The `polynomial` is **NOT** copied.

    '''
    class HideContainer(object): pass
    hidden = HideContainer()
    hidden.expolist = polynomial.expolist[:, -count:]
    hidden.polysymbols = polynomial.polysymbols[-count:]
    hidden.coeffs = polynomial.coeffs

    polynomial.number_of_variables -= count
    polynomial.polysymbols = polynomial.polysymbols[:-count]
    polynomial.expolist = polynomial.expolist[:, :-count]
    polynomial.coeffs = np.ones_like(polynomial.coeffs)
    return polynomial, hidden

def unhide(polynomial1, polynomial2):
    '''
    Undo the operation :func:`.hide`; i.e.
    ``unhide(*hide(polynomial))`` is equal
    to ``polynomial``.

    .. seealso::
        :func:`.hide`

    .. warning::
        `polynomial1` is modified in place.

    '''
    polynomial1.number_of_variables += len(polynomial2.polysymbols)
    polynomial1.polysymbols += polynomial2.polysymbols
    polynomial1.expolist = np.hstack([polynomial1.expolist, polynomial2.expolist])
    polynomial1.coeffs = polynomial2.coeffs
    return polynomial1

# -------------------------------------------------------------

"""The Sector class"""

from ..algebra import Polynomial, PolynomialProduct

class Sector(object):
    '''
    Container class for sectors that arise during the
    sector decomposition.

    :param cast:
        iterable of :class:`.algebra.PolynomialProduct` or
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

        # The elements of `cast` may be of type `Polynomial` or `PolynomialProduct`
        try:
            assert len(cast[0].factors) == 2, "Every `PolynomialProduct` must have exactly two factors" # raises AttributeError if type is `Polynomial`
            poly = cast[0].factors[1]
        except AttributeError:
            poly = cast[0]
        self.number_of_variables = poly.expolist.shape[1]

        initial_monomial_factor = Polynomial([ [0]*self.number_of_variables ], [1])

        if Jacobian is not None:
            assert len(Jacobian.coeffs) == 1, "`Jacobian` must be a monomial"
        else:
            Jacobian = initial_monomial_factor


        for item in cast + other + [Jacobian]:
            # explanation see above
            try:
                assert len(item.factors) == 2, "Every `PolynomialProduct` must have exactly two factors" # raises AttributeError if type is `Polynomial`
                # if control reaches this point, assume that `item` has type `PolynomialProduct`
                poly = item.factors[1]
                assert len(item.factors[0].coeffs) == 1, 'The first factor of every `PolynomialProduct` must be a monomial'
            except AttributeError:
                poly = item
            # number of variables must match for all input polynomials
            assert self.number_of_variables == poly.expolist.shape[1], "The number of variables must be equal for all input polynomials"

        self.Jacobian = Jacobian.copy()
        self.other = [poly.copy() for poly in other]

        self.cast = []
        for item in cast:
            if hasattr(item, 'factors'): # expect type `PolynomialProduct`
                self.cast.append(item.copy())
            else: # expect type `Polynomial`
                self.cast.append(PolynomialProduct(initial_monomial_factor, item))

    def __repr__(self):
        return 'Sector(Jacobian=%s, cast=%s, other=%s)' % (self.Jacobian, self.cast, self.other)

    def __str__(self):
        return 'Sector:\nJacobian=%s\ncast=%s\nother=%s' % (self.Jacobian, self.cast, self.other)

    def copy(self):
        "Return a copy of a :class:`.Sector`."
        return Sector(self.cast, self.other, self.Jacobian)

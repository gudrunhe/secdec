"""
Algebra
-------

Implementation of a simple computer algebra system.

"""

from .misc import argsort_2D_array, argsort_ND_array, doc, \
                  cached_property, sympify_expression
import numpy as np
import sympy as sp

class _Expression(object):
    '''
    Abstract base class for all expressions in this
    computer algebra system.

    '''
    # keep track if immutable expression types are simplified
    simplified = False

    # delete default hash function
    __hash__ = None

    def __add__(self, other):
        if not isinstance(other, _Expression):
            other = Polynomial(np.zeros([1,self.number_of_variables], dtype=int), np.array([sympify_expression(other)]), self.symbols, copy=False)
        return Sum(self, other)
    __radd__ = __add__

    def __neg__(self):
        minus_one = Polynomial(np.zeros([1,self.number_of_variables], dtype=int), np.array([-1]), self.symbols, copy=False)
        return Product(minus_one, self)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        if not isinstance(other, _Expression):
            other = Polynomial(np.zeros([1,self.number_of_variables], dtype=int), np.array([sympify_expression(other)]), self.symbols, copy=False)
        return Product(self, other)
    __rmul__ = __mul__

    def __pow__(self, other):
        if not isinstance(other, _Expression):
            other = Polynomial(np.zeros([1,self.number_of_variables], dtype=int), np.array([sympify_expression(other)]), self.symbols, copy=False)
        return Pow(self, other)

    def __rpow__(self, other):
        if not isinstance(other, _Expression):
            other = Polynomial(np.zeros([1,self.number_of_variables], dtype=int), np.array([sympify_expression(other)]), self.symbols, copy=False)
        return Pow(other, self)

    def clear_cache(self):
        'Clear cached `str`.'
        try:
            del self.__dict__['str']
        except KeyError:
            pass

    docstring_of_replace = \
        '''
        Replace a variable in an expression by a number or a
        symbol.
        The entries in all ``expolist`` of the underlying
        :class:`.Polynomial` are set to zero. The coefficients
        are modified according to `value` and the powers
        indicated in the ``expolist``.

        :param expression:
            :class:`._Expression`;
            The expression to replace the variable.

        :param index:
            integer;
            The index of the variable to be replaced.

        :param value:
            number or sympy expression;
            The value to insert for the chosen variable.

        :param remove:
            bool;
            Whether or not to remove the replaced
            parameter from the ``parameters`` in the
            `expression`.

        '''

class Function(_Expression):
    '''
    Symbolic function that can take care of
    parameter transformations.
    It keeps track of all taken derivatives:
    When :meth:`.derive` is called, save the
    multiindex of the taken derivative.

    The derivative multiindices are the keys in
    the dictionary ``self.derivative_tracks``. The
    values are lists with two elements: Its first
    element is the index to derive the derivative
    indicated by the multiindex in the second
    element by, in order to abtain the derivative
    indicated by the key:

    >>> from pySecDec.algebra import Polynomial, Function
    >>> x = Polynomial.from_expression('x', ['x','y'])
    >>> y = Polynomial.from_expression('y', ['x','y'])
    >>> poly = x**2*y + y**2
    >>> func = Function('f', x, y)
    >>> ddfuncd0d1 = func.derive(0).derive(1)
    >>> func
    Function(f( + (1)*x, + (1)*y), derivative_tracks = {(1, 0): [0, (0, 0)], (1, 1): [1, (1, 0)]})
    >>> func.derivative_tracks
    {(1, 0): [0, (0, 0)], (1, 1): [1, (1, 0)]}
    >>> func.compute_derivatives(poly)
    {(1, 0):  + (2)*x*y, (1, 1):  + (2)*x}

    :param symbol:
        string;
        The symbol to be used to represent
        the `Function`.

    :param arguments:
        arbitrarily many :class:`._Expression`;
        The arguments of the `Function`.

    :param copy:
        bool;
        Whether or not to copy the `arguments`.

    '''
#   hidden keyword arguments:
#   :param differentiated_args:
#       2D-array with shape ``(len(self.arguments), self.number_of_variables)``;
#       This is used to cache the derivatives of the arguments.
#
#   :param derivatives:
#       1D-array of length ``self.number_of_variables``;
#       This is used to cache the derivatives.
#
#   :param derivative_multiindex:
#       vector-like array;
#       The derivatives that have been taken. This is used
#       to sort the derivatives.
#
#   :param basename:
#       string;
#       The name of the function without the "d"'s for possibly taken
#       derivatives.
#
#   :param derivative_tracks:
#       dict;
#       The paths to the taken derivatives.

    def __init__(self, symbol, *arguments, **kwargs):
        copy = kwargs.get('copy', True)

        self.symbol = symbol
        self.number_of_arguments = len(arguments)
        self.number_of_variables = arguments[0].number_of_variables
        self.arguments = []
        for arg in arguments:
            assert arg.number_of_variables == self.number_of_variables, 'Must have the same number of variables in all arguments.'
            self.arguments.append(arg.copy() if copy else arg)

        self.derivative_tracks = kwargs.pop('derivative_tracks', {})

        self.basename = kwargs.pop('basename', self.symbol)
        self.derivative_multiindex = kwargs.pop('derivative_multiindex', tuple([0] * self.number_of_variables))

        self.differentiated_args = kwargs.pop('differentiated_args', np.empty((len(self.arguments), self.number_of_variables), dtype=object))
        self.derivatives = kwargs.pop('derivatives', np.empty(self.number_of_variables, dtype=object))
        self.derivative_symbols = kwargs.pop('derivative_symbols', set())
        self.derivative_symbols.add(self.symbol)

    @cached_property
    def str(self):
        outstr_template = self.symbol + '(%s)'
        str_args = ','.join(str(arg) for arg in self.arguments)
        return outstr_template % str_args

    def __str__(self):
        return self.str

    def __repr__(self):
        out = 'Function('
        out += str(self)
        out += ', derivative_tracks = ' + repr(self.derivative_tracks) + ')'
        return out

    def compute_derivatives(self, expression=None):
        '''
        Compute all derivatives of ``expression`` that
        are mentioned in ``self.derivative_tracks``.
        The purpose of this function is to avoid
        computing the same derivatives multiple times.

        :param expression:
            :class:`._Expression`, optional;
            The expression to compute the derivatives of.
            If not provided, the derivatives are shown
            as in terms of the `function`'s derivatives
            ``dfd<index>``.

        '''
        def update(repository, multiindex, recipies, expression):
            '''
            Recursively compute all intermediate derivatives
            of `expression` indicated by `multiindex` using
            `recipies` and add them to the `repository`.
            Return the derivative indicated by `multiindex`.

            '''
            try:
                # nothing to do if already computed
                return repository[multiindex]
            except KeyError: # not computed yet --> compute
                # stopping criterion: the multiindex (0,...,0) is `expression` itself
                if sum(multiindex) == 0:
                    return expression
                else:
                    # must compute the derivative from a lower derivative
                    index, lower_multiindex = recipies[multiindex]
                    required_derivative = update(repository, lower_multiindex, recipies, expression).derive(index).simplify()
                    repository[multiindex] = required_derivative
                    return required_derivative

        if expression is None:
            expression = self

        derivatives = {}
        for multiindex in self.derivative_tracks.keys():
            update(derivatives, multiindex, self.derivative_tracks, expression.simplify())

        return derivatives

    def copy(self):
        "Return a copy of a :class:`.Function`."
        return type(self)(self.symbol, *(arg.copy() for arg in self.arguments), copy=False,
                          differentiated_args=self.differentiated_args, derivatives=self.derivatives,
                          derivative_symbols=self.derivative_symbols, basename=self.basename,
                          derivative_multiindex=self.derivative_multiindex,
                          derivative_tracks=self.derivative_tracks)

    def simplify(self):
        'Simplify the arguments.'
        if self.simplified:
            return self
        self.clear_cache()
        self.arguments = [arg.simplify() for arg in self.arguments]
        self.simplified = True
        return self

    @property
    def symbols(self):
        return self.arguments[0].symbols

    def derive(self, index):
        '''
        Generate the derivative by the parameter indexed `index`.
        The derivative of a function with `symbol` ``f`` by
        some `index` is denoted as ``dfd<index>``.

        :param index:
            integer;
            The index of the paramater to derive by.

        '''
        derivative = self.derivatives[index]
        if derivative is not None:
            return derivative.copy()

        summands = []
        for argindex, arg in enumerate(self.arguments):
            if self.differentiated_args[argindex, index] is None:
                differentiated_arg = arg.derive(index)
                self.differentiated_args[argindex, index] = differentiated_arg
            else:
                differentiated_arg = self.differentiated_args[argindex, index]

            # catch ``differentiated_arg == 0``
            if type(differentiated_arg) is Polynomial and (differentiated_arg.coeffs == 0).all():
                continue

            # generate the multiindex of the requested derivative
            old_multiindex = self.derivative_multiindex
            new_multiindex = list(old_multiindex)
            new_multiindex[argindex] += 1
            new_multiindex = tuple(new_multiindex)
            self.derivative_tracks[new_multiindex] = [argindex, old_multiindex]

            derivative_symbol = 'd' * np.sum(new_multiindex) + self.basename + \
                ''.join( ('d%i' %argindex) * howmany for argindex,howmany in enumerate(new_multiindex) )

            self.derivative_symbols.add(derivative_symbol)
            summands.append(
                            ProductRule(    # chain rule
                                            differentiated_arg,
                                            type(self)(derivative_symbol, *(arg.copy() for arg in self.arguments),
                                                       differentiated_args=self.differentiated_args, copy=False,
                                                       derivative_symbols=self.derivative_symbols, basename=self.basename,
                                                       derivative_multiindex=new_multiindex, derivative_tracks=self.derivative_tracks),
                                            copy=False
                                       )
                           )
        if summands:
            derivative = self.derivatives[index] = Sum(*summands, copy=False)
        else: # return zero packed into a `Polynomial`
            derivative = self.derivatives[index] = Polynomial(np.zeros([1,self.number_of_variables], dtype=int), np.array([0]), self.symbols, copy=False)

        return derivative

    @doc(_Expression.docstring_of_replace)
    def replace(expression, index, value, remove=False):
        arguments = [arg.replace(index, value, remove) for arg in expression.arguments]
        return type(expression)(
                                    expression.symbol, *arguments, copy=False,
                                    derivative_symbols=expression.derivative_symbols,
                                    basename=expression.basename,
                                    derivative_multiindex=expression.derivative_multiindex,
                                    derivative_tracks=expression.derivative_tracks
                               )

class Polynomial(_Expression):
    '''
    Container class for polynomials.
    Store a polynomial as list of lists counting the powers of
    the variables. For example the polynomial "x1**2 + x1*x2" is
    stored as [[2,0],[1,1]].

    Coefficients are stored in a separate list of strings, e.g.
    "A*x0**2 + B*x0*x1" <-> [[2,0],[1,1]] and ["A","B"].

    :param expolist:
        iterable of iterables;
        The variable's powers for each term.

        .. hint::

            Negative powers are allowed.

    :param coeffs:
        1d array-like with numerical or sympy-symbolic
        (see http://www.sympy.org/) content, e.g. [x,1,2]
        where x is a sympy symbol;
        The coefficients of the polynomial.

    :param polysymbols:
        iterable or string, optional;
        The symbols to be used for the polynomial variables
        when converted to string. If a string is passed, the
        variables will be consecutively numbered.

        For example: expolist=[[2,0],[1,1]] coeffs=["A","B"]
         * polysymbols='x' (default) <-> "A*x0**2 + B*x0*x1"
         * polysymbols=['x','y']     <-> "A*x**2 + B*x*y"

    :param copy:
        bool;
        Whether or not to copy the `expolist` and the `coeffs`.

        .. note::
            If copy is ``False``, it is assumed that the
            `expolist` and the `coeffs` have the correct
            type.

    '''
    def __init__(self, expolist, coeffs, polysymbols='x', copy=True):
        if copy:
            self.expolist = np.array(expolist)
            assert len(self.expolist.shape) == 2, 'All entries in `expolist` must have the same length'
            if not np.issubdtype(self.expolist.dtype, np.integer):
                raise TypeError('All entries in `expolist` must be integer.')
        else:
            self.expolist = expolist

        if copy:
            self.coeffs = np.array(coeffs)
            assert len(self.expolist) == len(self.coeffs), \
                '`expolist` (length %i) and `coeffs` (length %i) must have the same length.' %(len(self.expolist),len(self.coeffs))
            assert len(self.coeffs.shape) == 1, '`coeffs` must be one-dimensional'
        else:
            self.coeffs = coeffs

        self.number_of_variables = self.expolist.shape[1]

        if copy and not np.issubdtype(self.coeffs.dtype, np.number):
            parsed_coeffs = []
            for coeff in self.coeffs:
                if isinstance(coeff,_Expression):
                    assert coeff.number_of_variables == self.number_of_variables, 'Must have the same number of variables as the `Polynomial` for all coeffs'
                    parsed_coeffs.append(coeff.copy() if copy else coeff)
                else:
                    parsed_coeffs.append(sympify_expression(coeff))
            self.coeffs = np.array(parsed_coeffs)

        if not copy:
            # always copy the symbols but do not check them if copy is False
            self.polysymbols = list(polysymbols)
        elif isinstance(polysymbols, str):
            self.polysymbols = sympify_expression([polysymbols + str(i) for i in range(self.number_of_variables)])
            for symbol in self.polysymbols:
                assert symbol.is_Symbol, 'All `polysymbols` must be symbols'
        else:
            self.polysymbols = []
            for item in polysymbols:
                if isinstance(item, sp.Symbol):
                    self.polysymbols.append(item)
                else:
                    item = sympify_expression(item)
                    assert item.is_Symbol, 'All `polysymbols` must be symbols'
                    self.polysymbols.append(item)

    @staticmethod
    def from_expression(expression, polysymbols):
        '''
        Alternative constructor.
        Construct the polynomial from an algebraic expression.

        :param expression:
            string or sympy expression;
            The algebraic representation of the polynomial, e.g.
            "5*x1**2 + x1*x2"

        :param polysymbols:
            iterable of strings or sympy symbols;
            The symbols to be interpreted as the polynomial variables,
            e.g. "['x1','x2']".

        '''
        polysymbols = list(polysymbols)

        if not polysymbols:
            raise TypeError("`polysymbols` must contain at least one symbol")

        polysymbols = sympify_expression(polysymbols)

        for symbol in polysymbols:
            if not symbol.is_Symbol:
                raise TypeError("'%s' is not a symbol" % symbol)

        polysymbols_pm = polysymbols + [1/symbol for symbol in polysymbols]

        sympy_poly = sp.poly(sympify_expression(expression), polysymbols_pm)
        expolist = np.array(sympy_poly.monoms())
        expolist = np.subtract(*np.split(expolist, 2, axis=1))
        coeffs = sympy_poly.coeffs()
        return Polynomial(expolist, coeffs, polysymbols)


    def __repr__(self):
        outstr = ''
        for coeff,expolist in zip(self.coeffs,self.expolist):
            outstr += (" + (%s)" % coeff)
            for i,(power,symbol) in enumerate(zip(expolist,self.polysymbols)):
                if power == 0:
                    continue
                elif power == 1:
                    outstr += "*%s" % symbol
                else:
                    outstr += "*%s**%i" %(symbol,power)
        return outstr

    __str__ = __repr__

    def copy(self):
        "Return a copy of a :class:`.Polynomial` or a subclass."
        return type(self)(self.expolist.copy(), self.coeffs.copy(), self.polysymbols, copy=False)

    def has_constant_term(self, indices=None):
        '''
        Return True if the polynomial can be written as:

        .. math::
            const + ...

        Otherwise, return False.

        :param indices:
            list of integers or None;
            The indices of the `polysymbols` to consider.
            If ``None`` (default) all indices are taken
            into account.

        '''
        relevant_exponents = self.expolist if indices is None else self.expolist[:,indices]
        return (relevant_exponents == 0).all(axis=1).any()

    def becomes_zero_for(self, zero_params):
        '''
        Return True if the polynomial becomes zero if the
        parameters passed in `zero_params` are set to zero.
        Otherwise, return False.

        :param zero_params:
            iterable of integers;
            The indices of the parameters to be checked.

        '''
        return (self.expolist > 0)[:,tuple(zero_params)].any(axis=1).all()

    def derive(self, index):
        '''
        Generate the derivative by the parameter indexed `index`.

        :param index:
            integer;
            The index of the paramater to derive by.

        '''
        # derivative by ``x`` --> have ``x`` coded in ``expolist`` and can have ``x`` in coeffs
        # product rule: derivative(<coeff> * x**k) = <coeff> * k * x**(k-1) + derivative(<coeff>) * x**k

        # summand1 = <coeff> * k * x**(k-1)
        # remove terms that are going to be multiplied by zero
        nonzero_coeffs = np.where(self.expolist[:,index] != 0)
        summand1_coeffs = self.coeffs[nonzero_coeffs] * self.expolist[nonzero_coeffs][:,index]

        summand1_expolist = self.expolist[nonzero_coeffs].copy()
        summand1_expolist[:,index] -= 1

        # need at least one term
        if len(summand1_coeffs) == 0:
            summand1 = None
        else:
            summand1 = Polynomial(summand1_expolist, summand1_coeffs, self.polysymbols, copy=False)


        # summand2 = derivative(<coeff>) * x**k
        summand2 = False
        if not np.issubdtype(self.coeffs.dtype, np.number):
            summand2_coeffs = np.zeros_like(self.coeffs)
            for i, coeff in enumerate(self.coeffs):
                if isinstance(coeff, _Expression):
                    summand2_coeffs[i] = coeff.derive(index)
                    summand2 = True


        if summand2:
            summand2 = Polynomial(self.expolist.copy(), summand2_coeffs, self.polysymbols, copy=False)
            if summand1 is None:
                return summand2
            else:
                return summand1 + summand2 # implicit simplify
        else:
            if summand1 is None:
                return Polynomial(np.zeros([1,self.number_of_variables], dtype=int), np.array([0]), self.polysymbols, copy=False)
            else:
                return summand1

    @property
    def symbols(self):
        return self.polysymbols

    def __add__(self, other):
        'addition operator'
        return self._sub_or_add(other, False)
    __radd__ = __add__

    def __sub__(self, other):
        'subtraction operator'
        return self._sub_or_add(other, True)

    def __rsub__(self,other):
        'other - self'
        return (-self) + other

    def _sub_or_add(self, other, sub):
        '''
        Function that implements addition and subtraction.
        The coefficients of `other` are negated if `sub`
        is `True`.

        '''
        if  type(other) is Polynomial:
            assert self.number_of_variables == other.number_of_variables, 'Number of varibales must be equal for both polynomials in +'

            sum_expolist = np.vstack([self.expolist, other.expolist])
            sum_coeffs = np.hstack([self.coeffs, -other.coeffs if sub else other.coeffs])

            return Polynomial(sum_expolist, sum_coeffs, self.polysymbols, copy=False).simplify(deep=False)

        elif np.issubdtype(type(other), np.number) or isinstance(other, sp.Expr):
            new_expolist = np.vstack([[0]*self.number_of_variables, self.expolist])
            new_coeffs = np.append(-other if sub else other, self.coeffs)
            return Polynomial(new_expolist, new_coeffs, self.polysymbols, copy=False).simplify(deep=False)

        else:
            return NotImplemented

    def __mul__(self, other):
        'multiplication operator'
        if  type(other) is Polynomial:
            assert self.number_of_variables == other.number_of_variables, 'Number of varibales must be equal for both factors in *'

            product_expolist = np.vstack([other.expolist + term for term in self.expolist])
            product_coeffs = np.hstack([other.coeffs * term for term in self.coeffs])

            return Polynomial(product_expolist, product_coeffs, self.polysymbols, copy=False).simplify(deep=False)

        elif np.issubdtype(type(other), np.number) or isinstance(other, sp.Expr):
            if other == 1:
                return self.copy()
            elif other == 0:
                return Polynomial(np.zeros([1,self.number_of_variables], dtype=int), np.array([0]), self.polysymbols, copy=False)
            else:
                return Polynomial(self.expolist.copy(), self.coeffs * other, self.polysymbols, copy=False)

        else:
            return NotImplemented

    __rmul__ = __mul__

    def __neg__(self):
        'arithmetic negation "-self"'
        return Polynomial(self.expolist.copy(), -self.coeffs, self.polysymbols, copy=False)

    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            return super(Polynomial, self).__pow__(exponent)

        if exponent == 0:
            return Polynomial(np.zeros([1,self.number_of_variables], dtype=int), np.array([1]), self.polysymbols, copy=False)

        def iterative_pow(polynomial, exponent):
            if exponent == 1:
                return polynomial
            half_exponent = exponent // 2
            outpoly = iterative_pow(polynomial, half_exponent)
            outpoly *= outpoly
            if 2 * half_exponent == exponent: # `exponent` is even
                return outpoly
            # else:
            assert 2 * half_exponent == exponent - 1 # `exponent` is odd
            return outpoly * self

        return iterative_pow(self, exponent)

    def simplify(self, deep=True):
        '''
        Combine terms that have the same exponents of
        the variables.

        :param deep:
            bool;
            If ``True`` (default) call the `simplify`
            method of the coefficients if they are
            of type :class:`._Expression`.

        '''
        # Sort the expolist first, such that identical entries are
        # grouped together
        sort_key = argsort_2D_array(self.expolist)
        self.expolist = self.expolist[sort_key]
        self.coeffs = self.coeffs[sort_key]

        if deep and not np.issubdtype(self.coeffs.dtype, np.number):
            for i in range(len(self.coeffs)):
                if isinstance(self.coeffs[i], _Expression):
                    self.coeffs[i] = self.coeffs[i].simplify()
                    # do not need to keep type `_Expression` if constant
                    if type(self.coeffs[i]) is Polynomial and len(self.coeffs[i].coeffs) == 1 and (self.coeffs[i].expolist == 0).all():
                        self.coeffs[i] = self.coeffs[i].coeffs[0]

        distance_to_nonzero_previous = 1
        for i in range(1,len(self.coeffs)):
            if self.coeffs[i] == 0:
                distance_to_nonzero_previous += 1
                continue
            # search `self.expolist` for the same term
            # since `self.expolist` is sorted, must only compare with the previous nonzero term
            if (self.expolist[i-distance_to_nonzero_previous] == self.expolist[i]).all():
                # add coefficients
                self.coeffs[i] += self.coeffs[i-distance_to_nonzero_previous]
                # mark previous term for removal by setting coefficient to zero
                self.coeffs[i-distance_to_nonzero_previous] = 0
            distance_to_nonzero_previous = 1

        # remove terms with zero coefficient
        nonzero_coeffs = np.where(self.coeffs != 0)
        self.coeffs = self.coeffs[nonzero_coeffs]
        self.expolist = self.expolist[nonzero_coeffs]

        # need at least one term
        if len(self.coeffs) == 0:
            self.coeffs = np.array([0])
            self.expolist = np.zeros([1,self.number_of_variables], dtype=int)

        return self

    @doc(_Expression.docstring_of_replace)
    def replace(expression, index, value, remove=False):
        # coeffs
        if value == 0: # <coeff> * 0**<exponent> = <coeff> if <exponent>==0 else 0
            new_coeffs = np.zeros_like(expression.coeffs)
            nonzero_indicator = np.where(expression.expolist[:,index] == 0)
            new_coeffs[nonzero_indicator] = expression.coeffs[nonzero_indicator]

        elif np.issubdtype(expression.coeffs.dtype, np.number):
            if value == 1: # <coeff> * 1**<something> = <coeff>
                new_coeffs = expression.coeffs.copy()
            else:
                new_coeffs = np.array([value**int(power) for power in expression.expolist[:,index]])
                new_coeffs *= expression.coeffs
        else:
            new_coeffs = np.empty_like(expression.coeffs)
            if value == 1: # <coeff> * 1**<something> = <coeff>
                for i,coeff in enumerate(expression.coeffs):
                    if isinstance(coeff, _Expression):
                        new_coeffs[i] = coeff.replace(index, value, remove)
                    else:
                        new_coeffs[i] = coeff
            else:
                for i,(coeff,power) in enumerate(zip(expression.coeffs, expression.expolist[:,index])):
                    if isinstance(coeff, _Expression):
                        new_coeffs[i] = coeff.replace(index, value, remove) * value**int(power)
                    else:
                        new_coeffs[i] = coeff * value**int(power)

        # exponent (if applicable)
        exponent = None
        if isinstance(expression,ExponentiatedPolynomial):
            # replace in exponent if it has a type `replace` can handle
            if isinstance(expression.exponent, _Expression):
                exponent = expression.exponent.replace(index, value, remove)
            else:
                exponent = expression.exponent

        # expolist
        if remove:
            new_expolist = np.delete(expression.expolist, index, axis=1)
            if index == -1:
                new_polysymbols = list(expression.polysymbols)
                new_polysymbols.pop()
            else:
                new_polysymbols = expression.polysymbols[:index] + expression.polysymbols[index+1:]
        else:
            new_expolist = expression.expolist.copy()
            new_expolist[:,index] = 0
            new_polysymbols = expression.polysymbols

        # final generation of the new polynomial
        if exponent is None:
            return type(expression)(new_expolist, new_coeffs, new_polysymbols, copy=False)
        else:
            return ExponentiatedPolynomial(new_expolist, new_coeffs, exponent, new_polysymbols, copy=False)

    def refactorize(self, *parameters):
        """
        Returns a product of the greatest factor that
        could be pulled out and the factorised polynomial.

        :param parameter:
            arbitrarily many integers;

        """
        prod = Product(Polynomial([np.zeros(len(self.symbols), dtype=int)],[1],self.symbols),self)
        refactorize(prod, *parameters)
        return prod

class ExponentiatedPolynomial(Polynomial):
    '''
    Like :class:`.Polynomial`, but with a global exponent.
    :math:`polynomial^{exponent}`

    :param expolist:
        iterable of iterables;
        The variable's powers for each term.

    :param coeffs:
        iterable;
        The coefficients of the polynomial.

    :param exponent:
        object, optional;
        The global exponent.

    :param polysymbols:
        iterable or string, optional;
        The symbols to be used for the polynomial variables
        when converted to string. If a string is passed, the
        variables will be consecutively numbered.

        For example: expolist=[[2,0],[1,1]] coeffs=["A","B"]
         * polysymbols='x' (default) <-> "A*x0**2 + B*x0*x1"
         * polysymbols=['x','y']     <-> "A*x**2 + B*x*y"

    :param copy:
        bool;
        Whether or not to copy the `expolist`, the `coeffs`,
        and the `exponent`.

        .. note::
            If copy is ``False``, it is assumed that the
            `expolist`, the `coeffs` and the `exponent` have
            the correct type.

    '''
    def __init__(self, expolist, coeffs, exponent=1, polysymbols='x', copy=True):
        Polynomial.__init__(self, expolist, coeffs, polysymbols, copy)
        if np.issubdtype(type(exponent), np.number):
            self.exponent = exponent
        elif isinstance(exponent,_Expression):
            self.exponent = exponent.copy() if copy else exponent
        else:
            self.exponent = sympify_expression(exponent) if copy else exponent

    @staticmethod
    def from_expression(*args,**kwargs):
        raise NotImplementedError('Cannot create `ExponentiatedPolynomial` from expression.')

    __pow__ = _Expression.__pow__
    __add__ = __radd__ = _Expression.__add__
    __neg__ = _Expression.__neg__
    __sub__ = _Expression.__sub__
    __rsub__ = _Expression.__rsub__
    __mul__ = __rmul__ = _Expression.__mul__

    def __repr__(self):
        if self.exponent == 1:
            return super(ExponentiatedPolynomial, self).__repr__()
        else:
            return '(' + super(ExponentiatedPolynomial, self).__repr__() \
                       + ')**(%s)' % self.exponent
    __str__ = __repr__

    def derive(self, index):
        '''
        Generate the derivative by the parameter indexed `index`.

        :param index:
            integer;
            The index of the paramater to derive by.

        '''
        # derive an expression of the form "poly**exponent"
        # derivative(poly**exponent) = poly**exponent*derivative(exponent)*log(poly) + poly**(exponent-1)*exponent*derivative(poly)

        if isinstance(self.exponent, _Expression):
            # summand0: poly**exponent*derivative(exponent)*log(poly)
            summand0_factors = [self.copy()]
            summand0_factors.append(self.exponent.derive(index))
            summand0_factors.append(LogOfPolynomial(self.expolist.copy(), self.coeffs.copy(), self.polysymbols, copy=False))
            summand0 = Product(*summand0_factors, copy=False)
        else:
            summand0 = None

        # summand1: poly**(exponent-1)*derivative(poly)
        # factor1 = "exponent*derivative(poly)"
        # catch ``derivative(poly) = 0``
        derivative_poly = Polynomial(self.expolist.copy(), self.coeffs.copy(), self.polysymbols, copy=False).derive(index)
        if (derivative_poly.coeffs == 0).all():
            summand1 = None
        else:
            factor1 = self.exponent * derivative_poly
            # factor0 = poly**(exponent-1)   -->   simplification: (...)**0 = 1
            # do not need factor 0 in that case
            new_exponent = self.exponent - 1
            if new_exponent == 0:
                # factor0 = 1 in this case
                summand1 = factor1
            else:
                factor0 = ExponentiatedPolynomial(self.expolist.copy(),
                                                  self.coeffs.copy(),
                                                  new_exponent,
                                                  self.polysymbols,
                                                  copy=False)
                summand1 = Product(factor0, factor1, copy=False)

        if summand0 is None:
            if summand1 is None:
                return Polynomial(np.zeros([1,self.number_of_variables], dtype=int), np.array([0]), self.polysymbols, copy=False)
            else:
                return summand1
        else:
            if summand1 is None:
                return summand0
            else:
                return Sum(summand0, summand1, copy=False)

    def copy(self):
        "Return a copy of a :class:`.Polynomial` or a subclass."
        exponent = self.exponent.copy() if isinstance(self.exponent, _Expression) else self.exponent
        return type(self)(self.expolist.copy(), self.coeffs.copy(), exponent, self.polysymbols, copy=False)

    def simplify(self):
        '''
        Apply the identity <something>**0 = 1 or
        <something>**1 = <something> or 1**<something> = 1
        if possible, otherwise call the simplify method of
        the base class. Convert ``exponent`` to symbol if
        possible.

        '''
        if isinstance(self.exponent, _Expression):
            self.exponent = self.exponent.simplify()
            if type(self.exponent) is Polynomial and (self.exponent.expolist == 0).all(): # exponent is a constant --> can take just the coefficient
                self.exponent = self.exponent.coeffs[0]
        if self.exponent == 0:
            self.coeffs = np.array([1])
            self.expolist = np.array([[0]*self.number_of_variables])
            self.exponent = 1
        if self.exponent == 1:
            return Polynomial(self.expolist, self.coeffs, self.polysymbols, copy=False)

        super(ExponentiatedPolynomial, self).simplify()

        # 1**<something> --> simplify to type `Polynomial`
        if len(self.coeffs) == 1 and self.coeffs[0] == 1 and (self.expolist == 0).all():
            return Polynomial(self.expolist, self.coeffs, self.polysymbols, copy=False)

        return self

    def refactorize(self, *parameters):
        """
        Returns a product of the greatest factor that
        could be pulled out and the factorised polynomial.

        :param parameter:
            arbitrarily many integers;

        """
        # remove global exponent
        polynomial_without_exponent = Polynomial(self.expolist, self.coeffs, self.symbols)
        # refactorize
        prod = Product(Polynomial([np.zeros(len(self.symbols), dtype=int)],[1],self.symbols),polynomial_without_exponent)
        refactorize(prod, *parameters)
        poly1,poly2 = prod.factors[0],prod.factors[1]
        # add exponent back on
        expo_poly1 = ExponentiatedPolynomial(poly1.expolist, poly1.coeffs, self.exponent, poly1.symbols)
        expo_poly2 = ExponentiatedPolynomial(poly2.expolist, poly2.coeffs, self.exponent, poly2.symbols)
        return Product(expo_poly1,expo_poly2)

class LogOfPolynomial(Polynomial):
    '''
    The natural logarithm of a :class:`.Polynomial`.

    :param expolist:
        iterable of iterables;
        The variable's powers for each term.

    :param coeffs:
        iterable;
        The coefficients of the polynomial.

    :param exponent:
        object, optional;
        The global exponent.

    :param polysymbols:
        iterable or string, optional;
        The symbols to be used for the polynomial variables
        when converted to string. If a string is passed, the
        variables will be consecutively numbered.

        For example: expolist=[[2,0],[1,1]] coeffs=["A","B"]
         * polysymbols='x' (default) <-> "A*x0**2 + B*x0*x1"
         * polysymbols=['x','y']     <-> "A*x**2 + B*x*y"

    '''
    @staticmethod
    def from_expression(expression, polysymbols):
        '''
        Alternative constructor.
        Construct the :class:`LogOfPolynomial` from an algebraic
        expression.

        :param expression:
            string or sympy expression;
            The algebraic representation of the polynomial, e.g.
            "5*x1**2 + x1*x2"

        :param polysymbols:
            iterable of strings or sympy symbols;
            The symbols to be interpreted as the polynomial variables,
            e.g. "['x1','x2']".

        '''
        poly = Polynomial.from_expression(expression,polysymbols)
        return LogOfPolynomial(poly.expolist, poly.coeffs, poly.polysymbols, copy=False)

    def __repr__(self):
        return 'log(%s)' % Polynomial.__repr__(self)

    __str__ = __repr__

    __pow__ = _Expression.__pow__
    __add__ = __radd__ = _Expression.__add__
    __neg__ = _Expression.__neg__
    __sub__ = _Expression.__sub__
    __rsub__ = _Expression.__rsub__
    __mul__ = __rmul__ = _Expression.__mul__

    def derive(self, index):
        '''
        Generate the derivative by the parameter indexed `index`.

        :param index:
            integer;
            The index of the paramater to derive by.

        '''
        # derive an expression of the form "log(poly)"
        # chain rule: poly**(-1) * derivative(poly)
        #   --> factor1 = "derivative(poly)"
        factor1 = Polynomial(self.expolist.copy(), self.coeffs.copy(), self.polysymbols, copy=False)
        factor1 = factor1.derive(index)
        # catch ``factor1 == 0``
        if (factor1.coeffs == 0).all():
            return Polynomial(np.zeros([1,len(self.polysymbols)], dtype=int), np.array([0]), self.polysymbols, copy=False)

        #   --> factor0 = "poly**(-1)"
        factor0 = ExponentiatedPolynomial(self.expolist.copy(),
                                          self.coeffs.copy(),
                                          -1,
                                          self.polysymbols,
                                          copy=False)

        return Product(factor0, factor1, copy=False)

    def simplify(self):
        '''
        Apply the identity ``log(1) = 0``, otherwise call
        the simplify method of the base class.

        '''
        super(LogOfPolynomial, self).simplify()
        if len(self.coeffs) == 1 and self.coeffs[0] == 1 and (self.expolist == 0).all():
            return Polynomial(np.zeros([1,self.number_of_variables], dtype=int), np.array([0]), self.polysymbols, copy=False)
        else:
            return self

class Sum(_Expression):
    r'''
    Sum of polynomials.
    Store one or polynomials :math:`p_i` to be interpreted as
    product :math:`\sum_i p_i`.

    :param summands:
        arbitrarily many instances of :class:`.Polynomial`;
        The summands :math:`p_i`.

    :param copy:
        bool;
        Whether or not to copy the `summands`.

    :math:`p_i` can be accessed with ``self.summands[i]``.

    Example:

    .. code-block:: python

        p = Sum(p0, p1)
        p0 = p.summands[0]
        p1 = p.summands[1]

    '''
    def __init__(self,*summands, **kwargs):
        copy = kwargs.get('copy', True)

        self.summands = [summand.copy() if copy else summand for summand in summands]
        assert self.summands, 'Must have at least one summand'

        self.number_of_variables = self.summands[0].number_of_variables

        for summand in self.summands:
            if summand.number_of_variables != self.number_of_variables:
                raise TypeError('Must have the same number of variables for all summands.')

    @cached_property
    def str(self):
        stringified_summands = []
        for summand in self.summands:
            stringified_summands.append( '(' + str(summand) + ')' )
        return ' + '.join(stringified_summands)

    def __repr__(self):
        return self.str

    __str__ = __repr__

    def simplify(self):
        '''
        If one or more of ``self.summands`` is a
        :class:`Sum`, replace it by its summands.
        If only one summand is present, return that summand.
        Remove zero from sums.

        '''
        if self.simplified:
            return self

        self.clear_cache()

        changed = True
        while changed:
            changed = False
            old_summands = self.summands
            self.summands = []
            for summand in old_summands:
                summand = summand.simplify()
                if isinstance(summand, Sum):
                    changed = True
                    self.summands.extend(summand.summands)
                elif isinstance(summand, Polynomial):
                    if (summand.coeffs == 0).all():
                        changed = True
                        zero = summand
                    else:
                        self.summands.append(summand)
                else:
                    self.summands.append(summand)
        if len(self.summands) == 1:
            return self.summands[0]
        elif len(self.summands) == 0:
            self.summands = [zero]
            return zero
        else:
            self.simplified = True
            return self

    @property
    def symbols(self):
        return self.summands[0].symbols

    def copy(self):
        "Return a copy of a :class:`.Sum`."
        return Sum(*(summand.copy() for summand in self.summands), copy=False)

    def derive(self, index):
        '''
        Generate the derivative by the parameter indexed `index`.

        :param index:
            integer;
            The index of the paramater to derive by.

        '''
        # derivative(p1 + p2 + ...) = derivative(p1) + derivative(p2) + ...
        return Sum(*(summand.derive(index) for summand in self.summands), copy=False)

    @doc(_Expression.docstring_of_replace)
    def replace(expression, index, value, remove=False):
        outsummands = []
        for summand in expression.summands:
            outsummands.append(summand.replace(index,value,remove))
        return Sum(*outsummands, copy=False)

class Product(_Expression):
    r'''
    Product of polynomials.
    Store one or polynomials :math:`p_i` to be interpreted as
    product :math:`\prod_i p_i`.

    :param factors:
        arbitrarily many instances of :class:`.Polynomial`;
        The factors :math:`p_i`.

    :param copy:
        bool;
        Whether or not to copy the `factors`.

    :math:`p_i` can be accessed with ``self.factors[i]``.

    Example:

    .. code-block:: python

        p = Product(p0, p1)
        p0 = p.factors[0]
        p1 = p.factors[1]


    '''
    def __init__(self,*factors, **kwargs):
        copy = kwargs.get('copy', True)

        self.factors = [factor.copy() if copy else factor for factor in factors]
        assert self.factors, 'Must have at least one factor'

        self.number_of_variables = self.factors[0].number_of_variables

        for factor in self.factors:
            if factor.number_of_variables != self.number_of_variables:
                raise TypeError('Must have the same number of variables for all factors.')

    @cached_property
    def str(self):
        stringified_factors = []
        for factor in self.factors:
            stringified_factors.append( '(' + str(factor) + ')' )
        return ' * '.join(stringified_factors)

    def __repr__(self):
        return self.str

    __str__ = __repr__

    def copy(self):
        "Return a copy of a :class:`.Product`."
        return Product(*(factor.copy() for factor in self.factors), copy=False)

    def simplify(self):
        '''
        If one or more of ``self.factors`` is a
        :class:`.Product`, replace it by its factors.
        If only one factor is present, return that factor.
        Remove factors of one and zero.

        '''
        if self.simplified:
            return self

        self.clear_cache()

        changed = True
        while changed:
            changed = False
            old_factors = self.factors
            self.factors = []
            for factor in old_factors:
                factor = factor.simplify()
                if isinstance(factor, Product):
                    changed = True
                    self.factors.extend(factor.factors)
                elif type(factor) is Polynomial:
                    if (factor.expolist == 0).all() and (factor.coeffs == 1).all():
                        one = factor
                        changed = True
                    elif (factor.coeffs == 0).all():
                        factor.simplify()
                        self.factors = [factor]
                        break
                    else:
                        self.factors.append(factor)
                else:
                    self.factors.append(factor)
        if len(self.factors) == 1:
            return self.factors[0]
        elif len(self.factors) == 0:
            self.factors = [one]
            return one
        else:
            self.simplified = True
            return self

    @property
    def symbols(self):
        return self.factors[0].symbols

    def derive(self, index):
        '''
        Generate the derivative by the parameter indexed `index`.
        Return an instance of the optimized :class:`.ProductRule`.

        :param index:
            integer;
            The index of the paramater to derive by.

        '''
        return ProductRule(*self.factors, copy=False).derive(index)

    @doc(_Expression.docstring_of_replace)
    def replace(expression, index, value, remove=False):
        outfactors = []
        for factor in expression.factors:
            outfactors.append(factor.replace(index,value,remove))
        return Product(*outfactors, copy=False)

class ProductRule(_Expression):
    r'''
    Store an expression of the form

    .. math::
        \sum_i c_i \prod_j \prod_k
        \left( \frac{d}{dx_k} \right)^{n_{ijk}}
        f_j \left( \lbrace x_k \rbrace \right)

    The main reason for introducing this class is a speedup
    when calculating derivatives. In particular, this class
    implements simplifications such that the number of
    terms grows less than exponentially (scaling of the
    naive implementation of the product rule) with the
    number of derivatives.

    :param expressions:
        arbitrarily many expressions;
        The expressions :math:`f_j`.

    '''
    def __init__(self, *expressions, **kwargs):
        # If keyword argument `internal_regenerate` is passed and True, use passed `factorlist`, `coeffs`, and `expressions` instead of regenerating
        copy = kwargs.get('copy', True)

        if kwargs.get('internal_regenerate', False):
            self.factorlist = kwargs['factorlist'].copy() if copy else kwargs['factorlist']
            self.coeffs = kwargs['coeffs'].copy() if copy else kwargs['coeffs']
            if copy:
                self.expressions = []
                for expression_derivatives in kwargs['expressions']:
                    this_expression = {}
                    for k,(key,value) in enumerate(expression_derivatives.items()):
                        this_expression[key] = value.copy()
                    self.expressions.append(this_expression)

                # ``value`` is still one of the derivatives from the for loop
                self.symbols = list(value.symbols)

            else:
                self.expressions = kwargs['expressions']
                for expression in self.expressions[0].values():
                    self.symbols = list(expression.symbols)
                    break

            self.number_of_variables = len(self.symbols)

        else:
            self.number_of_variables = expressions[0].number_of_variables

            # check input consistency
            for expression in expressions:
                if expression.number_of_variables != self.number_of_variables:
                    raise TypeError('Must have the same number of variables for all expressions.')

            # store the `expressions` and (later) its derivatives in ``self.expressions``
            # outer list: index ``j`` of the input expressions
            # inner dict: use ``n_k`` as keys
            # automatically simplify the cache
            self.expressions = [{
                                    tuple([0]*self.number_of_variables) :
                                    expression.copy().simplify() if copy else expression.simplify()
                                } for expression in expressions]

            # The `factorlist` is a 3 dimensional array. Its first
            # index denotes the terms of the sum :math:`i`. The
            # second index denotes the factor :math:`j`. The last
            # index denotes the varible to take the derivative with
            # respect to. The value of ``factorlist[i,j,k]`` denotes
            # how many derivatives :math:`n_{ijk}` are to be taken.
            # If not given, a product of the `expressions` without
            # derivatives are assumed.

            # initialize the `factorlist`:
            #  - all `n_ijk` zeros since we do not have derivatives yet --> np.zeros
            #  - shape ``(1,len(expressions),self.number_of_variables)`` --> one term in the sum with ``len(expressions)``
            #                                                                factors depending on ``self.number_of_variables``
            #                                                                variables
            self.factorlist = np.zeros((1,len(expressions),self.number_of_variables), dtype=int)

            # initialize the coeffs ``c_i`` with ones
            self.coeffs = np.array([1])

            self.symbols = list(expressions[0].symbols)

    @cached_property
    def str(self):
        outstr = ''
        for i,(coeff,n_jk) in enumerate(zip(self.coeffs,self.factorlist)):
            if coeff != 0:
                outstr += (" + (%i)" % coeff)
                for j,n_k in enumerate(n_jk):
                    outstr += ' * (%s)' % str(self.expressions[j][tuple(n_k)])
        return outstr if outstr else ' + (0)'

    def __repr__(self):
        return self.str

    __str__ = __repr__

    def derive(self, index):
        '''
        Generate the derivative by the parameter indexed `index`.
        Note that this class is particularly designed to hold
        derivatives of a product.

        :param index:
            integer;
            The index of the paramater to derive by.

        '''
        # product rule: derivative(<coeff> * x**k) = <coeff> * k * x**(k-1) + derivative(<coeff>) * x**k

        # generate new `factorlist`
        new_factorlist = []
        for j in range(len(self.expressions)):
            # from product rule: increase ``n_ijk`` for every ``i``
            factorlist = self.factorlist.copy()
            factorlist[:,j,index] += 1
            new_factorlist.append(factorlist)
        new_factorlist = np.vstack(new_factorlist)
        new_coeffs = np.hstack([self.coeffs]*len(self.expressions))

        # generate missing derivatives
        # do not make a copy since it does not hurt having the child point to the same ``expressions``
        for n, term in enumerate(new_factorlist):
            for derivative_multiindex, expression in zip(term, self.expressions):
                derivative_multiindex = tuple(derivative_multiindex)
                try:
                    expression[derivative_multiindex]
                except KeyError: # need a derivative that is not calculated yet
                    lower_derivative_multiindex = list(derivative_multiindex)
                    lower_derivative_multiindex[index] -= 1
                    lower_derivative_multiindex = tuple(lower_derivative_multiindex)
                    lower_derivative = expression[lower_derivative_multiindex]
                    expression[derivative_multiindex] = lower_derivative.derive(index).simplify() # automatically simplify cache
                # set the coefficent to zero if the derivative is zero, so that it doesn't get outputted in str()
                if new_coeffs[n] != 0 and isinstance(expression[derivative_multiindex], Polynomial) and not np.any(expression[derivative_multiindex].expolist):
                    if sympify_expression(expression[derivative_multiindex]).simplify() == 0:
                        new_coeffs[n] = 0


        return ProductRule(internal_regenerate=True, copy=False,
                           factorlist=new_factorlist, coeffs=new_coeffs,
                           expressions=self.expressions)

    def copy(self):
        "Return a copy of a :class:`.ProductRule`."
        return ProductRule(copy=True, internal_regenerate=True,
                           factorlist=self.factorlist,
                           coeffs=self.coeffs, expressions=self.expressions)

    def simplify(self):
        '''
        Combine terms that have the same derivatives
        of the `expressions`.

        '''
        if self.simplified:
            return self

        self.clear_cache()

        # Sort the factorlist first, such that identical entries are
        # grouped together
        sort_key = argsort_ND_array(self.factorlist)
        self.factorlist = self.factorlist[sort_key]
        self.coeffs = self.coeffs[sort_key]

        # find zeros in `self.expressions`
        for j, derivative_multiindex in enumerate(self.factorlist[0]):
            derivative_multiindex = tuple(derivative_multiindex)
            expression = self.expressions[j][derivative_multiindex]
            if type(expression) is Polynomial and (expression.coeffs == 0).all():
                self.coeffs[0] = 0

        for i in range(1,len(self.coeffs)):
            if self.coeffs[i] == 0: continue
            previous_term = self.factorlist[i-1]
            # search `self.factorlist` for the same term
            # since `self.factorlist` is sorted, must only compare with the previous term
            if (previous_term == self.factorlist[i]).all():
                # add coefficients
                self.coeffs[i] += self.coeffs[i-1]
                # mark previous term for removal by setting coefficient to zero
                self.coeffs[i-1] = 0

            # find zeros in `self.expressions`
            for j, derivative_multiindex in enumerate(self.factorlist[i]):
                derivative_multiindex = tuple(derivative_multiindex)
                expression = self.expressions[j][derivative_multiindex]
                if type(expression) is Polynomial and (expression.coeffs == 0).all():
                    self.coeffs[i] = 0

        # remove terms with zero coefficient
        nonzero_coeffs = np.where(self.coeffs != 0)
        self.coeffs = self.coeffs[nonzero_coeffs]
        self.factorlist = self.factorlist[nonzero_coeffs]

        # need at least one term
        if len(self.coeffs) == 0:
            self.coeffs = np.array([0])
            self.factorlist = np.zeros((1,self.factorlist.shape[1],self.number_of_variables), dtype=int)
            # simplify to zero if change of type possible
            return Polynomial(np.zeros([1,len(self.symbols)], dtype=int), np.array([0]), self.symbols, copy=False)

        self.simplified = True
        return self

    def to_sum(self):
        'Convert the :class:`.ProductRule` to :class:`.Sum`'
        summands = []
        for coeff, term in zip(self.coeffs, self.factorlist):
            coeff = Polynomial(np.zeros([1,self.number_of_variables], dtype=int), np.array([int(coeff)]), self.symbols, copy=False)
            factors = []
            for j, derivative_multiindex in enumerate(term):
                factors.append(self.expressions[j][tuple(derivative_multiindex)])
            summands.append(Product(coeff, *factors, copy=False))
        return Sum(*summands, copy=False)

    @doc(_Expression.docstring_of_replace)
    def replace(self, index, value, remove=False):
        summands = []
        for coeff, term in zip(self.coeffs, self.factorlist):
            factors = []
            for j, derivative_multiindex in enumerate(term):
                factors.append(self.expressions[j][tuple(derivative_multiindex)].replace(index, value, remove))
            summands.append(Product(*factors, copy=False) * int(coeff))
        return Sum(*summands, copy=False)

class Pow(_Expression):
    r'''
    Exponential.
    Store two expressions ``A`` and ``B`` to be interpreted
    as the exponential ``A**B``.

    :param base:
        :class:`._Expression`;
        The base ``A`` of the exponential.

    :param exponent:
        :class:`._Expression`;
        The exponent ``B``.

    :param copy:
        bool;
        Whether or not to copy `base` and `exponent`.

    '''
    def __init__(self, base, exponent, copy=True):
        if base.number_of_variables != exponent.number_of_variables:
            raise TypeError('Must have the same number of variables for `base` and `exponent`.')

        self.number_of_variables = exponent.number_of_variables
        self.base = base.copy() if copy else base
        self.exponent = exponent.copy() if copy else exponent

    @cached_property
    def str(self):
        return '(' + str(self.base) + ') ** (' + str(self.exponent) + ')'

    def __repr__(self):
        return self.str

    __str__ = __repr__

    def copy(self):
        "Return a copy of a :class:`.Pow`."
        return Pow(self.base.copy(), self.exponent.copy(), copy=False)

    @property
    def symbols(self):
        return self.base.symbols

    def simplify(self):
        '''
        Apply the identity <something>**0 = 1
        or <something>**1 = <something> or
        1**<something> = 1 if possible. Convert
        to :class:`.ExponentiatedPolynomial` or
        :class:`.Polynomial` if possible.

        '''
        if self.simplified:
            return self

        self.clear_cache()

        self.base = self.base.simplify()

        if type(self.base) is Polynomial: # need exact type `Polynomial` for this, not subtype
            return ExponentiatedPolynomial(self.base.expolist, self.base.coeffs, self.exponent, self.base.polysymbols, copy=False).simplify()

        self.exponent = self.exponent.simplify()

        if type(self.exponent) is Polynomial:
            if (self.exponent.coeffs==0).all():
                symbols = self.exponent.polysymbols
                return Polynomial(np.zeros([1,len(symbols)], dtype=int), np.array([1]), symbols, copy=False)
            elif len(self.exponent.coeffs)==1 and (self.exponent.coeffs==1).all() and (self.exponent.expolist==0).all():
                return self.base

        self.simplified = True
        return self

    def derive(self, index):
        '''
        Generate the derivative by the parameter indexed `index`.

        :param index:
            integer;
            The index of the paramater to derive by.

        '''
        # derive an expression of the form "base**exponent"
        # derivative(base**exponent) = base**exponent*derivative(exponent)*log(base) + base**(exponent-1)*exponent*derivative(poly)

        # summand0: base**exponent*derivative(exponent)*log(base)
        summand0_factors = [self.copy()]
        summand0_factors.append(self.exponent.derive(index))
        summand0_factors.append(Log(self.base.copy(), copy=False))
        summand0 = Product(*summand0_factors, copy=False)

        # summand1: base**(exponent-1)*derivative(base)
        # factor0 = base**(exponent-1)   -->   simplification: (...)**0 = 1
        # do not need factor 0 in that case
        new_exponent = self.exponent - 1
        if new_exponent == 0:
            factor0 = None
        else:
            factor0 = Pow(self.base.copy(), new_exponent, copy=False)
        # factor1 = "exponent*derivative(poly)"
        derivative_base = self.base.derive(index)
        factor1 = self.exponent * derivative_base
        if factor0 is None:
            summand1 = factor1
        else:
            summand1 = Product(factor0, factor1, copy=False)

        return Sum(summand0,summand1, copy=False)

    @doc(_Expression.docstring_of_replace)
    def replace(expression, index, value, remove=False):
        new_base = expression.base.replace(index,value,remove)
        new_exponent = expression.exponent.replace(index,value,remove)
        return Pow(new_base, new_exponent, copy=False)

class Log(_Expression):
    r'''
    The (natural) logarithm to base e (2.718281828459..).
    Store the expressions ``log(arg)``.

    :param arg:
        :class:`._Expression`;
        The argument of the logarithm.

    :param copy:
        bool;
        Whether or not to copy the `arg`.

    '''
    def __init__(self, arg, copy=True):
        self.number_of_variables = arg.number_of_variables
        self.arg = arg.copy() if copy else arg

    @cached_property
    def str(self):
        return 'log(' + str(self.arg) + ')'

    def __repr__(self):
        return self.str

    __str__ = __repr__

    def copy(self):
        "Return a copy of a :class:`.Log`."
        return Log(self.arg.copy(), copy=False)

    def simplify(self):
        'Apply ``log(1) = 0``.'
        if self.simplified:
            return self

        self.clear_cache()

        self.arg = self.arg.simplify()
        if type(self.arg) is Polynomial and len(self.arg.coeffs) == 1 and self.arg.coeffs[0] == 1 and (self.arg.expolist == 0).all():
            return Polynomial(np.zeros([1,len(self.arg.polysymbols)], dtype=int), np.array([0]), self.arg.polysymbols, copy=False)
        else:
            self.simplified = True
            return self

    @property
    def symbols(self):
        return self.arg.symbols

    def derive(self, index):
        '''
        Generate the derivative by the parameter indexed `index`.

        :param index:
            integer;
            The index of the paramater to derive by.

        '''
        # derivative(log(arg)) = 1/arg * derivative(arg) = arg**-1 * derivative(arg)
        symbols = self.arg.symbols
        minus_one = Polynomial(np.zeros([1,len(symbols)], dtype=int), np.array([-1]), symbols, copy=False)

        return Product(Pow(self.arg.copy(), minus_one, copy=False), self.arg.derive(index), copy=False)

    @doc(_Expression.docstring_of_replace)
    def replace(expression, index, value, remove=False):
        return Log( expression.arg.replace(index,value,remove) , copy=False )

def Expression(expression, polysymbols, follow_functions=False):
    '''
    Convert a sympy expression to an expression
    in terms of this module.

    :param expression:
        string or sympy expression;
        The expression to be converted

    :param polysymbols:
        iterable of strings or sympy symbols;
        The symbols to be stored as ``expolists``
        (see :class:`.Polynomial`) where
        possible.

    :param follow_functions:
        bool, optional (default = ``False``);
        If true, return the converted expression
        and a list of :class:`.Function`
        that occur in the `expression`.

    '''
    parsed_polysymbols = []
    for item in polysymbols:
        item = sympify_expression(item)
        assert item.is_Symbol, 'All `polysymbols` must be symbols'
        parsed_polysymbols.append(item)
    polysymbols = parsed_polysymbols

    sympy_function_calls = []
    functions = []

    def recursive_call(expression): # ``polysymbols`` from nonlocal context above
        try:
            return Polynomial.from_expression(expression, polysymbols)

        except:
            if expression.is_Mul:
                return Product(*(recursive_call(e) for e in expression.args))

            if expression.is_Pow:
                assert len(expression.args) == 2
                exponent = recursive_call(expression.args[1])
                try:
                    poly = Polynomial.from_expression(expression.args[0], polysymbols)
                    return ExponentiatedPolynomial(poly.expolist, poly.coeffs, exponent=exponent, polysymbols=polysymbols, copy=False)
                except:
                    return Pow(recursive_call(expression.args[0]), exponent)

            if expression.is_Add:
                return Sum(*(recursive_call(e) for e in expression.args))

            if isinstance(expression, sp.log):
                # make sure to have the natural log
                assert len(expression.args) == 1
                try:
                    return LogOfPolynomial.from_expression(expression.args[0], polysymbols)
                except:
                    return Log(recursive_call(expression.args[0]))

            if expression.is_Function:
                func = Function(expression.__class__.__name__, *(recursive_call(e) for e in expression.args))
                if follow_functions and expression not in sympy_function_calls:
                    sympy_function_calls.append(expression)
                    functions.append(func)
                return func

        raise ValueError('Could not parse the expression')

    parsed_expression = recursive_call(expression) if isinstance(expression, sp.Expr) else recursive_call( sympify_expression(str(expression)) )
    return (parsed_expression, functions) if follow_functions else parsed_expression

def refactorize(polyprod, *parameters):
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

    if not parameters:
        factorizable_powers = expolist_poly.min(axis=0)
        expolist_mono[:] += factorizable_powers
        expolist_poly[:] -= factorizable_powers
    else:
        for parameter in parameters:
            factorizable_power = expolist_poly[:,parameter].min()
            expolist_mono[:,parameter] += factorizable_power
            expolist_poly[:,parameter] -= factorizable_power

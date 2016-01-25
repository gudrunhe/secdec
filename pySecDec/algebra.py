"""Implementation of a simple computer algebra system"""

from .misc import argsort_2D_array
import numpy as np
import sympy as sp

class _Expression(object):
    '''
    Abstract base class for all expressions in this
    computer algebra system.

    '''
    def __add__(self, other):
        if not isinstance(other, _Expression):
            other = Polynomial([[0]*self.number_of_variables], [sp.sympify(other)], self.symbols)
        return Sum(self, other)
    __radd__ = __add__

    def __neg__(self):
        minus_one = Polynomial([[0]*self.number_of_variables], [-1], self.symbols)
        return Product(minus_one, self)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        if not isinstance(other, _Expression):
            other = Polynomial([[0]*self.number_of_variables], [sp.sympify(other)], self.symbols)
        return Product(self, other)

    def __pow__(self, other):
        if not isinstance(other, _Expression):
            other = Polynomial([[0]*self.number_of_variables], [sp.sympify(other)], self.symbols)
        return Pow(self, other)

    def __rpow__(self, other):
        if not isinstance(other, _Expression):
            other = Polynomial([[0]*self.number_of_variables], [sp.sympify(other)], self.symbols)
        return Pow(other, self)

class Function(_Expression):
    '''
    Symbolic function that can take care of
    parameter transformations.

    :param symbol:
        string;
        The symbol to be used to represent
        the `Function`.

    :param arguments:
        arbitrarily many :class:`._Expression`;
        The arguments of the `Function`.

    '''
    def __init__(self, symbol, *arguments):
        self.symbol = symbol
        self.number_of_arguments = len(arguments)
        self.number_of_variables = arguments[0].number_of_variables
        self.arguments = []
        for arg in arguments:
            assert arg.number_of_variables == self.number_of_variables, 'Must have the same number of variables in all arguments.'
            self.arguments.append(arg.copy())

    def __repr__(self):
        outstr_template = self.symbol + '(%s)'
        str_args = ','.join(str(arg) for arg in self.arguments)
        return outstr_template % str_args

    __str__ = __repr__

    def copy(self):
        "Return a copy of a :class:`.Function`."
        return Function(self.symbol, *self.arguments)

    def simplify(self):
        'Simplify the arguments.'
        self.arguments = [arg.simplify() for arg in self.arguments]
        return self

    @property
    def symbols(self):
        return self.arguments[0].symbols

    def derive(self, index):
        '''
        Generate the derivative by the parameter indexed `index`.
        The derivative of a function with `symbol` ``f`` by
        some `index` is denoted as ``df_d<index>``.

        :param index:
            integer;
            The index of the paramater to derive by.

        '''
        summands = []
        for argindex, arg in enumerate(self.arguments):
            summands.append(
                                Product(    # chain rule
                                            arg.derive(index),
                                            Function('d%s_d%i'%(self.symbol,argindex), *self.arguments)
                                       )
                           )
        return Sum(*summands).simplify()

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

    '''
    def __init__(self, expolist, coeffs, polysymbols='x'):
        self.expolist = np.array(expolist)
        assert len(self.expolist.shape) == 2, 'All entries in `expolist` must have the same length'
        if not np.issubdtype(self.expolist.dtype, np.integer):
            raise TypeError('All entries in `expolist` must be integer.')
        self.coeffs = np.array(coeffs)
        assert len(self.expolist) == len(self.coeffs), \
            '`expolist` (length %i) and `coeffs` (length %i) must have the same length.' %(len(self.expolist),len(self.coeffs))
        assert len(self.coeffs.shape) == 1, '`coeffs` must be one-dimensional'
        self.number_of_variables = self.expolist.shape[1]
        if not np.issubdtype(self.coeffs.dtype, np.number):
            parsed_coeffs = []
            for coeff in self.coeffs:
                if isinstance(coeff,_Expression):
                    assert coeff.number_of_variables == self.number_of_variables, 'Must have the same number of variables as the `Polynomial` for all coeffs'
                    parsed_coeffs.append(coeff.copy())
                else:
                    parsed_coeffs.append(sp.sympify(coeff))
            self.coeffs = np.array([coeff.copy() if isinstance(coeff,_Expression) else sp.sympify(coeff) for coeff in self.coeffs])
        if isinstance(polysymbols, str):
            self.polysymbols = sp.sympify([polysymbols + str(i) for i in range(self.number_of_variables)])
            for symbol in self.polysymbols:
                assert symbol.is_Symbol, 'All `polysymbols` must be symbols'
        else:
            self.polysymbols = []
            for item in polysymbols:
                if isinstance(item, sp.Symbol):
                    self.polysymbols.append(item)
                else:
                    item = sp.sympify(item)
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

        expression, polysymbols = sp.sympify((expression, polysymbols))

        for symbol in polysymbols:
            if not symbol.is_Symbol:
                raise TypeError("'%s' is not a symbol" % symbol)

        sympy_poly = sp.poly(expression, polysymbols)
        expolist = sympy_poly.monoms()
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
        return type(self)(self.expolist, self.coeffs, self.polysymbols)

    def has_constant_term(self):
        '''
        Return True if the polynomial can be written as:

        .. math::
            const + ...

        Otherwise, return False.

        '''
        return (self.expolist == 0).all(axis=1).any()

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
        summand1 = self.copy()
        summand1.expolist[:,index] -= 1
        summand1.coeffs *= self.expolist[:,index]
        summand1 = summand1.simplify()

        # summand2 = derivative(<coeff>) * x**k
        summand_2_coeffs = []
        need_summand2 = False
        for coeff in self.coeffs:
            if isinstance(coeff, _Expression):
                summand_2_coeffs.append(coeff.derive(index))
                need_summand2 = True
            else:
                summand_2_coeffs.append(0)

        if need_summand2:
            summand2 = Polynomial(self.expolist, summand_2_coeffs, self.polysymbols).simplify()
            return summand1 + summand2
        else:
            return summand1

    @property
    def symbols(self):
        return self.polysymbols

    def __add__(self, other):
        'addition operator'
        return self._sub_or_add(other, False)

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

            result = Polynomial(sum_expolist, sum_coeffs, self.polysymbols)
            result.simplify()
            return result

        elif np.issubdtype(type(other), np.number) or isinstance(other, sp.Expr):
            new_expolist = np.vstack([[0]*self.number_of_variables, self.expolist])
            new_coeffs = np.append(-other if sub else other, self.coeffs)
            outpoly = Polynomial(new_expolist, new_coeffs, self.polysymbols)
            outpoly.simplify()
            return outpoly

        else:
            return NotImplemented

    def __mul__(self, other):
        'multiplication operator'
        if  type(other) is Polynomial:
            assert self.number_of_variables == other.number_of_variables, 'Number of varibales must be equal for both factors in *'

            product_expolist = np.vstack([other.expolist + term for term in self.expolist])
            product_coeffs = np.hstack([other.coeffs * term for term in self.coeffs])

            result = Polynomial(product_expolist, product_coeffs, self.polysymbols)
            result.simplify()
            return result

        elif np.issubdtype(type(other), np.number) or isinstance(other, sp.Expr):
            new_coeffs = self.coeffs * other
            return Polynomial(self.expolist, new_coeffs, self.polysymbols)

        else:
            return NotImplemented

    __rmul__ = __mul__
    __radd__ = __add__

    def __neg__(self):
        'arithmetic negation "-self"'
        return Polynomial(self.expolist, [-coeff for coeff in self.coeffs], self.polysymbols)

    def __pow__(self, exponent):
        if not isinstance(exponent, int):
            return NotImplemented
        if exponent < 0:
            raise ValueError("The exponent must be nonnegative")
        if exponent == 0:
            return Polynomial([[0]*self.number_of_variables], [1], self.polysymbols)

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

    def simplify(self):
        '''
        Combine terms that have the same exponents of
        the variables.

        '''
        # Sort the expolist first, such that identical entries are
        # grouped together
        sort_key = argsort_2D_array(self.expolist)
        self.expolist = self.expolist[sort_key]
        self.coeffs = self.coeffs[sort_key]

        for i in range(1,len(self.coeffs)):
            # do not have to consider terms with zero coefficient
            if self.coeffs[i] == 0: continue

            previous_exponents = self.expolist[i-1]
            # search `self.expolist` for the same term
            # since `self.expolist` is sorted, must only compare with the previous term
            if (previous_exponents == self.expolist[i]).all():
                # add coefficients
                self.coeffs[i] += self.coeffs[i-1]
                # mark previous term for removal by setting coefficient to zero
                self.coeffs[i-1] = 0

        # remove terms with zero coefficient
        nonzero_coeffs = np.where(self.coeffs != 0)
        self.coeffs = self.coeffs[nonzero_coeffs]
        self.expolist = self.expolist[nonzero_coeffs]

        # need to have at least one term
        if len(self.coeffs) == 0:
            self.coeffs = np.array([0])
            self.expolist = np.array([[0]*self.number_of_variables])

        return self

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

    '''
    def __init__(self, expolist, coeffs, exponent=1, polysymbols='x'):
        Polynomial.__init__(self, expolist, coeffs, polysymbols)
        if np.issubdtype(type(exponent), np.number):
            self.exponent = exponent
        elif isinstance(exponent,_Expression):
            self.exponent = exponent.copy()
        else:
            self.exponent = sp.sympify(exponent)

    @staticmethod
    def from_expression(*args,**kwargs):
        raise NotImplementedError('Cannot create `ExponentiatedPolynomial` from expression.')

    __pow__ = _Expression.__pow__
    __radd__ = __add__ = _Expression.__add__
    __neg__ = _Expression.__neg__
    __sub__ = _Expression.__sub__
    __rsub__ = _Expression.__rsub__
    __rmul__ =__mul__ = _Expression.__mul__

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
            summand0_factors.append(LogOfPolynomial(self.expolist, self.coeffs, self.polysymbols))
            summand0 = Product(*summand0_factors).simplify()
        else:
            summand0 = None

        # summand1: poly**(exponent-1)*derivative(poly)
        # factor0 = poly**(exponent-1)   -->   simplification: (...)**0 = 1
        # do not need factor 0 in that case
        new_exponent = self.exponent - 1
        if new_exponent == 0:
            factor0 = None
        else:
            factor0 = ExponentiatedPolynomial(self.expolist,
                                              self.coeffs,
                                              new_exponent,
                                              self.polysymbols)
        # factor1 = "exponent*derivative(poly)"
        derivative_poly = Polynomial(self.expolist, self.coeffs, self.polysymbols).derive(index)
        factor1 = self.exponent * derivative_poly
        if factor0 is None:
            summand1 = factor1
        else:
            summand1 = Product(factor0, factor1)

        if summand0 is None:
            return summand1
        else:
            return Sum(summand0,summand1).simplify()

    def copy(self):
        "Return a copy of a :class:`.Polynomial` or a subclass."
        return type(self)(self.expolist, self.coeffs, self.exponent, self.polysymbols)

    def simplify(self):
        '''
        Apply the identity <something>**0 = 1 or
        <something>**1 = <somethng> if possible,
        otherwise call the simplify method of the base class.

        '''
        if isinstance(self.exponent, _Expression):
            self.exponent = self.exponent.simplify()
        if self.exponent == 0 or (isinstance(self.exponent, Polynomial) and (self.exponent.coeffs==0).all()):
            self.coeffs = np.array([1])
            self.expolist = np.array([[0]*self.number_of_variables])
            self.exponent = 1
        elif self.exponent == 1 or (isinstance(self.exponent, Polynomial) and len(self.exponent.coeffs)==1 and (self.exponent.coeffs==1).all() and (self.exponent.expolist==0).all()):
            return Polynomial(self.expolist, self.coeffs, self.polysymbols)
        else:
            super(ExponentiatedPolynomial, self).simplify()

        return self

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
        return LogOfPolynomial(poly.expolist, poly.coeffs, poly.polysymbols)

    def __repr__(self):
        return 'log(%s)' % Polynomial.__repr__(self)

    __str__ = __repr__

    __pow__ = _Expression.__pow__
    __radd__ = __add__ = _Expression.__add__
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
        #   --> factor0 = "poly**(-1)"
        factor0 = ExponentiatedPolynomial(self.expolist,
                                          self.coeffs,
                                          -1,
                                          self.polysymbols)

        #   --> factor1 = "derivative(poly)"
        factor1 = Polynomial(self.expolist, self.coeffs, self.polysymbols)
        factor1 = factor1.derive(index)

        return Product(factor0, factor1)

    def simplify(self):
        '''
        Apply the identity ``log(1) = 0``, otherwise call
        the simplify method of the base class.

        '''
        super(LogOfPolynomial, self).simplify()
        if len(self.coeffs) == 1 and self.coeffs[0] == 1 and (self.expolist == 0).all():
            return Polynomial([[0]*len(self.polysymbols)], [0], self.polysymbols)
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

    :math:`p_i` can be accessed with ``self.summands[i]``.

    Example:

    .. code-block:: python

        p = Sum(p0, p1)
        p0 = p.summands[0]
        p1 = p.summands[1]

    '''
    def __init__(self,*summands):
        self.summands = [summand.copy() for summand in summands]
        assert self.summands, 'Must have at least one summand'

        self.number_of_variables = self.summands[0].number_of_variables

        for summand in self.summands:
            if summand.number_of_variables != self.number_of_variables:
                raise TypeError('Must have the same number of variables for all summands.')

    def __repr__(self):
        stringified_summands = []
        for summand in self.summands:
            stringified_summands.append( '(' + str(summand) + ')' )
        return ' + '.join(stringified_summands)

    __str__ = __repr__

    def simplify(self):
        '''
        If one or more of ``self.summands`` is a
        :class:`Sum`, replace it by its summands.
        If only one summand is present, return that summand.
        Remove zero from sums.

        '''
        changed = True
        while changed:
            changed = False
            old_summands = self.summands
            self.summands = []
            for summand in old_summands:
                if isinstance(summand, Product):
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
            return zero
        else:
            return self

    @property
    def symbols(self):
        return self.summands[0].symbols

    def copy(self):
        "Return a copy of a :class:`.Sum`."
        return Sum(*self.summands)

    def derive(self, index):
        '''
        Generate the derivative by the parameter indexed `index`.

        :param index:
            integer;
            The index of the paramater to derive by.

        '''
        # derivative(p1 + p2 + ...) = derivative(p1) + derivative(p2) + ...
        return Sum(*(summand.derive(index) for summand in self.summands)).simplify()

class Product(_Expression):
    r'''
    Product of polynomials.
    Store one or polynomials :math:`p_i` to be interpreted as
    product :math:`\prod_i p_i`.

    :param factors:
        arbitrarily many instances of :class:`.Polynomial`;
        The factors :math:`p_i`.

    :math:`p_i` can be accessed with ``self.factors[i]``.

    Example:

    .. code-block:: python

        p = Product(p0, p1)
        p0 = p.factors[0]
        p1 = p.factors[1]


    '''
    def __init__(self,*factors):
        self.factors = [factor.copy() for factor in factors]
        assert self.factors, 'Must have at least one factor'

        self.number_of_variables = self.factors[0].number_of_variables

        for factor in self.factors:
            if factor.number_of_variables != self.number_of_variables:
                raise TypeError('Must have the same number of variables for all factors.')

    def __repr__(self):
        stringified_factors = []
        for factor in self.factors:
            stringified_factors.append( '(' + str(factor) + ')' )
        return ' * '.join(stringified_factors)

    __str__ = __repr__

    def copy(self):
        "Return a copy of a :class:`.Product`."
        return Product(*self.factors)

    def simplify(self):
        '''
        If one or more of ``self.factors`` is a
        :class:`.Product`, replace it by its factors.
        If only one factor is present, return that factor.
        Remove factors of one and zero.

        '''
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
                elif isinstance(factor, Polynomial):
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
            return one
        else:
            return self

    @property
    def symbols(self):
        return self.factors[0].symbols

    def derive(self, index):
        '''
        Generate the derivative by the parameter indexed `index`.

        :param index:
            integer;
            The index of the paramater to derive by.

        '''
        # product rule: derivative(p1 * p2 * ...) = derivative(p1) * p2 * ... + p1 * derivative(p2) * ...
        summands = []
        factors = list(self.factors) # copy
        for i,factor in enumerate(self.factors):
            factors[i] = factor.derive(index)
            summands.append(Product(*factors).simplify())
            factors[i] = factor
        return Sum(*summands).simplify()

class Pow(_Expression):
    r'''
    Exponential.
    Store two expressions ``A`` and ``B`` to be interpreted
    the exponential ``A**B``.

    :param base:
        :class:`._Expression`;
        The base ``A`` of the exponential.

    :param exponent:
        :class:`._Expression`;
        The exponent ``B``.

    '''
    def __init__(self, base, exponent):
        if base.number_of_variables != exponent.number_of_variables:
            raise TypeError('Must have the same number of variables for `base` and `exponent`.')

        self.number_of_variables = exponent.number_of_variables
        self.base = base.copy()
        self.exponent = exponent.copy()

    def __repr__(self):
        return '(' + str(self.base) + ') ** (' + str(self.exponent) + ')'

    __str__ = __repr__

    def copy(self):
        "Return a copy of a :class:`.Pow`."
        return Pow(self.base, self.exponent)

    @property
    def symbols(self):
        return self.base.symbols

    def simplify(self):
        '''
        Apply the identity <something>**0 = 1
        or <something>**1 = something if possible.

        '''
        self.base = self.base.simplify()
        self.exponent = self.exponent.simplify()

        if isinstance(self.exponent, Polynomial):
            if (self.exponent.coeffs==0).all():
                symbols = self.exponent.polysymbols
                return Polynomial([[0]*len(symbols)], [1], symbols)
            elif len(self.exponent.coeffs)==1 and (self.exponent.coeffs==1).all() and (self.exponent.expolist==0).all():
                return self.base
            else:
                return self
        else:
            return self

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
        summand0_factors.append(Log(self.base))
        summand0 = Product(*summand0_factors).simplify()

        # summand1: base**(exponent-1)*derivative(base)
        # factor0 = base**(exponent-1)   -->   simplification: (...)**0 = 1
        # do not need factor 0 in that case
        try:
            new_exponent = self.exponent - 1
        except:
            symbols = self.symbols
            new_exponent = Sum( self.exponent.copy(), Polynomial([[0]*len(symbols)],[-1],symbols) )
        if new_exponent == 0:
            factor0 = None
        else:
            factor0 = Pow(self.base.copy(), new_exponent)
        # factor1 = "exponent*derivative(poly)"
        derivative_base = self.base.derive(index)
        try:
            factor1 = self.exponent * derivative_base
        except:
            factor1 = Product(self.exponent.copy(), derivative_base)
        if factor0 is None:
            summand1 = factor1
        else:
            summand1 = Product(factor0, factor1)

        return Sum(summand0,summand1).simplify()

class Log(_Expression):
    r'''
    The (natural) logarithm to base e (2.718281828459..).
    Store the expressions ``log(arg)``.

    :param arg:
        :class:`._Expression`;
        The argument of the logarithm.

    '''
    def __init__(self, arg):
        self.number_of_variables = arg.number_of_variables
        self.arg = arg.copy()

    def __repr__(self):
        return 'log(' + str(self.arg) + ')'

    __str__ = __repr__

    def copy(self):
        "Return a copy of a :class:`.Log`."
        return Log(self.arg)

    def simplify(self):
        'Apply ``log(1) = 0``.'
        self.arg = self.arg.simplify()
        if isinstance(self.arg, Polynomial) and len(self.arg.coeffs) == 1 and self.arg.coeffs[0] == 1 and (self.arg.expolist == 0).all():
                return Polynomial([[0]*len(self.arg.polysymbols)], [0], self.arg.polysymbols)
        else:
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
        minus_one = Polynomial([[0]*len(symbols)], [-1], symbols)

        return Product(Pow(self.arg, minus_one), self.arg.derive(index)).simplify()

# TODO: extensively test this function
def make_expr(expression, polysymbols):
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

    '''
    parsed_polysymbols = []
    for item in polysymbols:
        item = sp.sympify(item)
        assert item.is_Symbol, 'All `polysymbols` must be symbols'
        parsed_polysymbols.append(item)
    polysymbols = parsed_polysymbols

    try:
        return Polynomial.from_expression(expression, polysymbols)

    except sp.PolynomialError:
        if expression.is_Mul:
            return Product(*(make_expr(e, polysymbols) for e in expression.args))

        if expression.is_Pow:
            assert len(expression.args) == 2
            exponent = make_expr(expression.args[1], polysymbols)
            try:
                poly = Polynomial.from_expression(expression.args[0], polysymbols)
                return ExponentiatedPolynomial(poly.expolist, poly.coeffs, exponent=exponent, polysymbols=polysymbols)
            except sp.PolynomialError:
                return Pow(make_expr(expression.args[0], polysymbols), exponent)

        if expression.is_Add:
            return Sum(*(make_expr(e, polysymbols) for e in expression.args))

        if isinstance(expression, sp.log):
            # make sure to have the natural log
            assert len(expression.args) == 1
            try:
                return LogOfPolynomial.from_expression(expression.args[0], polysymbols)
            except sp.PolynomialError:
                return Log(make_expr(expression.args[0], polysymbols))

        if expression.is_Function:
            return Function(expression.__class__.__name__, *(make_expr(e, polysymbols) for e in expression.args))

    raise ValueError('Could not parse the expression')

# TODO: replace should be an instancemethod for each subclass of `_Expression`
def replace(expression, index, value, remove=False):
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
    if isinstance(expression,Polynomial):
        outpoly = expression.copy()
        # act on coefficients first
        # nothing todo if the coeffs are just numbers
        if not np.issubdtype(outpoly.coeffs.dtype, np.number):
            for i,coeff in enumerate(outpoly.coeffs):
                try:
                    outpoly.coeffs[i] = replace(coeff, index, value, remove)
                except TypeError:
                    pass
        if isinstance(expression,ExponentiatedPolynomial):
            # replace in exponent if it has a type `replace` can handle
            try:
                outpoly.exponent = replace(outpoly.exponent, index, value, remove)
            except TypeError:
                pass
        if value != 1: # nothing to do if ``value==1`` since <coeff> * 1**<something> = <coeff>
            powers = expression.expolist[:,index]
            outpoly.coeffs *= np.array([value**int(power) for power in powers])
        if remove:
            outpoly.number_of_variables -= 1
            outpoly.expolist = np.delete(outpoly.expolist, index, axis=1)
            if index == -1:
                outpoly.polysymbols.pop()
            else:
                outpoly.polysymbols = outpoly.polysymbols[:index] + outpoly.polysymbols[index+1:]
        else:
            outpoly.expolist[:,index] = 0
        outpoly.simplify()
        return outpoly
    elif isinstance(expression, Product):
        outfactors = []
        for factor in expression.factors:
            outfactors.append(replace(factor,index,value,remove))
        return Product(*outfactors)
    elif isinstance(expression, Sum):
        outsummands = []
        for summand in expression.summands:
            outsummands.append(replace(summand,index,value,remove))
        return Sum(*outsummands)
    elif isinstance(expression, Pow):
        new_base = replace(expression.base,index,value,remove)
        new_exponent = replace(expression.exponent,index,value,remove)
        return Pow(new_base, new_exponent)
    elif isinstance(expression, Log):
        return Log( replace(expression.arg,index,value,remove) )
    elif isinstance(expression, Function):
        outfunction = expression.copy()
        outfunction.arguments = [replace(arg, index, value, remove) for arg in expression.arguments]
        if remove:
            outfunction.number_of_variables -= 1
        return outfunction
    else:
        raise TypeError('Can only operate on `Polynomial`, `Product`, `Sum`, `Pow`, `Log`, and `Function`, not `%s`' % type(expression))

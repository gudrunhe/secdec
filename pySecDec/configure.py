"""
Configure pySecDec:

The string output can be globally configured using
the functions in this module.

"""

# ********** default values **********
_default_powsymbol = _powsymbol = '**'
_default_coeffs_in_parentheses = _coeffs_in_parentheses = True
_default_polysymbol = _polysymbol = 'x'


# ***** getter/setter functions *****
def reset():
    'Reset to default configuration.'
    global _powsymbol; _powsymbol = _default_powsymbol
    global _coeffs_in_parentheses; _coeffs_in_parentheses = _default_coeffs_in_parentheses
    global _polysymbol; _polysymbol = _default_polysymbol

def polysymbol(symbol=None):
    '''
    Get or set the symbol used for the varibales
    of a :class:`pySecDec.polynomial.Polynomial`.
    Default: %s

    :param symbol:
        string, optional;
        If present, use this string to indicate the
        polynomial variables. If not present, return
        the string that is currently used.


    ''' % _default_polysymbol
    global _polysymbol
    if symbol is None:
        return _polysymbol
    else:
        assert issubclass(type(symbol),str), '`symbol` must be a string'
        _polysymbol = symbol

def powsymbol(symbol=None):
    '''
    Get or set the symbol used for "raise to power".
    Default: %s

    :param symbol:
        string, optional;
        If present, use this string to indicate powers
        in algebraic expressions. If not present, return
        the string that is currently used.

    ''' % _default_powsymbol
    global _powsymbol
    if symbol is None:
        return _powsymbol
    else:
        assert issubclass(type(symbol),str), '`symbol` must be a string'
        _powsymbol = symbol

def coeffs_in_parentheses(value=None):
    '''
    Get or set if the coefficients of a
    :class:`pySecDec.polynomial.Polynomial` should be decorated
    with parentheses.

    .. note::
        Setting this to true (default) is highly recommended.

    :param value:
        bool, optional;
        Whether to set parentheses or not.

    '''
    global _coeffs_in_parentheses
    if value is None:
        return _coeffs_in_parentheses
    else:
        assert issubclass(type(value),bool), '`value` must be a bool'
        _coeffs_in_parentheses = value

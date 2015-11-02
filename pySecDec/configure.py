"""
Configure pySecDec:

The string output can be globally configured using
the functions in this module.

"""

# ********** default values **********
_powsymbol = '**'
_coeffs_in_parentheses = True


# ***** getter/setter functions *****
def powsymbol(symbol=None):
    '''
    Get or set the symbol used for "raise to power".
    Default: **

    :param symbol:
        string, optional;
        If present, use this string to indicate powers
        in algebraic expressions. If not present, return
        the string that is currently used.

    '''
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

"""
Configure pySecDec:

The string output can be globally configured using
the functions in this module.

"""

# ********** default values **********
_powsymbol = '**'


# ***** getter/setter functions *****
def powsymbol(symbol=None):
    '''
    Get or set the symbol used for "raise to power".

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


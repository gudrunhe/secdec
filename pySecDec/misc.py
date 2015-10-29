"""miscellaneous routines"""

from itertools import chain,combinations

def powerset(iterable,exclude_empty=False):
    """
    Return an iterator over the powerset of a given set.
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)

    :param iterable:
        iterable;
        The set to generate the powerset for.

    :param exclude_empty:
        bool, optional;
        If True, skip the empty set in the powerset.
        Default is False.

    """
    # taken from python's own documentation
    s = list(iterable)
    powerset_iterator = iter(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))
    if exclude_empty:
        # The first element of the iterator is the empty set -> discard
        next(powerset_iterator)
    return powerset_iterator

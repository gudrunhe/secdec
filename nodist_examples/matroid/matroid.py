#! /usr/bin/env python
import numpy as np
from pySecDec.loop_integral import LoopIntegralFromGraph
from pySecDec.matrix_sort import Pak_sort
from pySecDec.algebra import Polynomial
from pySecDec.misc import flatten, argsort_ND_array

def canonicalize(internal_lines, external_lines, onshell_conditions, masses):
    '''
    Compute a canonical representation of the loop integral

    :param internal_lines:
        iterable of internal line specification, consisting
        of string or sympy expression for mass and a pair
        of strings or numbers for the vertices, e.g.
        [['m', [1,2]], ['0', [2,1]]].

    :param external_lines:
        iterable of external line specification, consisting
        of string or sympy expression for external momentum
        and a strings or number for the vertex, e.g.
        [['p1', 1], ['p2', 2]].

    :param onshell_conditions:
        list of strings;
        The value of each external momentum squared.

    :param masses:
        list of strings;
        The names of the masses.

    :return:
        dict;
        A canonical representation of the loop integral

    '''

    def drop_empty_cols(expolist, indices):
        '''
        Remove empty columns from `expolist` considering
        only those columns with index in `indices`.

        :param expolist:
            2D array of integers;
            The exponent list.

        :param indices:
            list of integers;
            The list of column indices to consider.

        :return:
            tuple (expolist, dropped_cols);
            Where `expolist` has zero columns dropped and
            `dropped_cols` is the indices of the dropped
            columns.

        '''
        indices = sorted(indices)
        dropped_cols = []
        for index in indices:
            shifted_index = index - len(dropped_cols)
            if not any(expolist[:, shifted_index]):
                expolist = np.hstack((expolist[:, :shifted_index], expolist[:, shifted_index + 1:]))
                dropped_cols.append(index)
        return expolist, dropped_cols

    def get_shortest_candidates(candidates, indices):
        '''
        Consider a list `candidates` of candidate expolists,
        drop zero columns present in `indices` and return
        a list of tuples containing the candidates with
        the most columns dropped and the indices of the
        dropped columns.

        :param candidates:
            list of tuples (expolist, dropped_cols);
            The candidates to consider and the columns
            already dropped.

        :param indices:
            list of integers;
            The columns to consider.

        :return:
            list of tuples (expolist, dropped_cols);
            The shortest candidates and the dropped columns.

        '''
        likely_candidates = []
        max_dropped_parameters = 0
        for candidate in candidates:
            expolist, dropped_cols = drop_empty_cols(candidate[0], indices)
            if len(dropped_cols) > max_dropped_parameters:
                likely_candidates = [(expolist, candidate[1] + dropped_cols)]
                max_dropped_parameters = len(dropped_cols)
            elif len(dropped_cols) == max_dropped_parameters:
                likely_candidates.append((expolist, candidate[1] + dropped_cols))
        return likely_candidates

    external_momenta = [ line[0] for line in external_lines]
    assert len(external_momenta)==len(onshell_conditions),\
           'Number of onshell conditions does not match the number of external momenta'

    # Generate replacement_rules
    replacement_rules = [('(' + mom + ')**2', '(' + onshell_condition + ')**2')
                         for mom, onshell_condition in zip(external_momenta, onshell_conditions)]

    # Generate N possible ways of imposing momentum conservation
    momentum_conservation = [(i, '-' + '-'.join(external_lines[j][0] for j in range(len(external_lines)) if i != j))
                             for i in range(len(external_lines))]

    # Generate a loop integral with no external_momenta eliminated
    li = LoopIntegralFromGraph(
        internal_lines=internal_lines,
        external_lines=external_lines,
        replacement_rules=replacement_rules
    )

    # Make li_f a polynomial also in the external_momenta and masses
    symbols = li.F.polysymbols + external_momenta + masses
    poly_coeffs = [Polynomial.from_expression(coeff, symbols) for coeff in li.F.coeffs]
    li_f = Polynomial(
        np.hstack(
            (li.F.expolist, np.zeros((len(li.F.expolist), len(external_momenta) + len(masses)), dtype=int))
        ), poly_coeffs, symbols)
    li_f = flatten(li_f).simplify()

    all_li_f = []
    if not momentum_conservation:

        # No external lines, so no momentum conservation to impose
        new_li_f = li_f.copy()
        new_li_f.coeffs = new_li_f.coeffs.astype(int)
        all_li_f.append(np.hstack((new_li_f.coeffs.reshape((-1, 1)), new_li_f.expolist)))

    else:

        # Impose momentum conservation on f polynomial in every possible way
        for momentum_rule in momentum_conservation:

            # Insert momentum conservation
            momentum_index = len(li.F.polysymbols) + momentum_rule[0]
            replacement = Polynomial.from_expression(momentum_rule[1], symbols)
            new_li_f = li_f.replace(momentum_index, replacement)

            new_li_f = flatten(new_li_f)
            new_li_f.coeffs = new_li_f.coeffs.astype(int)

            # Apply onshell conditions
            for index, repl in enumerate(onshell_conditions):
                offset_index = len(li.F.polysymbols) + index
                squared_mom = new_li_f.expolist[:, offset_index] == 2
                new_li_f.expolist[:, offset_index][squared_mom] = 0
                if repl == '0':
                    new_li_f.coeffs[squared_mom] = 0
                else:
                    mass_index = len(li.F.polysymbols) + len(external_momenta) + masses.index(repl)
                    new_li_f.expolist[:, mass_index][squared_mom] += 2

            new_li_f.simplify()

            all_li_f.append(np.hstack((new_li_f.coeffs.reshape((-1, 1)), new_li_f.expolist)))

    # Build u
    li_u = Polynomial(
        np.hstack((li.U.expolist, np.zeros((len(li.U.expolist), len(external_momenta) + len(masses)), dtype=int))),
        li.U.coeffs.copy(), symbols, copy=False).simplify()
    li_u.coeffs = li_u.coeffs.astype(int)
    li_u = np.hstack((li_u.coeffs.reshape((-1, 1)), li_u.expolist))
    labelled_li_u = np.hstack((np.zeros((len(li_u), 1), dtype=int), li_u))

    # Label u,f and prepend u to f for each candidate li_f
    all_li = []
    for li_f in all_li_f:
        labelled_li_f = np.hstack((np.ones((len(li_f), 1), dtype=int), li_f))
        li_uf = np.vstack((labelled_li_u, labelled_li_f))
        all_li.append(li_uf)

    # Create candidates
    candidates_li = [(candidate, []) for candidate in all_li]

    # Drop Feynman parameters missing in u and f polynomials
    index_of_first_feynman = 2  # Skip tag and coefficient of polynomials
    feynman_indices = range(index_of_first_feynman, index_of_first_feynman + len(li.F.polysymbols))
    candidates_li = get_shortest_candidates(candidates_li, feynman_indices)
    # Update feynman indices
    num_dropped = len(candidates_li[0][1])
    feynman_indices = range(index_of_first_feynman, index_of_first_feynman + len(li.F.polysymbols) - num_dropped)

    # Drop momenta which do not appear in f
    index_of_first_momentum = 2 + len(li.F.polysymbols) - len(candidates_li[0][1])
    momenta_indices = range(index_of_first_momentum, index_of_first_momentum + len(external_momenta))
    candidates_li = get_shortest_candidates(candidates_li, momenta_indices)
    # Update momenta indices
    num_dropped = len(candidates_li[0][1]) - num_dropped
    momenta_indices = range(index_of_first_momentum, index_of_first_momentum + len(external_momenta) - num_dropped)

    # Drop masses which do not appear in f
    index_of_first_mass = 2 + len(li.F.polysymbols) + len(external_momenta) - len(candidates_li[0][1])
    mass_indices = range(index_of_first_mass, index_of_first_mass + len(masses))
    candidates_li = get_shortest_candidates(candidates_li, mass_indices)
    # Update mass indices
    num_dropped = len(candidates_li[0][1]) - num_dropped
    mass_indices = range(index_of_first_mass, index_of_first_mass + len(masses) - num_dropped)

    # Pick polynomials with fewest terms
    len_shortest_poly = np.array([len(candidate[0]) for candidate in candidates_li]).min()
    candidates_li = [(candidate[0], candidate[1]) for candidate in candidates_li
                     if len(candidate[0]) == len_shortest_poly]

    # Remove duplicates
    all_li = np.array([candidate[0] for candidate in candidates_li])
    sort_key = argsort_ND_array(all_li)
    candidates_li = [candidates_li[key] for key in sort_key]
    unique_li = [candidates_li[0]]
    for candidate in candidates_li[1:]:
        if not np.array_equal(candidate[0], unique_li[-1][0]):
            unique_li.append(candidate)

    # Canonicalize u,f
    canonicals_li = []
    for li_uf in unique_li:
        Pak_sort(li_uf[0], feynman_indices, momenta_indices, mass_indices)
        canonicals_li.append(li_uf[0])

    sort_key = argsort_ND_array(canonicals_li)
    canonical_li = canonicals_li[sort_key[0]]

    expo_f = Polynomial.from_expression(li.exponent_F, ('eps',)).simplify()
    expo_u = Polynomial.from_expression(li.exponent_U, ('eps',)).simplify()

    # Output canonical representation
    representation = {}
    representation["uf"] = canonical_li.tolist()
    representation["expo_u_expolist"] = expo_u.expolist.tolist()
    representation["expo_u_coeffs"] = expo_u.coeffs.tolist()
    representation["expo_f_expolist"] = expo_f.expolist.tolist()
    representation["expo_f_coeffs"] = expo_f.coeffs.tolist()
    representation["num_feynman"] = len(feynman_indices)
    representation["num_momenta"] = len(momenta_indices)
    representation["num_mass"] = len(mass_indices)

    return representation

def dump_dict(dictionary, file):
    '''
    Write a dictionary to a file, sorting by keys.

    :param representation:
        dict;
        The dictionary to write to a file.

    :param file:
        str;
        The name of the file to write to.

    '''
    with open(file, 'w') as f:
        for key in sorted(dictionary.keys()):
            f.write(str(key) + ': ' + str(dictionary[key]) + '\n')

if __name__ == '__main__':

    # Matroid 1a
    internal_lines1a = [[0, [1, 2]], [0, [2, 3]], [0, [2, 3]], [0, [3, 4]], [0, [4, 1]], [0, [4, 1]]]
    external_lines1a = []
    onshell_conditions1a = []
    masses1a = []

    # Matroid 1b
    internal_lines1b = [[0, [1, 2]], [0, [2, 3]], [0, [2, 3]], [0, [3, 4]], [0, [3, 4]], [0, [4, 1]]]
    external_lines1b = []
    onshell_conditions1b = []
    masses1b = []

    # Matroid 2a
    internal_lines2a = [[0, [1, 2]], [0, [2, 3]], [0, [2, 3]], [0, [3, 4]], [0, [4, 1]], [0, [4, 1]]]
    external_lines2a = [['p1', 1], ['p2', 3]]
    onshell_conditions2a = ['0', '0']
    masses2a = []

    # Matroid 2b
    internal_lines2b = [[0, [1, 2]], [0, [2, 3]], [0, [2, 3]], [0, [3, 4]], [0, [3, 4]], [0, [4, 1]]]
    external_lines2b = [['p1', 1], ['p2', 3]]
    onshell_conditions2b = ['0', '0']
    masses2b = []

    # Matroid 3a
    internal_lines3a = [['m1', [1, 2]], ['m1', [2, 3]], ['m1', [2, 3]], ['m1', [3, 4]], ['m1', [4, 1]], ['m1', [4, 1]]]
    external_lines3a = []
    onshell_conditions3a = []
    masses3a = ['m1']

    # Matroid 3b
    internal_lines3b = [['m1', [1, 2]], ['m1', [2, 3]], ['m1', [2, 3]], ['m1', [3, 4]], ['m1', [3, 4]], ['m1', [4, 1]]]
    external_lines3b = []
    onshell_conditions3b = []
    masses3b = ['m1']

    # Matroid 4a
    internal_lines4a = [[0, [1, 2]], [0, [1, 4]], [0, [1, 4]], [0, [4, 3]], [0, [3, 2]], [0, [3, 2]]]
    external_lines4a = [['p1', 1], ['p2', 2], ['p3', 4]]
    onshell_conditions4a = ['m1', 'm2', 'm3']
    masses4a = ['m1', 'm2', 'm3']

    # Matroid 4b
    internal_lines4b = [[0, [1, 2]], [0, [1, 4]], [0, [1, 4]], [0, [4, 3]], [0, [4, 3]], [0, [3, 2]]]
    external_lines4b = [['p1', 1], ['p2', 2], ['p3', 4]]
    onshell_conditions4b = ['m1', 'm2', 'm3']
    masses4b = ['m1', 'm2', 'm3']

    # Matroid 5a
    internal_lines5a = [[0, [1, 5]], [0, [5, 2]], [0, [1, 4]], [0, [1, 4]], [0, [4, 3]], [0, [3, 2]], [0, [3, 2]]]
    external_lines5a = [['p1', 1], ['p2', 2], ['p3', 4], ['p4', 5]]
    onshell_conditions5a = ['m1', '0', '0', '0']
    masses5a = ['m1']

    # Matroid 5b
    internal_lines5b = [[0, [1, 5]], [0, [5, 2]], [0, [1, 4]], [0, [1, 4]], [0, [4, 3]], [0, [4, 3]], [0, [3, 2]]]
    external_lines5b = [['p1', 1], ['p2', 2], ['p3', 4], ['p4', 5]]
    onshell_conditions5b = ['m1', '0', '0', '0']
    masses5b = ['m1']

    call1a = canonicalize(internal_lines1a, external_lines1a, onshell_conditions1a, masses1a)
    call1b = canonicalize(internal_lines1b, external_lines1b, onshell_conditions1b, masses1b)

    call2a = canonicalize(internal_lines2a, external_lines2a, onshell_conditions2a, masses2a)
    call2b = canonicalize(internal_lines2b, external_lines2b, onshell_conditions2b, masses2b)

    call3a = canonicalize(internal_lines3a, external_lines3a, onshell_conditions3a, masses3a)
    call3b = canonicalize(internal_lines3b, external_lines3b, onshell_conditions3b, masses3b)

    call4a = canonicalize(internal_lines4a, external_lines4a, onshell_conditions4a, masses4a)
    call4b = canonicalize(internal_lines4b, external_lines4b, onshell_conditions4b, masses4b)

    call5a = canonicalize(internal_lines5a, external_lines5a, onshell_conditions5a, masses5a)
    call5b = canonicalize(internal_lines5b, external_lines5b, onshell_conditions5b, masses5b)

    dump_dict(call1a,'call1a.txt')
    dump_dict(call1b,'call1b.txt')

    dump_dict(call2a,'call2a.txt')
    dump_dict(call2b,'call2b.txt')

    dump_dict(call3a,'call3a.txt')
    dump_dict(call3b,'call3b.txt')

    dump_dict(call4a,'call4a.txt')
    dump_dict(call4b,'call4b.txt')

    dump_dict(call5a,'call5a.txt')
    dump_dict(call5b,'call5b.txt')

    print(call1a == call1b)
    print(call2a == call2b)
    print(call3a == call3b)
    print(call4a == call4b)
    print(call5a == call5b)

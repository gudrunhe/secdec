"""Routines to Feynman parametrize a loop integral given the graph."""

from .common import LoopIntegral
from ..algebra import Polynomial, ExponentiatedPolynomial
from ..misc import missing, cached_property, sympify_symbols, assert_degree_at_most_max_degree, sympify_expression, rec_subs
from itertools import combinations
import sympy as sp
import numpy as np

class LoopIntegralFromGraph(LoopIntegral):
    __doc__ = '''
    Construct the Feynman parametrization of a
    loop integral from the graph using the cut construction method.

    Example:

    >>> from pySecDec.loop_integral import *
    >>> internal_lines = [['0',[1,2]], ['m',[2,3]], ['m',[3,1]]]
    >>> external_lines = [['p1',1],['p2',2],['-p12',3]]
    >>> li = LoopIntegralFromGraph(internal_lines, external_lines)
    >>> li.exponentiated_U
    ( + (1)*x0 + (1)*x1 + (1)*x2)**(2*eps - 1)
    >>> li.exponentiated_F
    ( + (m**2)*x2**2 + (2*m**2 - p12**2)*x1*x2 + (m**2)*x1**2 + (m**2 - p1**2)*x0*x2 + (m**2 - p2**2)*x0*x1)**(-eps - 1)

    :param internal_lines:
        iterable of internal line specification, consisting
        of string or sympy expression for mass and a pair
        of strings or numbers for the vertices, e.g.
        ``[['m', [1,2]], ['0', [2,1]]]``.

    :param external_lines:
        iterable of external line specification, consisting
        of string or sympy expression for external momentum
        and a strings or number for the vertex, e.g.
        ``[['p1', 1], ['p2', 2]]``.

    ''' + LoopIntegral.common_properties_doc

    def __init__(self, internal_lines, external_lines, replacement_rules=[], Feynman_parameters='x', \
                 regulators=None, regulator=None, dimensionality='4-2*eps', powerlist=[]):

        # sympify and store internal lines
        assert len(internal_lines) > 0, \
            "To define a loop integral please input a graph with at least one internal line."
        self.internal_lines=[]
        for line in internal_lines:
            assert len(line)==2 and len(line[1])==2, \
                "Internal lines must have the form [mass, [vertex, vertex]]."
            mass = sympify_symbols([line[0]], "Names of internal masses must be symbols or numbers.", \
                                   allow_number=True)[0]
            vertices = sympify_symbols(line[1], "Names of vertices must be symbols or numbers.", \
                                       allow_number=True)
            self.internal_lines.append([mass,vertices])
        self.P = len(self.internal_lines)

        # sympify and store external lines and create list of external momenta
        self.external_lines=[]
        self.external_momenta=[]
        for line in external_lines:
            assert len(line)==2, "External lines must have the form [momentum, vertex]."
            extmom = sympify_expression(line[0])
            vertex = sympify_symbols([line[1]], "Names of vertices must be symbols or numbers.", \
                                     allow_number=True)[0]
            self.external_lines.append([extmom,vertex])
            if extmom.is_Symbol and extmom not in self.external_momenta:
                self.external_momenta.append(extmom)

        for line in self.external_lines:
            assert_degree_at_most_max_degree(line[0], self.external_momenta, 1, "The first element of the external line specifications must be at most linear combinations of external momenta.")

        # store properties shared between derived classes
        regulators = \
                regulators if regulators is not None else \
                [regulator] if regulator is not None else \
                ["eps"]
        self.all_momenta = self.external_momenta
        self.set_common_properties(replacement_rules, Feynman_parameters, regulators,
                                   dimensionality, powerlist)

        # remove `internal_lines` and `Feynman_parameters` that are set to zero by the `powerlist`
        for i in range(self.P-1,-1,-1): # traverse backwards to stay consistent with the indexing
            if self.powerlist[i] == 0:
                vertex1, vertex2 = self.internal_lines[i][1]
                self.internal_lines = self.internal_lines[:i] + self.internal_lines[i+1:]
                self.P -= 1
                self.powerlist = self.powerlist[:i] + self.powerlist[i+1:]
                self.derivativelist = self.derivativelist[:i] + self.derivativelist[i+1:]
                self.Feynman_parameters = self.Feynman_parameters[:i] + self.Feynman_parameters[i+1:]

                # re-connect graph after removing line -> pinch
                for line in self.internal_lines:
                    for j in range(2):
                        if line[1][j] == vertex1:
                            line[1][j] = vertex2
                for line in self.external_lines:
                    if line[1] == vertex1:
                        line[1] = vertex2

        # calculate number of loops from the relation #loops = #internal lines - (#vertices - 1)
        self.V = len(self.intverts)
        self.L = self.P - (self.V - 1)
        assert self.L > 0, \
            "To define a loop integral please input a graph with at least one closed loop."

        # no support for tensor integrals in combination with cutconstruct for now
        self.highest_rank = 0
        self.preliminary_numerator = Polynomial(np.zeros([1,len(self.Feynman_parameters)+2], dtype=int), \
                                                np.array([1]), self.Feynman_parameters+sympify_expression(['U','F']), \
                                                copy=False)

    @cached_property
    def intverts(self):
        # creates a list of all internal vertices and indexes them
        # returns a dictionary that relates the name of a vertex to its index
        lines = self.internal_lines
        vertices=set()
        for line in lines:
            vertices = vertices | set(line[1])
        vertices=list(vertices)
        vertnames = {}

        for i in range(len(vertices)):
            vertnames[vertices[i]] = i
        return vertnames

    @cached_property
    def vertmatrix(self):  # create transition matrix representation of underlying graph

        # each vertex is trivially connected to itself, so start from unit matrix:
        numvert = self.V+len(self.external_lines)
        M = np.identity(numvert, dtype=int)

        # for each internal propagator connecting two vertices add an entry in the matrix
        for i in range(self.P):
            start = self.intverts[self.internal_lines[i][1][0]]
            end = self.intverts[self.internal_lines[i][1][1]]
            M[start,end] += 1
            M[end,start] += 1

        # for each external line add a vertex and an entry in the matrix
        for i in range(len(self.external_lines)):
            start = self.V + i
            end = self.intverts[self.external_lines[i][1]]
            M[start,end] += 1
            M[end,start] += 1

        return M

    @cached_property
    def preliminary_U(self):

        expolists=[]
        coeffs=[]

        # iterate over all possible L-fold cuts
        for cut in combinations(range(self.P), self.L):
            # find uncut propagators
            uncut = missing(range(self.P),cut)

            # Define transition matrix for cut graph by removing cut propagators
            newmatrix = np.matrix(self.vertmatrix) # copies data by default
            for i in cut:
                start = self.intverts[self.internal_lines[i][1][0]]
                end = self.intverts[self.internal_lines[i][1][1]]
                newmatrix[start,end] -= 1
                newmatrix[end,start] -= 1

            # Check if cut graph is connected
            numvert = self.V + len(self.external_lines)
            if(0 in newmatrix**numvert):
                # not connected if exponentiated matrix has a zero
                continue

            # construct monomial of Feynman parameters of cut propagators to be added to U
            expolist=[0]*self.P
            for i in cut:
                expolist[i]=1
            expolists.append(expolist)
            coeffs.append(1)

        return Polynomial(expolists, coeffs, polysymbols=self.Feynman_parameters)

    @cached_property
    def preliminary_F(self):

        expolists=[]
        coeffs=[]

        # iterate over all possible (L+1)-fold cuts
        for cut in combinations(range(self.P), self.L+1):
            # find uncut propagators
            uncut = missing(range(self.P),cut)

            # Define transition matrix for cut graph by removing cut propagators
            newmatrix = np.matrix(self.vertmatrix) # copies data by default
            for i in cut:
                start = self.intverts[self.internal_lines[i][1][0]]
                end = self.intverts[self.internal_lines[i][1][1]]
                newmatrix[start,end] -= 1
                newmatrix[end,start] -= 1

            numvert = self.V + len(self.external_lines)
            newmatrix = newmatrix**numvert

            # find all internal vertices *not* connected to vertex 0 (arbitrary choice)
            intnotconnectedto0 = []
            for i in range(self.V):
                if newmatrix[0,i]==0:
                    intnotconnectedto0.append(i)

            # find all external vertices *not* connected to vertex 0
            extnotconnectedto0 = []
            for i in range(self.V,numvert):
                if newmatrix[0,i]==0:
                    extnotconnectedto0.append(i)

            # check if all vertices not connected to 0 are connected to each other
            valid2tree = True
            notconnectedto0 = intnotconnectedto0 + extnotconnectedto0
            for i,j in combinations(notconnectedto0, 2):
                if newmatrix[i,j]==0:
                    # there are more than two disconnected components -> not a valid two-tree cut
                    valid2tree = False
                    break
            # drop current cut if it is not a valid two-tree cut
            if not valid2tree:
                continue

            # find momementa running through the two cut lines
            # choose either all the external momenta connected to vertex 0 or the complement
            cutmomenta = []
            if(len(extnotconnectedto0) <= len(self.external_lines)-len(extnotconnectedto0)):
                cutmomenta = extnotconnectedto0
            else:
                cutmomenta = missing(range(self.V,numvert),extnotconnectedto0)

            # construct monomial of Feynman parameters of cut propagators to be added to F,
            # if the momentum flow through the cuts is non-zero
            if cutmomenta:
                expolist=[0]*self.P
                for i in cut:
                    expolist[i]=1
                # sum the momenta flowing through the cuts, square it, and use replacement rules
                sumsqr = sum(self.external_lines[i-self.V][0] for i in cutmomenta)**2
                sumsqr = sumsqr.expand().subs(self.replacement_rules)
                expolists.append(expolist)
                coeffs.append(-sumsqr)

        if expolists:
            F0 = Polynomial(expolists, coeffs, polysymbols=self.Feynman_parameters)
        else:
            F0 = 0

        # construct terms proportial to the squared masses
        expolists=[]
        coeffs=[]
        for i in range(self.P):
            coeff = self.internal_lines[i][0]
            if coeff != 0:
                coeff = rec_subs((coeff**2), self.replacement_rules)
                coeffs.append(coeff)
                expolist = [0]*self.P
                expolist[i] = 1
                expolists.append(expolist)

        if expolists:
            Fm = Polynomial(expolists, coeffs, polysymbols=self.Feynman_parameters)
        else:
            Fm = 0

        return F0 + self.preliminary_U*Fm

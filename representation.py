from sage.rings.integer import Integer
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.rational_field import QQ
from copy import deepcopy

class gl:
    def __init__(self,N,asymptotic=False,base=QQ):
        self.N = N
        if asymptotic:
            self.base_ring = LaurentPolynomialRing(base,'h')
            self.hbar = self.base_ring.gen()
        if not asymptotic:
            self.base_ring = base
            self.hbar = Integer(1)

    def _lie_names(self):
        return ['E_' + str(i+1) + '_' + str(j+1) for i in range(self.N) for j in range(self.N)]
    
    def _lie_bracket(self):
        br = {}
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    for l in range(self.N):
                        if not i==j==k==l:
                            br[('E_' + str(i+1) + '_' + str(j+1),'E_' + str(k+1) + '_' + str(l+1))] = {'E_' + str(i+1) + '_' + str(l+1) : self.hbar*Integer(j==k), 'E_' + str(k+1) + '_' + str(j+1) : -self.hbar*Integer(l==i)}
        return br

    

class TensorPowerVectorRepresentation:
    def __init__(self,N,asymptotic,prefixes,base):
        self.gl = gl(N,asymptotic,base)
        self.base_ring = self.gl.base_ring
        self.N = self.gl.N
        self.hbar = self.gl.hbar
        self.prefixes = prefixes

    def _vector_names(self):
        ls = []
        for prefix in self.prefixes:
            ls += [prefix + '_' + str(k+1) for k in range(self.N)]
        return ls

    def _vector_rep_bracket(self,prefixes):
        br = {}
        for prefix in prefixes:
            for i in range(self.N):
                for j in range(self.N):
                    for k in range(self.N):
                        br[('E_' + str(i+1) + '_' + str(j+1),prefix + '_' + str(k+1))] = {prefix + '_' + str(i+1) : self.hbar*Integer(j==k)}
        return br
    
    def _all_names(self):
        return self.gl._lie_names() + self.vector_names(self.prefixes)
    
    def _all_brackets(self):
        d = deepcopy(self._vector_rep_bracket(self.prefixes))
        d.update(self.gl._lie_bracket())
        return d
    
    def _triangular_basis(self):
        d = []
        for diff in range(1,self.N):
            for i in range(1+diff,self.N+1):
                d.append('E_' + str(i) + '_' + str(i-diff))
        for i in range(1,self.N+1):
            d.append('E_' + str(i) + '_' + str(i))
        for diff in range(1,self.N):
            for i in range(1+diff,self.N+1):
                d.append('E_' + str(i-diff) + '_' + str(i))
        return d + self._vector_names()


from sage.rings.integer import Integer
from sage.algebras.lie_algebras.lie_algebra import LieAlgebra
from sage.arith.all import factorial
from numpy import sum
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import matrix

from representation import TensorPowerVectorRepresentation

class Whittaker:
    
    """
    This is the class encoding the Harish-Chandra bimodule U(gl_N) \otimes (Sym(C^N))^{\otimes m} and computing
    the Whittaker vector with the help of the Kirillov projector.

    Parameters
    ----------
    N : int
        Dimension
    asymptotic : bool
        Determines if constructions involve formal parameter \hbar. If True, then self.hbar is an element of the Laurent
        polynomial ring; if False, then self.hbar=1.
    prefixes : list
        Variable names in the copies of Sym(C^N). 
    base : Field (Sage)
        Base field, by default - rational numbers.

    Attributes 
    ----------
    t : PolynomialRing_field.element_class
        Formal variable primarly used for quantum minors.
    hbar : LaurentPolynomial_univariate or Integer
        Formal asymptotic variable (if asymptotic=False, then =1)
    G : FiniteFamily
        The algebra generators of the universal enveloping algebra corresponding to the matrix unit basis of gl_N.
        Example: G['E_1_2']
    """


    def __init__(self,N,asymptotic=False,prefixes=['u','v'],base=QQ):
        self.N = N
        self.asymptotic = asymptotic

        representation = TensorPowerVectorRepresentation(self.N,self.asymptotic,prefixes,base)
        self._base_ring = representation.base_ring
        _base = PolynomialRing(representation.base_ring,'t')

        self.t = _base.gen()
        self.hbar = representation.hbar

        self._pbw = LieAlgebra(representation.base_ring,
                              representation._all_brackets(),
                              representation._triangular_basis()).pbw_basis(prefix='G')
        self.G = self._pbw.algebra_generators()

        self._nilpotent_dict = self._regular_nil_dict(self.N)
        self._capelli_matrix = self._mat()


    def noe(self,i,j,el,shift=0):
        """
            Computes el*P_{ij}^{\psi}(u).

            Parameters
            ----------
            i,j : int
                Correspond to i,j above.
            el : PoincareBirkhoffWittBasis.element_class
                Element on which we act by this power series.
            shift : int
                Corresponds to parameter u.
        """

        indices = list(range(j,i))
        A = Integer(-1)**(i+j)*self.quotient(self.minor(indices,indices))
        sum = Integer(0)
        aux_el = el
        K = Integer(1)
        k = 0
        while aux_el != Integer(0):
            sum += self.hbar**(-k)*(1/factorial(k))*aux_el*K
            k += 1
            aux_el = self.quotient(aux_el*(self._nilpotent_dict['E_' + str(i) + '_' + str(j)]-self.G['E_' + str(i) + '_' + str(j)]))
            K *= A(-self.hbar*(Integer(k-1)-shift))
        return self.quotient(sum)
    
    def proj(self,el,truncation=0,**kwargs):
        """
            Computes the action of the Kirillov projector P_{\m_N}^{\psi}(u_1,\ldots,u_{N-1}).

            Parameters
            ----------
            el : PoincareBirkhoffWittBasis.element_class
            truncation : int
                Corresponds to `i` in the partial projection along gl_{N-i} embedded as lower right corner.
            **kwargs:
                u_1,u_2,... : int
                Correspond to same-named variables in P_{\m_N}^{\psi}(u_1,\ldots,u_{N-1}). By default, all set to zero.
        """

        new_el = el
        for j in reversed(range(1+truncation,self.N)):
            if kwargs.get('u_' + str(j+1)) != None:
                s = kwargs['u_' + str(j+1)]
            else:
                s = 0
            for i in range(j+1,self.N+1):
                new_el = self.noe(i,j,new_el,shift=s)
        return new_el
    
    def vector(self,prefix):
        """
            Returns a list with the generators of the symmetric algebra corresponding to `prefix`.
        """
        return [Integer(0)] + [self.G[prefix + '_' + str(i+1)] for i in range(self.N)]
    
    
    def quotient(self,expr):
        """
            Computes the image of `expr` under the quotient map corresponding to the non-degenerate Whittaker character
            \psi(E_{ij}) = \delta_{i,j+1}.
        """
        new_element = Integer(0)
        for i,c in (self.t**0*expr).dict().items():
            if c == self._pbw.one():
                new_element += self.t**i
            else:
                for monom,k in c._monomial_coefficients.items():
                    new_monom = k
                    for key,deg in monom._sorted_items():
                        if self._nilpotent_dict.get(key) != None:
                            new_multiplier = self._nilpotent_dict[key]**deg
                            new_monom *= new_multiplier
                        else:
                            new_multiplier = self.G[key]**deg
                            new_monom *= new_multiplier
                    new_element += new_monom*self.t**i
        return new_element
        
    def isInv(self,el):
        """
            Computes if `el` is Whittaker.
        """
        val = 1
        if el != Integer(0):
            for n in self._nilpotent_dict.keys():
                psi = self._nilpotent_dict[n]
                commutator = self.quotient(el*self.G[n] - psi*el) 
                val *= (commutator == Integer(0))
                if not commutator == Integer(0):
                    print(n + ': ',commutator)
        return val == 1
    
    def minor(self,rows,cols,shift=0):
        """
            Computes the quantum minor with rows=`rows` and columns=`cols`. 
            Parameter `shift` corresponds to a shift of variable t -> t+shift.
        """
        if len(rows) != len(cols):
            print("Number of rows should be equal to the number of colums")
        else:
            rows = [a-1 for a in rows]
            cols = [a-1 for a in cols]
            pmts_l = self._permutations(rows)
            pmts_l = list(zip(pmts_l[0],pmts_l[1]))

            sum = Integer(0)

            for (pmt,sgn) in pmts_l:
                prod = Integer(1)
                for i,ind in enumerate(pmt):
                    prod = prod*(self._capelli_matrix[ind,cols[i]])(self.t - self.hbar*(Integer(i)-shift))
                sum = sum + Integer(sgn)*prod

            return sum

    
    def _regular_nil_dict(self,N):
        d = {}
        for j in range(N-1):
            for i in range(j+1,N):
                d['E_' + str(i+1) + '_' + str(j+1)] = Integer(i==j+1)
        return d 
    
    
    def _mat(self):
        A = [[Integer(0)]*(self.N) for i in range(self.N)]
        for i in range(self.N):
            for j in range(self.N):
                if j == i:
                    A[i][j] = self.G['E_' + str(i+1) + '_' + str(j+1)] + self.t
                else:
                    A[i][j] = self.G['E_' + str(i+1) + '_' + str(j+1)]
        return matrix(A)
    
    def _permutations(self,xs):
        N = len(xs)
        if N == 0 or N == 1:
            return ([xs],[1])
        else:
            final = []
            final_sgns = []
            for i in range(N):
                a = xs[i]
                loc_xs = xs[:i] + xs[i+1:]
                sub_permutations, sub_signs = self._permutations(loc_xs) 
                final = final + [[a] + q for q in sub_permutations]
                final_sgns = final_sgns + [(-1)**i*s for s in sub_signs]
            return (final,final_sgns)
# reference:
# H.M.Edwards "Fermat's Last Theorem" section 4.5

from CYC import *
from sympy.abc import x

class Period:
    """ cyclotomic periods
    represented by array of int of length e
    """
    @staticmethod
    def init(f):
        """ set length f of period
        Assume CYC.init(n) has been executed.
        Assume f divides n-1
        """
        n1 = CYC.n1
        e,r = divmod(n1,f) # e = number of conjugates
        if r:
            raise RuntimeError("wrong f in Period")
        r = range(e)
        i = [CYC.g_pow[np.r_[i:n1:e]] for i in r]
        j = [np.roll(np.r_[:e], -i) for i in r]
        Period.e, Period.f = e,f
        Period.index = np.asarray(i)
        Period.cjind = np.asarray(j)

        # multiplication table
        w = np.zeros((e,e), dtype='object')
        for i in r:
            for j in range(f):
                k = (1 + Period.index[i,j])%CYC.n
                if k: w[i, CYC.log_g[k]%e] += 1
                else: w[i] -= f

        w = [np.roll(w,i,axis=(0,1)) for i in r]
        Period.w = np.asarray(w, dtype='object')
        # eta_i*eta_j = sum_{k=0}^{e-1} w[i,j,k]*eta_k

        # minimal polynomial
        u,a,y = [], 1, -Period.basis(0)
        for i in r:
            u.append(y)
            for j in range(i,-1,-1):
                if j<i: u[j] *= y
                if j: u[j] += u[j-1]
            y = y.conj()

        for c in u[::-1]: a = a*x + c.to_int()
        Period.MPoly = sp.expand(a)
        # (x-eta_0)(x-eta_1)...(x-eta_{e-1})

    @staticmethod
    def basis(i=0):
        """ return eta_i (i=0,...,e-1) """
        return Period(np.roll([1]+[0]*(Period.e-1), i))

    def __init__(a,c):
        """ set coefficients c of periods
        as sum_{j=0}^{e-1} c_j eta_j
        """
        a.c = np.asarray(c, dtype='object')
        a.c.resize(a.e)

    def __neg__(a): return Period(-a.c) # -a

    def __add__(a,b): # a+b, b must be Period object
        return Period(a.c + b.c)

    def __mul__(a,b): # a*b, b must be Period object
        c = np.dot(a.c, np.dot(b.c, b.w))
        return Period(c)

    # conjugate by cyclic permutation: eta_j -> eta_{j+i}
    def conj(a, i=1): return Period(a.c[a.cjind[i]])

    # product of all conjugates of a
    def norm(a):
        b = [a.conj(i) for i in range(a.e)]
        return np.prod(b, axis=0).to_int()

    # -(coeff of eta_0)
    def to_int(a): return -a.c[0]

    def isRational(a):
        """ test if a is rational.
        return None if a is not rational,
        return a.to_int() if a is rational.
        """
        if np.all(a.c == a.c[0]):
            return a.to_int()
        else: return None

    def to_CYC(a):
        """ transform a to cyclotomic integer
        as a polynomial in omega = exp(2pi*i/n)
        """
        y = np.zeros(CYC.n, dtype='object')
        y[a.index] = a.c[:,np.newaxis]
        return CYC(y - y[-1]) # normalize

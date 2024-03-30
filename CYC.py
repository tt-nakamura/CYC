# reference:
# H.M.Edwards "Fermat's Last Theorem" section 4.2

import numpy as np
import sympy as sp
from math import gcd

class CYC:
    """ cyclotomic integers """
    @staticmethod
    def init(n):
        """ set degree n of cyclotomic integer
        assume n is odd prime
        """
        if n<3 or not sp.isprime(n):
            raise RuntimeError('wrong n')
        n1,g = n-1, sp.primitive_root(n)
        CYC.n, CYC.n1 = n,n1
        CYC.g_pow = np.empty(n1, dtype=np.int)
        CYC.log_g = np.empty(n, dtype=np.int)
        CYC.cjind = np.empty((n1,n), dtype=np.int)
        CYC.order = np.empty(n, dtype=np.int)
        a = 1; CYC.order[0] = 1
        for i in range(n1):
            CYC.g_pow[i] = a
            CYC.log_g[a] = i
            CYC.order[a] = n1//gcd(i,n1)
            CYC.cjind[-i] = np.r_[:a*n:a] % n
            a = a * g % n

    @staticmethod
    def ord(a): # order of a mod n; ord(0)==1
        return CYC.order[a%CYC.n]

    def __init__(a, c=0):
        """ set coefficients c of cyclotomic integer
           as a polynomial of degree n-1.
        c: 1d-array or scalar of int
           in order of increasing power of omega
           where omega = exp(2*pi*i/n)
        """
        a.c = np.asarray(c, dtype='object')
        a.c.resize(a.n)

    def __hash__(a): return hash(sum(a.c[::2]))
    def __repr__(a): return str(np.trim_zeros(a.c, 'b'))
    def __neg__(a): return CYC(-a.c) # -a

    # make degree of a less than n-1
    def normalize(a): return CYC(a.c - a.c[-1])

    # conjugate by replacing omega with omega^{g^i}
    def conj(a, i=1): return CYC(a.c[a.cjind[i]])

    # product of all conjugates of a
    def norm(a):
        b = [a.conj(i) for i in range(a.n1)]
        return np.prod(b, axis=0).c[0]

    # return content and primitive part of a
    def primitive(a):
        d = 0
        for c in a.c: d = gcd(c,d)
        return d, CYC(a.c//d)

    def __add__(a,b): # a+b
        if isinstance(b,CYC): return CYC(a.c + b.c)
        else: return CYC(a.c + b)

    def __sub__(a,b): # a-b
        if isinstance(b,CYC): return CYC(a.c - b.c)
        else: return CYC(a.c - b)

    def __mul__(a,b): # a*b
        # b is either CYC or int
        if not isinstance(b,CYC):
            return CYC(b * a.c)
        c = np.convolve(a.c, b.c)
        c[:a.n1] += c[a.n:]
        return CYC(c[:a.n] - c[a.n1]) # normalize

    def __truediv__(a,b):
        """ test if b divides a.
        return None if b doesn't divide a,
        return b/a if b divides a.
        """
        if isinstance(b,CYC):
            c = [b.conj(i) for i in range(1,b.n1)]
            c = np.prod(c, axis=0)
            a *= c; b = (b*c).c[0]
        if np.any(a.c%b): return None
        else: return CYC(a.c//b)

    def to_int(a): # const term afeter normalizing
        return a.c[0] - a.c[-1]

    def isRational(a):
        """ test if a is rational
        return None if a is not rational
        return a.to_int() if a is raitonal
        """
        if np.all(a.c[1:] == a.c[1]):
            return a.to_int()
        else: return None

    def __eq__(a,b): # equality testing
        # b is either CYC or int
        if isinstance(b,CYC):
            return np.all(a.c - a.c[-1] ==
                          b.c - b.c[-1])
        else:
            a = a.isRational()
            if not a: return False
            else: return a == b

    def __pow__(a,e):
        """ a^e, assume e>=0 """
        if e==0: return CYC(1)
        m = (1<<(e.bit_length()-1))>>1
        b = a
        while m:
            b *= b
            if e&m: b *= a
            m>>=1
        return b

    def isUnit(a): # test if a is unit
        return a.norm() == 1

    def isAssoc(a,b):
        """ test if a and b are associates """
        q = a/b
        if q: return q.isUnit()
        else: return False

    def __radd__(a,b): return CYC(b + a.c)
    def __rsub__(a,b): return CYC(b - a.c)
    def __rmul__(a,b): return CYC(b * a.c)

    @staticmethod
    def random(b=1):
        """ random cyclotomic integer
        with |coefficients| < b
        """
        return CYC(np.random.randint(-b+1, b, CYC.n1))

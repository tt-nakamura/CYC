# reference:
# H.M.Edwards "Fermat's Last Theorem" sections 4.4, 4.7

from Period import *
from LLL import LLL
from linsolve import linsolve,InvMod

def SolveCongruence(u0,p):
    """ given u0 == eta_0 (mod p)
    get u_j == eta_j (mod p) (j=1,...,e-1)
    and return u_j (j=0,...,e-1)
    """
    u = -u0*Period.w[0,1:,0]
    A = Period.w[0,1:,1:].copy()
    for i in range(len(A)): A[i,i] -= u0
    return linsolve(A,u,p)

def FactorPrime(p, f=None):
    """ factoring prime number p
    p: prime number (int)
    f: order of p mod n
    return a factor pi of p such that norm(pi)==p^f.
    raise RuntimeError if factor couldn't be found.
    Assume CYC.init(n) has been executed.
    Assume n<=19 is prime.
    If f is none, f is set to ord_n(p).
    """
    if p == CYC.n: return CYC([1,-1])
    if f is None: f = CYC.ord(p)
    if f == CYC.n1: return CYC(p)
    Period.init(f)
    F = sp.factor_list(Period.MPoly, modulus=p)
    for q in F[1]:
        u0 = int(-q[0].coeff(x,0))%p
        try: u = SolveCongruence(u0,p)
        except: continue
        break

    u = (p-u) * InvMod(u0,p) % p
    a,e1 = p, Period.e-1
    for b in LLL(np.c_[u, np.eye(e1)]):
        y = Period(b)
        N = abs(y.norm())
        if N==p: return y.to_CYC()
        g = sp.factorint(N//p)
        m = max(g.keys())
        if m<a: a,h,s = m,g,y

    if a==p:
        raise RuntimeError("factor not found")
    q = s.to_CYC()
    for k in h:
        j = CYC.ord(k)
        s = FactorPrime(k,j)
        m = h[k]*f
        while m:
            t = q/s # trial division
            if t: m,q = m-j,t
            else: s = s.conj()
    return q

def factor(a):
    """ factoring cyclotomic integer a
    a: CYC object or ordinary integer
    return dictionary of {factor: exponent} pair
       such that product of factor^exponent
       is associate of a
    Assume CYC.init(n) has been executed.
    Assume n<=19 is prime.
    """
    if CYC.n>=23:
        raise RuntimeError("l>=23 in factor")
    if not isinstance(a,CYC): a = CYC(a)
    F = {}
    if a==0: return F
    d,a = a.primitive()
    G = sp.factorint(a.norm()) # primitive part
    H = sp.factorint(d) # and content
    n,n1 = CYC.n, CYC.n1

    k = G.pop(n,0) + n1*H.pop(n,0)
    if k: F[CYC([1,-1])] = k # ramify
    for k in H:
        f = CYC.ord(k)
        p = FactorPrime(k,f)
        e = n1//f
        for _ in range(e):
            F[p] = H[k]
            p = p.conj()
        if k in G:
            while G[k]:
                b = a/p # trial division
                if b:
                    F[p] += 1
                    G[k] -= f
                    a = b
                else: p = p.conj()
            G.pop(k)

    for k in G:
        f = CYC.ord(k)
        p = FactorPrime(k,f)
        while G[k]:
            b = a/p # trial division
            if b:
                if p in F: F[p] += 1
                else: F[p] = 1
                G[k] -= f
                a = b
            else: p = p.conj()

    return F

def GenPrime(l, f=None, NTRY=1000):
    """ generate random cyclotomic prime integer pi
    l: bit length of prime number p below pi
    f: order of p mod n
    NRY: number of trials to generate p of order f
    return random pi such that norm(pi)=p^f (except p=n).
    Assume CYC.init(n) has been executed.
    Assume n<=19 is prime.
    If f is None, pi of any f is returned.
    """
    if f is None:
        p = sp.randprime(1<<(l-1), 1<<l)
    elif CYC.n1%f:
        raise RuntimeError("f must divide n-1")
    else:
        for _ in range(NTRY):
            p = sp.randprime(1<<(l-1), 1<<l)
            if CYC.ord(p) == f: break
        else:
            raise RuntimeError("NTRY failed")
    return FactorPrime(p,f)

import numpy as np

def LLL(A):
    """ Lattice reduction by Lenstra-Lenstra-Lovasz
    A: row vectors (shape(m,n)).
    return reduced vectors (shape(m,n)).
    raise exception if row vectors are dependent.
    Assume m<=n.
    Floating point arithmetic (long double) is used.
    referece: H.Cohen
      "A Course in Computational Algebraic Number Thoery"
       Algorithm 2.6.3
    """
    B = np.array(A, dtype=np.longdouble)
    m,k,kmax = B.shape[0],1,0
    C = np.empty_like(B)
    c = np.empty(m, dtype=np.longdouble)
    M = np.eye(m, dtype=np.longdouble)
    C[0],c[0] = B[0], np.dot(B[0], B[0])
    while k < m:
        if k > kmax:
            kmax = k
            M[k,:k] = np.dot(C[:k], B[k])/c[:k]
            C[k] = B[k] - np.dot(M[k,:k], C[:k])
            c[k] = np.dot(C[k], C[k])
            if np.isclose(c[k],0):
                raise RuntimeError("singular matrix in LLL")

        reduce(k,k-1,B,M)
        u = M[k,k-1]
        d = c[k] + u**2*c[k-1]
        if d >= 0.75*c[k-1]:
            for l in range(2,k+1): reduce(k,k-l,B,M)
            k += 1
        else:
            B[[k-1,k]] = B[[k,k-1]]
            M[[k-1,k],:k-1] = M[[k,k-1],:k-1]
            M[k,k-1],c[k] = u*c[k-1]/d, c[k]/d
            C[k-1],C[k] = C[k] + u*C[k-1],\
                c[k]*C[k-1] - M[k,k-1]*C[k]
            c[k],c[k-1] = c[k-1]*c[k], d
            M[k+1:,[k-1,k]] = M[k+1:,[k,k-1]]
            M[k+1:,k] -= u*M[k+1:,k-1]
            M[k+1:,k-1] += M[k,k-1]*M[k+1:,k]
            if k>1: k -= 1
    B = [[int(b) for b in r] for r in B]
    return np.asarray(B, dtype='object')

def reduce(k,l,B,M):
    q = np.round(M[k,l])
    B[k] -= q*B[l]
    M[k,:l+1] -= q*M[l,:l+1]

import numpy as np

def linsolve(A,b,p):
    """ solve linear system Ax == b (mod p), where
    A: square matrix (shape(m,m)) of int,
    b: RHS vector (shape(m,...)) of int.
    return x (shape(m,...)) mod p.
    raise exception if A is singular.
    """
    B = np.asarray(np.c_[A,b], dtype='object')
    m = B.shape[0]
    for k in range(m):
        for j in range(k,m):
            if B[j,k]%p: break
        else:
            raise RuntimeError("singular matrix in linsolve")
        if j>k: B[[k,j],k:] = B[[j,k],k:]
        B[k,k:] *= InvMod(B[k,k], p)
        B[k,k:] %= p
        B[k+1:,k:] -= B[k+1:,k,np.newaxis] * B[k,k:] % p
        B[:k,k:] -= B[:k,k,np.newaxis] * B[k,k:] % p
    return np.squeeze(B[:,m:]%p)

def InvMod(a,p):
    """ return x such that ax==1 (mod p) """
    s,u = 1,0
    while p:
        q,r = divmod(a,p)
        a,p,s,u = p,r,u,s-q*u
    return s

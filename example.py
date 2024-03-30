from CYCFactoring import CYC, GenPrime, factor

CYC.init(5)
for _ in range(10):
    a = CYC.random(5)
    b = CYC.random(5)
    c = a*b
    f = factor(c); print(c,f)
    d = 1
    for k in f: d *= k**f[k]
    if not d.isAssoc(c):
        raise RuntimeError("wrong")

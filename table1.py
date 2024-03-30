from CYC import CYC
from CYCFactoring import GenPrime

for n in [5,7,11,13,17,19]:
    CYC.init(n)
    p = GenPrime(20,1)
    print(n, p.norm(), p)

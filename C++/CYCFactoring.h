// uses NTL
//   http://www.shoup.net/ntl

#ifndef __CYCFactroing_h__
#define __CYCFactroing_h__

#include "CYC.h"
#include<NTL/Pair.h>

void FactorPrime(CYC& q, const NTL::ZZ& p, long f=0);
// q = cyclotomic prime integer that divides prime number p
// f = order of p mod CYC::n
// Assume CYC::init(n) has been executed.
// If f<=0, f is set to CYC::ord(p)

void factor(NTL::Vec<NTL::Pair<CYC, long> >& F, const CYC& a);
// F = factorization of a into cycrotomic primes
// each element of F is a pair of prime and its exponent
// such that product of prime^{exponent} is associate of a
// factors of norm l are inserted first in F (if any)
// other factors are sorted by norm in increasing order
// Assume CYC::init(n) has been executed.
// Assume n<=19 is prime.

void mul(CYC& a, const NTL::Vec<NTL::Pair<CYC, long> >& q);
// a = product of (cyclotomic prime integer)^{exponent} in q
// each element of q is a pair of integer and exponent

void GenPrime(CYC&pi, long l, long f=0, long NTRY=1000);
// generate random cyclotomic prime integer pi
// l: bit length of prime number p below pi
// f: order of p mod n
// NRY: number of trials to generate p of order f
// return random pi such that norm(pi)==p^f
// If f==0, pi of any f is output.

#endif // __CYCFactroing_h__
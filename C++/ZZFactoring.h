#ifndef __ZZFactoring_h__
#define __ZZFactoring_h__

#include<NTL/vec_ZZ.h>
#include<NTL/pair.h>

void factor(NTL::Vec<NTL::Pair<NTL::ZZ, long> >& f, const NTL::ZZ& n);
// n = integer
// f = prime factorization of |n|
//     vector of (prime, exponent) pair
//     in increasing order of primes

void factor(NTL::Vec<NTL::Pair<long, long> >& f, long n);
// n = integer, |n| < 2^{60}
// f = prime factorization of |n| by trial division
//     vector of (prime, exponent) pair
//     in increasing order of primes

void divisor(NTL::vec_ZZ& d, const NTL::ZZ& n);
// d = vector of positive divisors of n
// d[0] = 1 and d[i] increases

void conductor(NTL::ZZ& f, NTL::ZZ& d, const NTL::ZZ& D);
// D = discriminant, D==0 or 1 (mod 4)
// return f,d such that f**2 divide D,
//   d = D/f**2 == 0 or 1 (mod 4) is fundamental disc.
// if d==1 (mod 4), d is square-free
// if d==0 (mod 4), d/4 is square-free and
//                  d/4 == 2 or 3 (mod 4)

long IsFundDisc(const NTL::ZZ& D);
// test if D is fundamental discriminant

#endif // __ZZFactoring_h__
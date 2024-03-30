// uses NTL
//   http://www.shoup.net/ntl

#ifndef __period_h__
#define __period_h__

#include<NTL/mat_ZZ.h>
#include<NTL/vec_ZZ_p.h>
#include "CYC.h"

struct Period : NTL::vec_ZZ
// cyclotomic periods
// represented by array of ZZ of length e
{
    static long f;// length of Period
    static long e;// number of Periods
    static NTL::Mat<long> W;// multiplication table
    // eta_0*eta_i = sum_{j=0}^{e-1} W[i][j]*eta_j
    static NTL::ZZX MPoly;// minimal polynomial
    static void init(long f);// set length of period
    static long MulTab(int i, int j, int k)// multiplication table
    // eta_i*eta_j = sum_{k=0}^{l-1} MulTab(i,j,k)*eta_k
    { return W[j<i ? j-i+e:j-i][k<i ? k-i+e:k-i]; }
    Period() { SetLength(e); }
};

inline void set(Period& a)// a = 1
{ for(long i=0; i<a.e; i++) a[i] = -1; }

void mul(Period& c, const Period& a, const Period& b);// c = a*b

Period& operator*=(Period& b, const Period& a);// b *= a

long IsRational(NTL::ZZ& x, const Period& a);
// test if a is rational
// return 1 and set x=-a[0] if a is rational
// return 0 if a is not rational (x is unchanged)

inline void conv(NTL::ZZ& x, const Period& a)
// set x=-a[0] without testing if a is rational
{ negate(x, a[0]); }

void conv(CYC& b, const Period& a);
// transform a to cyclotomic integer b
// as a polynomial in omega = exp(2pi*i/n)

#endif // __period_h__
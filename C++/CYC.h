// uses NTL
//   http://www.shoup.net/ntl

#ifndef __CYC_h__
#define __CYC_h__

#include<NTL/ZZX.h>
#include<NTL/matrix.h>

struct CYC : NTL::ZZX // cyclotomic integers
// as polynomial in omega = exp(2*pi*i/n)
// represented by array of ZZ of length n
// in order of increasing power of omega
{
    static long n; // degree + 1
    static long g; // generator (primitive root)
    static NTL::Vec<long> g_pow; // g^j % n (j=0...n-1)
    static NTL::Vec<long> log_g; // ind_g a (a=1...n-1)
    static NTL::Vec<long> order; // order of a mod n (a=1...n-1)
    static NTL::Mat<long> cjind; // j*g^i%n (i=0...n-1,j=0...n)
    static NTL::ZZX modulus;// cyclotomic polynomial
    static void init(long n_);// set degree
    static long ord(const NTL::ZZ& a)
    { return order[a%n]; }// order of a mod n; ord(0)==1
    CYC() {;}
    CYC(long a) : NTL::ZZX(a) {;} // (*this)[0] = a
    CYC& operator=(const NTL::ZZ& a) { this->ZZX::operator=(a); }
    CYC& operator=(long a) { this->ZZX::operator=(a); }
};

long IsUnit(const CYC& a);
// test if norm(a)==1

long IsAssoc(const CYC& a, const CYC& b);
// test if a/b is unit

void SetCoeff(CYC& a, long i, const NTL::ZZ& c);// a[i] = c
void SetCoeff(CYC& a, long i, long c);// a[i] = c
void SetCoeff(CYC& a, long i);// a[i] = 1

void conv(CYC& a, const NTL::vec_ZZ& c);
// a[i]=c[i] (i=0...,n-1)

void conj(CYC& b, const CYC& a, long i=1);
// b[j*g^i%n] = a[j] (j=0...n-1)

inline void mul(CYC& c, const CYC& a, const CYC& b)
{ NTL::MulMod(c, a, b, c.modulus); } // c = a*b

inline void sqr(CYC& b, const CYC& a)
{ NTL::SqrMod(b, a, b.modulus); } // b = a*a

CYC& operator*=(CYC& b, const CYC& a); // b *= a

inline void norm(NTL::ZZ& n, const CYC& a)
{ NTL::NormMod(n, a, a.modulus); } // n = norm(a)

long divide(const CYC& a, const CYC& b);
// test if b divides a and
// n = norm(b) is optionally supplied.

long divide(CYC& q, const CYC& a, const CYC& b);
// test if b divides a and
// set q = a/b if b divides a, else q is unchanged.
// n = norm(b) is optionally supplied.

void power(CYC& b, const CYC& a, long n);
// b = a^n; assume n>=0

void RandomBnd(CYC& x, const NTL::ZZ& n);
inline void RandomBnd(CYC& x, long n) { RandomBnd(x, NTL::ZZ(n)); }
// x = random cyclotomic integer
//     with |coefficients| less than n

#endif // __CYC_h__
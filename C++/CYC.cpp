// reference:
// H.M.Edwards "Fermat's Last Theorem" section 4.2

// uses NTL
//   http://www.shoup.net/ntl

#include "CYC.h"
#include "ZZFactoring.h"
#include<exception>
using namespace NTL;

long CYC::n; // degree + 1
long CYC::g; // primitive root mod n
ZZX CYC::modulus; // cyclotomic polynomial
Vec<long> CYC::g_pow; // g^j % n (j=0...n-1)
Vec<long> CYC::log_g; // ind_g a (a=1...n-1)
Vec<long> CYC::order; // order of a mod n (a=1...n-1)
Mat<long> CYC::cjind; // j*g^i%n (i=0...n-1,j=0...n)

long primitive_root(long);

void CYC::init(long n_)
// set degree n of cyclotomic integer
// assume n is odd prime
{
    if(n_<3 || !ProbPrime(n_))
        throw std::runtime_error("wrong n in CYC::init");
    n = n_;
    long i,j, n1(n-1), a(1);
    for(i=0; i<n; i++)
        SetCoeff(modulus, i);

    g = primitive_root(n);
    g_pow.SetLength(n1);
    log_g.SetLength(n);
    cjind.SetDims(n1,n);
    order.SetLength(n); order[0] = 1;
    for(i=0; i<n1; i++) {
        g_pow[i] = a;// g^j mod n
        log_g[a] = i;// ind_g a
        order[a] = n1/GCD(i,n1);// order of a mod n
        cjind[i][0] = 0;
        for(j=1; j<n; j++)// indices of conjugates
            cjind[i][j] = AddMod(cjind[i][j-1], a, n);
        a = MulMod(a,g,n);
    }
}

long IsUnit(const CYC& a)
// test if norm(a)==1
{ ZZ n; norm(n,a); return IsOne(n); }

long IsAssoc(const CYC& a, const CYC& b)
// test if a/b is unit
{ CYC q; return divide(q,a,b) && IsUnit(q); }

void SetCoeff(CYC& a, long i, const ZZ& c)
{// a[i] = c
    SetCoeff((ZZX&)a, i%=a.n, c);
    if(i==a.n-1) a %= a.modulus;
}

void SetCoeff(CYC& a, long i, long c)
{// a[i] = c
    SetCoeff((ZZX&)a, i%=a.n, c);
    if(i==a.n-1) a %= a.modulus;
}

void SetCoeff(CYC& a, long i)
{// a[i] = 1
    SetCoeff((ZZX&)a, i%=a.n);
    if(i==a.n-1) a %= a.modulus;
}

void conv(CYC& a, const NTL::vec_ZZ& c)
{// a[i] = c[i] (i=0...n-1)
    conv((ZZX&)a, c);
    a %= a.modulus;
}

CYC& operator*=(CYC& b, const CYC& a)// b *= a
{ MulMod(b, b, a, b.modulus); return b; }

void conj(CYC& b, const CYC& a, long i)
// b = i-th conjugate of a
// replacing omega by omega^{g^i}
{
    if(&b==&a) { CYC c(a); conj(b,c); return; }
    long j;
    b.SetLength(a.n);
    for(j=0; j<=deg(a); j++)
        b[b.cjind[i][j]] = a[j];
    for(j=deg(a)+1; j<a.n; j++)
        clear(b[b.cjind[i][j]]);
    b %= b.modulus;
}

long divide(const CYC& a, const CYC& b)
// test if b divides a and
// n = norm(b) is optionally supplied.
{
    long i;
    CYC c,d,s(b),t(1);
    for(i=2; i<a.n; i++) { conj(s,s); t *= s; }
    mul(c,a,t);
    mul(d,b,t);
    for(i=0; i<=deg(c); i++)
        if(!divide(c[i], d[0])) return 0;
    return 1;
}

long divide(CYC& q, const CYC& a, const CYC& b)
// test if b divides a and
// set q = a/b if b divides a, else q is unchanged.
// n = norm(b) is optionally supplied.
{
    long i;
    ZZ u[q.n];
    CYC c,d,s(b),t(1);
    for(i=2; i<q.n; i++) { conj(s,s); t *= s; }
    mul(c,a,t);
    mul(d,b,t);
    for(i=0; i<=deg(c); i++)
        if(!divide(u[i], c[i], d[0])) return 0;
    clear(q);
    for(i=0; i<=deg(c); i++) SetCoeff(q, i, u[i]);
    return 1;
}

void power(CYC& b, const CYC& a, long n)
// b = a^n; assume n>=0
{
    if(n==0) { set(b); return; }
    if(&b==&a) { CYC c(a); power(b,c,n); return; }
    long m(1<<(NumBits(n)-1));
    b=a;
    for(m>>=1; m; m>>=1) {
        sqr(b,b);
        if(n&m) b*=a;
    }
}

void RandomBnd(CYC& x, const ZZ& n)
// x = random cyclotomic integer
//     with |coefficients| less than n
{
    x.SetLength(x.n-1);
    for(long i=x.n-2; i>=0; i--) {
        RandomBnd(x[i], n);
        if(RandomBits_long(1))
            negate(x[i], x[i]);
    }
    x.normalize();
}
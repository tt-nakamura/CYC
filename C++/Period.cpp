// reference:
// H.M.Edwards "Fermat's Last Theorem" section 4.5

// uses NTL
//   http://www.shoup.net/ntl

#include "Period.h"
#include<exception>
using namespace NTL;

long Period::f;// length of period
long Period::e;// number of conjugates
Mat<long> Period::W;// multiplication table
ZZX Period::MPoly;// minimal polynomial

void Period::init(long f_)
// set length f of period.
// Assume CYC::init(n) has been executed.
// Assume f divides n-1.
{
    long i,j,k,n(CYC::n),n1(n-1);
    if(n1%f_)
        throw std::runtime_error("wrong f in Period");
    f = f_;
    e = n1/f;

    // minimal polynomial
    ZZX eta;
    for(i=0; i<n1; i+=e)
        SetCoeff(eta, CYC::g_pow[i]);
    eta %= CYC::modulus;
    MinPolyMod(MPoly, eta, CYC::modulus);

    // multiplication table
    // eta_0*eta_i = sum_{j=0}^{e-1} W[0][j]*eta_j
    W.SetDims(e,e);
    for(i=0; i<e; i++) {
        for(j=0; j<e; j++) W[i][j] = 0;
        for(j=i; j<n1; j+=e) {
            k = (1 + CYC::g_pow[j])%n;
            if(k) W[i][CYC::log_g[k]%e]++;
            else for(k=0; k<e; k++) W[i][k] -= f;
        }
    }
}

void conv(CYC& b, const Period& a)
// transform a to cyclotomic integer b
// as a polynomial in omega = exp(2pi*i/n)
{
    long i,j,n1(CYC::n-1);
    b.SetLength(CYC::n);
    clear(b[0]);
    for(i=0; i<a.e; i++)
        for(j=i; j<n1; j+=a.e)
            b[CYC::g_pow[j]] = a[i];
    b %= CYC::modulus;
}

void mul(Period& c, const Period& a, const Period& b)
// c = a*b
{
    if(&c==&a) { Period d(a); mul(c,d,b); return; }
    if(&c==&b) { Period d(b); mul(c,a,d); return; }
    long i,j,k;
    ZZ u;
    clear(c);
    for(i=0; i<a.e; i++) {
        for(j=0; j<b.e; j++) {
            mul(u, a[i], b[j]);
            for(k=0; k<c.e; k++)
                MulAddTo(c[k], u, c.MulTab(i,j,k));
        }
    }
}

Period& operator*=(Period& b, const Period& a)
{ mul(b,b,a); return b; } // b *= a

long IsRational(ZZ& x, const Period& a)
// test if a is rational
// return 1 and set x=a[0]-a[1] if a is rational
// return 0 if a is not rational (x is unchanged)
// 
{
    for(long i=0; i<a.e; i++)
        if(a[i] != a[0]) return 0;
    negate(x, a[0]); return 1;
}
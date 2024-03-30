// uses NTL
//   http://www.shoup.net/ntl

#include "ZZFactoring.h"
#include<exception>
using namespace NTL;

#define TRYDIV_BOUND (1<<16)
#define MR_NUM_TRIAL 20
#define RHO_TIME_OUT 5

long IsPrimePower(ZZ& p, const ZZ& n, long N)
// input:
//   n = odd integer, n>=3
//   N = number of times to perform Miller-Rabin test
// output:
//   p = prime factor of n if n = p^k (k>=1)
//       with error probability < 4^{-N}
// return:
//   k if n = p^k (k>=1)
//   0 otherwise
{
    long i;
    ZZ a,d;
    for(p=n;; p=d) {
        for(i=0; i<N; i++) {
            RandomBnd(a,p);
            if(MillerWitness(p,a)) break;
        }
        if(i==N) break;
        PowerMod(d,a,p,p);
        d -= a;
        GCD(d,d,p);
        if(IsOne(d) || d==p) return 0;
    }
    d=n;
    for(i=0; divide(d,d,p); i++);
    return (IsOne(d) ? i:0);
}

long brent_rho(ZZ&, const ZZ&, double);
long mpqs(ZZ&, const ZZ&);

void factor_(Vec<Pair<ZZ, long> >& f, const ZZ& n)
// input:
//   n = odd, integer, n>=3
// output:
//   f = prime factorization of n (appended to f)
{
    long i,j,k(f.length());
    ZZ p,q;
    if(j = IsPrimePower(p, n, MR_NUM_TRIAL)) {
        f.SetLength(k+1);
        f[k].a = p;
        f[k].b = j;
        return;
    }
    Vec<Pair<ZZ, long> > g,h;
    if(brent_rho(p, n, RHO_TIME_OUT) == 0);
    else if(mpqs(p,n) == 0);
    else throw std::runtime_error("factor not found");
    div(q,n,p);
    factor_(g,p);
    factor_(h,q);
    for(i=j=0; i<g.length() || j<h.length(); k++) {
        f.SetLength(k+1);
        if(j==h.length() || i<g.length() && g[i].a < h[j].a)
            f[k] = g[i++];
        else if(i==g.length() || g[i].a > h[j].a)
            f[k] = h[j++];
        else {
            f[k] = g[i++];
            f[k].b += h[j++].b;
        }
    }
}

void factor(Vec<Pair<ZZ, long> >& f, const ZZ& n)
// input:
//   n = integer
// output:
//   f = prime factorization of |n|
//       vector of (prime, exponent) pair
//       in increasing order of primes
{
    long i(0),j,p;
    ZZ m;
    abs(m,n);
    f.SetLength(0);
    if(IsZero(m) || IsOne(m)) return;
    if(j = MakeOdd(m)) {
        f.SetLength(1);
        f[0].a = 2;
        f[0].b = j;
        if(IsOne(m)) return;
        i++;
    }
    PrimeSeq ps;
    ps.reset(3);
    while((p = ps.next()) <= TRYDIV_BOUND) {
        for(j=0; divide(m,m,p); j++);
        if(j==0) continue;
        f.SetLength(i+1);
        f[i].a = p;
        f[i].b = j;
        if(IsOne(m)) return;
        i++;
    }
    factor_(f,m);
}

void conductor(ZZ& f, ZZ& d, const ZZ& D)
// D = discriminant, D==0 or 1 (mod 4)
// return f,d such that
// f**2 divide D and d = D/f**2 == 0 or 1 (mod 4)
// if d==1 (mod 4), d is square-free
// if d==0 (mod 4), d/4 is square-free and
//                  d/4 == 2 or 3 (mod 4)
{
    if(IsZero(D)) { clear(f); clear(d); return; }
    if(&d==&D) { conductor(f,d,ZZ(D)); return; }
    Vec<Pair<ZZ, long> > p;
    factor(p,D);
    set(d);
    for(long i=0; i<p.length(); i++)
        if(p[i].b&1) d *= p[i].a;
    if(sign(D) < 0) negate(d,d);
    div(f,D,d);
    SqrRoot(f,f);
    if(d%4 > 1) { d<<=2; f>>=1; }
}

long IsFundDisc(const ZZ& D)
// test if D is fundamental discriminant
{
    if(IsZero(D) || IsOne(D) || D%4 > 1) return 0;
    ZZ f,d;
    conductor(f,d,D);
    return IsOne(f);
}

static int cmp(const void *a, const void *b)
{ return compare(*(const ZZ*)a, *(const ZZ*)b); }

void divisor(vec_ZZ& d, const ZZ& n)
// d = vector of positive divisors of n
// d[0] = 1 and d[i] increases
{
    long i,j,k,l,m;
    ZZ p;
    Vec<Pair<ZZ, long> > f;
    factor(f,n);
    d.SetLength(1);
    set(d[0]);
    for(i=0; i<f.length(); i++) {
        l = d.length();
        d.SetLength(l*(f[i].b + 1));
        set(p);
        for(j=l; j<d.length();) {
            p *= f[i].a;
            for(k=0; k<l; k++)
                mul(d[j++], d[k], p);
        }
    }
    qsort((void *)d.data(), d.length(), sizeof(ZZ), cmp);
}
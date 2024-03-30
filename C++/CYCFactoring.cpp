// reference:
// H.M.Edwards "Fermat's Last Theorem" sections 4.4, 4.7

// uses NTL
//   http://www.shoup.net/ntl

#include "CYCFactoring.h"
#include "ZZFactoring.h"
#include "Period.h"
#include<NTL/ZZ_pXFactoring.h>
#include<NTL/LLL.h>
#include<exception>
using namespace NTL;

long SolCongr(Vec<ZZ_p>& u, const ZZ_p& u0)
// given u0 == eta_0 (mod p)
// get u_j == eta_j (mod p) (j=1,...,e-1)
// and output u_j (j=0,...,e-1)
// return 0 if successful, else return -1
{
    long i,j,e1(Period::e-1);
    ZZ_p t;
    Mat<ZZ_p> A;
    Vec<ZZ_p> v;
    A.SetDims(e1,e1);
    v.SetLength(e1);
    for(i=0; i<e1; i++) {
        for(j=0; j<e1; j++)
            A[j][i] = Period::W[i+1][j+1];
        A[i][i] -= u0;
        mul(v[i], u0, -Period::W[i+1][0]);
    }
    solve(t,v,A,v);
    if(IsZero(t)) return -1;
    u.SetLength(Period::e);
    u[0] = u0;
    for(i=0; i<e1; i++) u[i+1] = v[i];
    return 0;
}

void FactorPrime(CYC& q, const ZZ& p, long f)
// q = cyclotomic prime integer that divides prime number p
// f = order of p mod CYC::n
// Assume CYC::init(n) has been executed
// If f<=0, f is set to CYC::ord(p)
{
    if(p == CYC::n) { q=1; SetCoeff(q,1,-1); return; }
    if(f<=0) f = CYC::ord(p);
    if(f == CYC::n-1) { q=p; return; }

    long i,j,e((CYC::n-1)/f),e1(e-1);
    ZZ_pPush p_(p);
    ZZ pf,N,a;
    ZZ_pX g;
    Mat<ZZ> B;
    Mat<ZZ_p> C;
    Vec<ZZ_p> r;
    Vec<Pair<ZZ,long> > G,H;
    CYC s;
    
    power(pf,p,f);
    Period::init(f);
    conv(g, Period::MPoly);
    FindRoots(r,g);
    C.SetDims(1,e);
    for(i=0; i<e; i++)
        if(SolCongr(C[0], r[i]) == 0) break;

    transpose(C,C);
    kernel(C,C);
    conv(B,C);
    // optional ***************************
    // extra rows to increase probability
    // of getting shorter vector in lattice
    B.SetDims(e1<<1, e);
    for(i=0; i<e1; i++) B[i+e1] = B[i];
    for(i=0; i<e1; i++) B[i+e1][0] -= p;
    //*************************************
    j = LLL(a,B);
    a = pf;
    for(i=B.NumRows()-j; i<B.NumRows(); i++) {
        conv(q, (Period&)B[i]);
        norm(N,q);
        if(N==pf) return;
        factor(G, N/=pf);
        N = G[G.length()-1].a;
        if(N<a) { a=N; H=G; s=q; }
    }
    if(a==pf)
        throw std::runtime_error("factor not found");
    q=s;
    for(i=0; i<H.length(); i++) {
        j = CYC::ord(H[i].a);
        FactorPrime(s, H[i].a, j);
        while(H[i].b)// trial division
            if(divide(q,q,s)) H[i].b -= j;
            else conj(s,s);
    }
}

void factor(Vec<Pair<CYC, long> >& F, const CYC& a)
// F = factorization of a into cycrotomic primes
// each element of F is a pair of prime and its exponent
// such that product of prime^{exponent} is associate of a.
// Assume CYC::init(n) has been executed
// Assume n<=19 is prime.
{
    int i(0),j(0),k(0),n(CYC::n),n1(n-1),e(0),f,g;
    if(n>=23) throw std::runtime_error("n>=23 in factor");
    ZZ s,t;
    CYC b,c;
    Vec<Pair<ZZ, long> > G,H;
    F.SetLength(0);
    if(IsZero(a)) return;
    content(t,a);
    div(b,a,t);
    norm(s,b);
    while(divide(s,s,n)) e++;
    while(divide(t,t,n)) e+=n1;
    if(e) {// ramify
        F.SetLength(k+1);
        FactorPrime(F[k].a, ZZ(n));
        F[k++].b = e;
    }
    factor(G,s);
    factor(H,t);
    while(i<G.length() || j<H.length()) {
        if(j==H.length() || i<G.length() && G[i].a < H[j].a) {
            f = CYC::ord(G[i].a);
            FactorPrime(c, G[i].a, f);
            while(G[i].b) {
                if(divide(b,b,c)) {// trial division
                    F.SetLength(k+1);
                    F[k].a = c;
                    F[k].b = 0;
                    do { F[k].b++; G[i].b -= f; }
                    while(divide(b,b,c));
                    k++;
                }
                conj(c,c);
            }
            i++;
        }
        else {
            f = CYC::ord(H[j].a);
            FactorPrime(c, H[j].a, f);
            e = n1/f;
            F.SetLength(k+e);
            for(g=k; g<F.length(); g++) {
                F[g].a = c;
                F[g].b = H[j].b;
                conj(c,c);
            }
            if(i<G.length() && G[i].a == H[j].a) {
                for(g=k; G[i].b; g++) {// trial division
                    while(divide(b,b,c)) {
                        F[g].b++;
                        G[i].b -= f;
                    }
                    conj(c,c);
                }
                i++;
            }
            j++; k+=e;
        }
    }
}

void mul(CYC& a, const Vec<Pair<CYC, long> >& q)
// a = product of (cycrotomic integer)^{exponent} in q
// each element of q is a pair of integer and exponent
{
    int i;
    CYC b;
    set(a);
    for(i=0; i<q.length(); i++) {
        power(b, q[i].a, q[i].b);
        a *= b;
    }
}

void GenPrime(CYC&pi, long l, long f, long NTRY)
// generate random cyclotomic prime integer pi
// l: bit length of prime number p below pi
// f: order of p mod n
// NRY: number of trials to generate p of order f.
// Output random pi such that norm(pi)=p^f (except p=n).
// Assume CYC::init(n) has been executed.
// Assume n<=19 is prime.
// If f<=0, pi of any f is output.
{
    ZZ p;
    if(f<=0) {
        GenPrime(p,l);
        f = CYC::ord(p);
    }
    else if((CYC::n-1)%f)
        throw std::runtime_error("wrong f in GenPrime");
    else {
        long i;
        for(i=0; i<NTRY; i++) {
            GenPrime(p,l);
            if(CYC::ord(p) == f) break;
        }
        if(i==NTRY)
            throw std::runtime_error("i==NTRY in GenPrime");
    }
    FactorPrime(pi, p, f);
}
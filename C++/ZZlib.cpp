// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/ZZ.h>
#include<NTL/pair.h>
using namespace NTL;

long Jacobi(long a, long b)
// input:
//   a,b = integers, 0 <= a < b, b odd
// return:
//   Jacobi symbol (a/b)
{
    long c,t(0);
    while(a) {
        for(c=0; (a&1)==0; c++) a>>=1;
        if(c&1) {
            c = b&7;
            if(c==3 || c==5) t ^= 1;
        }
        if(a&b&2) t ^= 1;
        c = b%a;
        b = a;
        a = c;
    }
    if(b!=1) return 0;
    else if(t) return -1;
    else return 1;
}

long SqrRootMod(long a, long p)
// input:
//   a = integer, 0 <= a < p
//   p = odd prime
// return:
//   x such that x^2 = a (mod p)
//   by Cipolla method
{
    if(a==0) return 0;
    long d, s((p+1)>>1);
    long b,c0(1),c1(0),x0,x1(1);
    do {
        b = rand()%p;
        d = (b*b-a)%p;
        if(d<0) d+=p;
    } while(Jacobi(d,p) >= 0);
    for(x0=b;;) {
        if(s&1) {
            b = (c0*x0 + c1*x1%p*d)%p;
            c1 = (c0*x1 + c1*x0)%p;
            c0 = b;
        }
        s>>=1;
        if(s==0) return b;
        b = (x0*x0 + x1*x1%p*d)%p;
        x1 = (x0*x1<<1)%p;
        x0 = b;
    }
}

void factor(Vec<Pair<long, long> >& f, long n)
// input:
//   n = integer, |n| < 2^{60}
// output:
//   f = prime factorization of |n| by trial division
//       vector of (prime, exponent) pair
//       in increasing order of primes
{
    long i(0),j,m,p;
    
    if(n<0) n=-n;
    if(n==0 || n==1) {
        f.SetLength(0);
        return;
    }
    m = SqrRoot(n);
    PrimeSeq ps;
    while((p = ps.next()) <= m) {
        if(p==0) Error("too large n");
        for(j=0; n%p==0; j++) n/=p;
        if(j==0) continue;
        f.SetLength(i+1);
        f[i].a = p;
        f[i].b = j;
        if(n==1) return;
        i++;
        m = SqrRoot(n);
    }
    f.SetLength(i+1);
    f[i].a = n;
    f[i].b = 1;
}

long primitive_root(long p)
// least primitive root mod p
// Assume p is odd prime
{
    long i,p1(p-1),g(2);
    Vec<Pair<long,long> > q;
    factor(q,p1);
    for(;; g++) {
        for(i=0; i<q.length(); i++)
            if(PowerMod(g, p1/q[i].a, p) == 1)
                break;
        if(i==q.length()) return g;
    }
}

long ord(long a, long p)
// order of a mod p
// Assume p is odd prime and 0<a<p
{
    long i,j,k,f(p-1);
    if(a==0) return -1;
    if(a==1) return 1;
    if(a==f) return 2;
    Vec<Pair<long,long> > q;
    factor(q,f);
    for(i=0; i<q.length(); i++)
        for(j=0; j<q[i].b; j++)
            if(PowerMod(a, k = f/q[i].a, p) == 1)
                f = k;
            else break;
    return f;
}
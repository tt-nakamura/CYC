#include "CYCFactoring.h"
using namespace NTL;

main() {
    CYC::init(13);
    CYC a,b;
    Vec<Pair<CYC, long> > q;
    for(long i=0; i<10; i++) {
        RandomBnd(a,5);
        RandomBnd(b,5);
        a *= b;
        factor(q,a);
        mul(b,q);
        std::cout << q << std::endl;
        if(!IsAssoc(a,b)) Error("wrong");
    }
}

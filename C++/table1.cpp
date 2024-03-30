#include "CYCFactoring.h"
using namespace NTL;

main() {
    long i;
    long n[] = {5, 7, 11, 13, 17, 19};
    CYC pi;
    ZZ N;
    for(i=0; i<6; i++) {
        CYC::init(n[i]);
        GenPrime(pi, 30, 1);
        norm(N,pi);
        std::cout << n[i] << ' ';
        std::cout << N << ' ';
        std::cout << pi << '\n';
    }
}
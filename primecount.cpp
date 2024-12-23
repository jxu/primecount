#include <vector>
#include <iostream>
#include <cmath>

const long long ALPHA = 2; // tuning parameter
const int C = 4;
const int Q = 2 * 3 * 5 * 7; // product first C primes
long long X, CBRTX, Z, ACBRTX;

std::vector<long long> mu_pmin, PRIMES, PRIME_COUNT;
std::vector<int> PHI_C(Q+1);

int sgn(long long x) 
{
    return (x > 0) - (x < 0);
}

void mu_prime_sieve(void) 
{
    long long i, j;
    long long pc = 0;  
    for (i = 0; i <= Z; ++i) {
        mu_pmin[i] = 1;
        PRIME_COUNT[i] = 0;
    }

    for (j = 2; j <= Z; ++j) {
        if (mu_pmin[j] == 1) {
            for (i = j; i <= Z; i += j) {
                mu_pmin[i] = (mu_pmin[i] == 1) ? -j : -mu_pmin[i];
            }
        }
    }

    for (j = 2; j <= Z; ++j) {
        if (mu_pmin[j] == -j) { // prime
            PRIMES[pc+1] = j; // 1-indexed
            ++pc;
        
            for (i = j*j; i <= Z; i += j*j) 
                mu_pmin[i] = 0;
        }

        PRIME_COUNT[j] = pc;
    }
}

void pre_phi_c(void) 
{
    int i, j; 
    // for the moment, PHI_C stores an indicator sieve
    
    for (i = 1; i <= Q; ++i) {
        PHI_C[i] = 1;
    }

    for (i = 1; i <= C; ++i) {
        int p = PRIMES[i]; // ith prime
        for (j = p; j <= Q; j += p) 
            PHI_C[j] = 0;
    }

    // accumulate
    for (i = 1; i <= Q; ++i) {
        PHI_C[i] += PHI_C[i-1];
        //printf("%d ", PHI_C[i]);
    }

}

long long phi_c(long long y) 
{
    return (y / Q) * PHI_C[Q] + PHI_C[y % Q];
}

long long exact_acbrt(long long x)
{
    long long z = ALPHA * cbrt(x);
    for(++z; z*z*z > ALPHA * ALPHA * ALPHA * X; --z) ;
    return z;
}

long long exact_sqrt(long long x)
{
    long long z = sqrt(x);
    for (++z; z*z > x; --z) ;
    return z;
}

long long primecount(void) 
{
    if (X <= Z) 
        return PRIME_COUNT[X];

    // exact floor of acbrtx
    long long iacbrtx = exact_acbrt(X);
    long long isqrtx = exact_sqrt(X);

    //printf("%lld %lld\n", iacbrtx, isqrtx);

    long long a = PRIME_COUNT[iacbrtx];
    long long a2 = PRIME_COUNT[isqrtx];

    printf("a=%lld a2=%lld\n", a, a2);

    // phi2(x,a)
    long long P2 = (a * (a-1) / 2) - (a2 * (a2-1) / 2);

    for (long long b = a + 1; b <= a2; ++b)
        P2 += PRIME_COUNT[X / PRIMES[b]];

    printf("P2=%lld \n", P2);

    long long S0 = 1225;

    long long S = 165;


    return S0 + S + a - P2 - 1;
}

int main() 
{
    X = 10000; // TODO user input

    CBRTX = std::cbrt(X); // integer approx

    Z = (CBRTX * CBRTX / ALPHA);  // approx

    ACBRTX = (ALPHA * CBRTX); // approx

    mu_pmin.resize(Z+1);
    PRIMES.resize(Z+1);
    PRIME_COUNT.resize(Z+1);
    
    mu_prime_sieve();
    pre_phi_c();

    long long pi = primecount();

    /*
    for (long long i = 0; i <= Z; ++i) {
        printf("%lld\t%lld\t%lld\n", mu_pmin[i], primes[i], PRIME_COUNT[i]);
    } 
    */

    printf("%lld\n", pi);
}

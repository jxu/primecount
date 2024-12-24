#include <vector>
#include <iostream>
#include <cmath>

#include "fenwick_tree.hpp"

using namespace std;

// global variables shared by all functions
const long long ALPHA = 2; // tuning parameter
const int C = 4;
const int Q = 2 * 3 * 5 * 7; // product first C primes
long long X, CBRTX, Z, ACBRTX;

vector<long long> MU_PMIN, PRIMES, PRIME_COUNT;
vector<int> PHI_C(Q+1);

// signum: returns -1, 0, or 1 
int sgn(long long x) 
{
    return (x > 0) - (x < 0);
}

void mu_prime_sieve(void) 
{
    long long i, j;
    long long pc = 0;  
    for (i = 0; i <= Z; ++i) {
        MU_PMIN[i] = 1;
        PRIME_COUNT[i] = 0;
    }

    for (j = 2; j <= Z; ++j) {
        if (MU_PMIN[j] == 1) {
            for (i = j; i <= Z; i += j) {
                MU_PMIN[i] = (MU_PMIN[i] == 1) ? -j : -MU_PMIN[i];
            }
        }
    }

    for (j = 2; j <= Z; ++j) {
        if (MU_PMIN[j] == -j) { // prime
            PRIMES[pc+1] = j; // 1-indexed
            ++pc;
        
            for (i = j*j; i <= Z; i += j*j) 
                MU_PMIN[i] = 0;
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

long long cube(long long n)
{
    return n * n * n;
}

long long primecount(void) 
{
    if (X <= Z) 
        return PRIME_COUNT[X];

    // exact floor of acbrtx
    long long iacbrtx = exact_acbrt(X);
    long long isqrtx = exact_sqrt(X);


    long long a = PRIME_COUNT[iacbrtx];
    long long a2 = PRIME_COUNT[isqrtx];

    printf("a=%lld a2=%lld\n", a, a2);

    // phi2(x,a)
    long long P2 = (a * (a-1) / 2) - (a2 * (a2-1) / 2);

    for (long long b = a + 1; b <= a2; ++b)
        P2 += PRIME_COUNT[X / PRIMES[b]];
        

    long long S0 = 0;

    for (long long n = 1; n <= iacbrtx; ++n) 
    {
        // pmin(1) = +inf
        if (n == 1 || abs(MU_PMIN[n]) > PRIMES[C])
        {
            S0 += sgn(MU_PMIN[n]) * phi_c(X / n);
        }
    }

    // phi(x,a) Special leaves
    // basic S computation
    long long S = 0;

    // sieve out first C primes, storing odd values
    vector<long long> sieve_ind(Z / 2, 1);

    for (long long i = 2; i <= C; ++i) // skip p = 2
    {
        long long p = PRIMES[i];
        for (long long j = p; j < Z; j += 2*p) 
        {
            sieve_ind[j/2] = 0;
        }
    }

    // phi(y,b) = tree sum up to index (y-1)/2 
    fenwick_tree dyn_sieve(sieve_ind);

    for (long long b = C; b < a-1; ++b) 
    {
        //for (auto x : sieve_ind)
        //    std::cout << x << " ";
        //std::cout << std::endl;
        //std::cout << "b " << b << " "; 
        
        long long pb1 = PRIMES[b+1];

        for (long long m = iacbrtx; m > pb1 && m*pb1 > iacbrtx; --m)
        {
            //std::cout << m << " "; 
            if (abs(MU_PMIN[m]) > pb1)
            {
                long long y = X / (m * pb1);
                long long phi_b = dyn_sieve.sum_to((y-1) / 2);

                //std::cout << sgn(MU_PMIN[m]) * phi_b << " "; 

                S -= sgn(MU_PMIN[m]) * phi_b;
            }
        }

        //std::cout << "\n";
        

        // sieve out (odd) multiples of p_(b+1) for next step
        for (long long i = pb1; i < Z; i += 2*pb1)
        {
            if (sieve_ind[i/2]) 
            {
                dyn_sieve.add_to(i/2, -1);
                sieve_ind[i/2] = 0; 
            }
        }
    }

    
    std::cout << S << std::endl;

    return S0 + S + a - 1 - P2;
}

int main() 
{
    X = 1e12; // TODO user input

    CBRTX = std::cbrt(X); // integer approx
    Z = (CBRTX * CBRTX / ALPHA);  // approx
    ACBRTX = (ALPHA * CBRTX); // approx

    std::cout << "X=" << X << "\tcbrtx=" << CBRTX;
    std::cout << "\tZ=" << Z << "\tACBRTX=" << ACBRTX << "\n"; 

    MU_PMIN.resize(Z+1);
    PRIMES.resize(Z+1);
    PRIME_COUNT.resize(Z+1);
    
    mu_prime_sieve();
    pre_phi_c();

    long long pi = primecount();

    std::cout << pi << '\n';
}

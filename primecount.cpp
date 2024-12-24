#include <vector>
#include <iostream>
#include <cmath>

#include "fenwick_tree.hpp"

// global variables shared by all functions
using namespace std; // whatever

const int64_t ALPHA = 2; // tuning parameter
const int C = 7; // precompute phi_c table size (Q)

int64_t Q, X, CBRTX, Z, ACBRTX;

vector<int64_t> MU_PMIN, PRIMES, PRIME_COUNT, PHI_C;

// signum: returns -1, 0, or 1
int sgn(int64_t x)
{
    return (x > 0) - (x < 0);
}

// precompute PRIMES, PRIME_COUNT, and MU_PMIN with a standard sieve
void mu_prime_sieve(void)
{
    int64_t i, j;
    int64_t pc = 0;
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
    Q = 1;
    for (int i = 1; i <= C; ++i)
        Q *= PRIMES[i];

    PHI_C.resize(Q+1);

    int i, j;
    for (i = 1; i <= Q; ++i)
        PHI_C[i] = 1; // indicator


    for (i = 1; i <= C; ++i) {
        int p = PRIMES[i]; // ith prime
        for (j = p; j <= Q; j += p)
            PHI_C[j] = 0;
    }

    // accumulate
    for (i = 1; i <= Q; ++i)
        PHI_C[i] += PHI_C[i-1];


}

// precomputed phi(x,c)
int64_t phi_c(int64_t y)
{
    return (y / Q) * PHI_C[Q] + PHI_C[y % Q];
}

int64_t exact_acbrt(int64_t x)
{
    int64_t z = ALPHA * cbrt(x);
    for(++z; z*z*z > ALPHA * ALPHA * ALPHA * X; --z) ;
    return z;
}

int64_t exact_sqrt(int64_t x)
{
    int64_t z = sqrt(x);
    for (++z; z*z > x; --z) ;
    return z;
}

int64_t comb2(int64_t n)
{
    return n * (n-1) / 2;
}

int64_t sq(int64_t n) 
{
    return n * n;
}

int64_t primecount(void)
{
    if (X <= Z)
        return PRIME_COUNT[X];

    // exact floor of acbrtx
    int64_t iacbrtx = exact_acbrt(X);
    int64_t isqrtx = exact_sqrt(X);


    int64_t a = PRIME_COUNT[iacbrtx];
    int64_t a2 = PRIME_COUNT[isqrtx];

    cout << "a = " << a << "\n";

    // phi2(x,a)
    int64_t P2 = comb2(a) - comb2(a2);

    for (int64_t b = a + 1; b <= a2; ++b) {
        P2 += PRIME_COUNT[X / PRIMES[b]];
    }

    cout << "P2 = " << P2 << "\n";

    int64_t S0 = 0;

    for (int64_t n = 1; n <= iacbrtx; ++n) {
        // pmin(1) = +inf
        if (n == 1 || abs(MU_PMIN[n]) > PRIMES[C]) {
            S0 += sgn(MU_PMIN[n]) * phi_c(X / n);
        }
    }

    cout << "S0 = " << S0 << "\n";

    // phi(x,a) Special leaves

    // sieve out first C primes, only storing odd values (skipping p_1 = 2)

    vector<int64_t> sieve_ind(Z / 2, 1);

    for (int64_t i = 2; i <= C; ++i) {
        int64_t p = PRIMES[i];
        for (int64_t j = p; j < Z; j += 2*p) {
            sieve_ind[j/2] = 0;
        }
    }

    int64_t S = 0;

    // phi(y,b) = tree sum up to index (y-1)/2
    fenwick_tree dyn_sieve(sieve_ind);

    for (int64_t b = C; b < a-1; ++b) {
        int64_t pb1 = PRIMES[b+1];

        // Case 1 leaves: Algorithm 1
        if (pb1*pb1 <= iacbrtx) {
            int64_t S1b = 0;

            for (int64_t m1b = iacbrtx; m1b * pb1 > iacbrtx; --m1b) {
                 if (abs(MU_PMIN[m1b]) > pb1) {
                    int64_t y = X / (m1b * pb1);
                    int64_t phi_b = dyn_sieve.sum_to((y-1) / 2);

                    S1b -= sgn(MU_PMIN[m1b]) * phi_b;
                }               
            }

            S += S1b;
        }

        // Case 2 leaves: Algorithm 2 hell
        else {
            int64_t xpb12 = X / (pb1*pb1);

            // number of trivial leaves is a + 1 - tb
            int64_t tb;

            if (xpb12 <= pb1) {
                tb = b + 2;
            } else if (xpb12 < iacbrtx) {
                tb = PRIME_COUNT[xpb12] + 1;    
            } else {
                tb = a + 1;
            }

            // step 1
            int64_t d2b = tb - 1; // largest d not considered yet
            int64_t S2b = a - d2b;
            char t = 0;

            while (d2b > b + 1) { // step 2
                int64_t y = X / (pb1 * PRIMES[d2b]);

                if (t == 0) { // step 3, clustered easy leaves
                    if (y >= iacbrtx) {
                        t = 2;
                    } else {
                        int64_t l = PRIME_COUNT[y] - b + 1;
                        int64_t d_ = PRIME_COUNT[X / (pb1 * PRIMES[b+l])];

                        // step 4 
                        if ((sq(PRIMES[d_+1]) <= X / pb1) || d_ <= b) {
                            t = 1; // goto step 6
                        } else {
                            S2b += l * (d2b - d_);
                            d2b = d_;
                            continue; // goto 2
                        }
                    }
                }

                if (t == 1) { // step 5, sparse easy leaves
                    if (y >= iacbrtx) {
                        t = 2;
                    } else {
                        int64_t l = PRIME_COUNT[y] - b + 1;
                        // step 6
                        S2b += l;
                        d2b -= 1;
                        continue; // goto 2
                    }
                } 

                if (t == 2) { // steps 7-9, hard leaves
                    S2b += dyn_sieve.sum_to((y-1)/2);
                    d2b -= 1;
                }
            }

            S += S2b;
        }


        // sieve out (odd) multiples of p_(b+1) for next step
        for (int64_t i = pb1; i < Z; i += 2*pb1) {
            if (sieve_ind[i/2]) {
                dyn_sieve.add_to(i/2, -1);
                sieve_ind[i/2] = 0;
            }
        }
    }

    cout << "S = " << S << "\n";

    return S0 + S + a - 1 - P2;
}

int main()
{
    X = 1e12; // TODO user input


    CBRTX = std::cbrt(X); // integer approx
    Z = (CBRTX * CBRTX / ALPHA);  // approx
    ACBRTX = (ALPHA * CBRTX); // approx

    cout << "X = " << X << "\n";
    cout << "acbrtx = " << ACBRTX << "\n";
    cout << "Z = " << Z << "\n";

    // stored for whole interval [1, z]
    // for large values of X, need to subdivide the interval
    MU_PMIN.resize(Z+1);
    PRIMES.resize(Z+1);
    PRIME_COUNT.resize(Z+1);

    mu_prime_sieve();

    pre_phi_c();

    int64_t pi = primecount();

    std::cout << pi << '\n';
}
